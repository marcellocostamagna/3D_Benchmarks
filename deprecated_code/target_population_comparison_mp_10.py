import os
import time
import numpy as np
import multiprocessing as mp
from collections import defaultdict, Counter
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from hsr import fingerprint as fp
from hsr import similarity as sim

os.environ["QT_QPA_PLATFORM"] = "xcb"  # For headless systems

### ---------- Utility Functions ---------- ###

def get_array_from_ccdcmol(ccdcmol):
    coords = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return coords - coords.mean(axis=0)

def fingerprint_key(fp_obj):
    return tuple(fp_obj.tolist() if isinstance(fp_obj, np.ndarray) else fp_obj)

def formula_signature(fragment):
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, formula)

def generate_fp_data(fragment):
    return {
        "sdf": fragment.to_string("sdf"),
        "fp": fp.generate_fingerprint_from_data(get_array_from_ccdcmol(fragment)),
        "formula": formula_signature(fragment),
        "n_atoms": len(fragment.atoms),
        "central_atom": fragment.atoms[0].atomic_symbol,
    }

# def generate_fp_data(fragment):
#     timings = {}

#     # Start total timer
#     t0 = time.time()

#     # SDF conversion
#     t1 = time.time()
#     sdf = fragment.to_string("sdf")
#     timings["SDF"] = time.time() - t1

#     # Coordinate array generation
#     t2 = time.time()
#     coords = get_array_from_ccdcmol(fragment)
#     timings["Array"] = time.time() - t2

#     # Fingerprint generation
#     t3 = time.time()
#     fingerprint = fp.generate_fingerprint_from_data(coords)
#     timings["Fingerprint"] = time.time() - t3

#     # Formula signature
#     t4 = time.time()
#     formula = formula_signature(fragment)
#     timings["Formula"] = time.time() - t4

#     # Remaining metadata
#     central_atom = fragment.atoms[0].atomic_symbol
#     n_atoms = len(fragment.atoms)
#     timings["Metadata"] = time.time() - t4  # Very quick usually

#     total_time = time.time() - t0
#     timings["Total"] = total_time

#     # Debug output (optional: comment/remove in production)
#     print(f"[{fragment.identifier}] SDF: {timings['SDF']:.2f}s | Array: {timings['Array']:.2f}s | "
#           f"FP: {timings['Fingerprint']:.2f}s | Formula: {timings['Formula']:.2f}s | "
#           f"Total: {total_time:.2f}s")

#     return {
#         "sdf": sdf,
#         "fp": fingerprint,
#         "formula": formula,
#         "n_atoms": n_atoms,
#         "central_atom": central_atom
#     }


### ---------- Fragmentation Functions ---------- ###

def component_of_interest(molecule):
    components = molecule.components
    if not components:
        return None
    props = [{
        "component": c,
        "is_organometallic": c.is_organometallic,
        "mw": sum(atom.atomic_weight for atom in c.atoms),
        "atom_count": len(c.atoms)
    } for c in components]
    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["atom_count"])
    for prop in props:
        if sum([
            prop["is_organometallic"],
            prop["component"] == heaviest["component"],
            prop["component"] == most_atoms["component"]
        ]) >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None

def create_fragment(central_atom):
    frag = Molecule(identifier=f"{central_atom.label}_frag")
    atom_map = {central_atom: frag.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = frag.add_atom(neighbor)
    for bond in central_atom.bonds:
        a1, a2 = bond.atoms
        if a1 in atom_map and a2 in atom_map:
            try:
                frag.add_bond(bond.bond_type, atom_map[a1], atom_map[a2])
            except:
                pass
    return frag

def get_fragments(mol):
    return [create_fragment(atom) for atom in mol.atoms]

### ---------- Target Processing ---------- ###

def process_target(entry_id, threshold=0.999):
    reader = io.EntryReader("CSD")
    mol = component_of_interest(reader.entry(entry_id).molecule)
    fragments = get_fragments(mol)

    grouped = defaultdict(list)
    for frag in fragments:
        data = generate_fp_data(frag)
        grouped[data["formula"]].append(data)

    unique_by_formula = defaultdict(list)
    for formula, frags in grouped.items():
        for frag in frags:
            if all(sim.compute_similarity_score(frag["fp"], other["fp"]) < threshold for other in unique_by_formula[formula]):
                unique_by_formula[formula].append(frag)

    return unique_by_formula

### ---------- Population Processing ---------- ###

def process_population_entry(entry_id):
    try:
        # t0 = time.time()
        reader = io.EntryReader("CSD")
        entry = reader.entry(entry_id)
        # t1 = time.time()

        mol = component_of_interest(entry.molecule)
        # t2 = time.time()

        fragments = get_fragments(mol)
        # t3 = time.time()

        processed = [generate_fp_data(frag) for frag in fragments]
        # t4 = time.time()

        # print(f"[{entry_id}] Load: {t1 - t0:.2f}s | Filter: {t2 - t1:.2f}s | Frag: {t3 - t2:.2f}s | FP: {t4 - t3:.2f}s")
        return processed
    except Exception as e:
        print(f"❌ Error processing {entry_id}: {e}")
        return []

def process_population(entry_ids):
    t0 = time.time()
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(process_population_entry, entry_ids)
    t1 = time.time()
    print(f"✅ Multiprocessing finished in {t1 - t0:.2f}s")

    flat = [item for sublist in results for item in sublist]

    grouped = defaultdict(list)
    for item in flat:
        grouped[item["formula"]].append(item)
    t2 = time.time()
    print(f"✅ Grouping finished in {t2 - t1:.2f}s")

    return grouped

### ---------- Parallel Similarity Comparison ---------- ###

def compare_group(args):
    target_group, pop_group, threshold = args

    # Prepare tracking for each target
    target_status = {
        fingerprint_key(t["fp"]): {
            "data": t,
            "top_matches": [],
            "matched": False
        }
        for t in target_group
    }

    # Track which targets are still unmatched
    unmatched_keys = set(target_status.keys())

    for pop in pop_group:
        pop_fp = pop["fp"]
        pop_sdf = pop["sdf"]

        to_remove = []
        for key in unmatched_keys:
            target = target_status[key]["data"]
            score = sim.compute_similarity_score(target["fp"], pop_fp)

            # Add to top matches
            target_status[key]["top_matches"].append((score, pop_sdf))
            target_status[key]["top_matches"].sort(key=lambda x: -x[0])
            target_status[key]["top_matches"] = target_status[key]["top_matches"][:3]

            # If match found, flag and schedule for removal
            if score >= threshold:
                target_status[key]["matched"] = True
                to_remove.append(key)

        # Remove matched targets from unmatched set
        for key in to_remove:
            unmatched_keys.remove(key)

        # Early exit: all targets matched
        if not unmatched_keys:
            break

    # Final formatting of results
    results = []
    for key, info in target_status.items():
        results.append((
            key,
            {
                "target_sdf": info["data"]["sdf"],
                "top_matches": info["top_matches"],
                "matched": info["matched"]
            }
        ))

    return results

def compare_fragments_parallel(target_fragments, population_fragments, threshold=0.999):
    tasks = []
    for formula, targets in target_fragments.items():
        if formula in population_fragments:
            tasks.append((targets, population_fragments[formula], threshold))

    all_results = {}
    with mp.Pool(mp.cpu_count()) as pool:
        for batch in pool.imap_unordered(compare_group, tasks):
            all_results.update(batch)
    return all_results

### ---------- Main Orchestration ---------- ###

def run_analysis(entry_id, population_file, idx, threshold=0.999):
    start = time.time()
    print(f"🔍 Target: {entry_id}")

    t0 = time.time()
    target_frags = process_target(entry_id, threshold)
    t1 = time.time()
    total = sum(len(v) for v in target_frags.values())
    print(f"✅ Unique target fragments: {total} (in {t1 - t0:.2f}s)")

    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f]

    print(f"⚙️ Population size: {len(pop_ids)}")
    t2 = time.time()
    pop_frags = process_population(pop_ids)
    t3 = time.time()
    print(f"✅ Population processed in {t3 - t2:.2f}s")

    print("🔗 Running similarity comparison...")
    t4 = time.time()
    comparisons = compare_fragments_parallel(target_frags, pop_frags, threshold)
    t5 = time.time()
    print(f"✅ Comparison done in {t5 - t4:.2f}s")

    matched = sum(1 for v in comparisons.values() if v["matched"])
    print(f"📊 Matched {matched}/{total} fragments.")

    t6 = time.time()
    with MoleculeWriter(f"{idx}_{entry_id}_target_unique_fragments.sdf") as w:
        for group in target_frags.values():
            for frag in group:
                w.write(Molecule.from_string(frag["sdf"], format="sdf"))

    for i, (fp_key, comp) in enumerate(comparisons.items(), start=1):
        with MoleculeWriter(f"{idx}_{entry_id}_frag{i}_matches.sdf") as writer:
            writer.write(Molecule.from_string(comp["target_sdf"], format="sdf"))
            for _, sdf in comp["top_matches"]:
                writer.write(Molecule.from_string(sdf, format="sdf"))
    t7 = time.time()

    print(f"📦 Writing SDFs took {t7 - t6:.2f}s")
    print(f"✅ Done in {time.time() - start:.2f}s")


### ---------- Script Entry ---------- ###

if __name__ == "__main__":
    start_all = time.time()
    targets = ['ABAHIW', 'ABAKIZ', 'ABADOX']
    for i, target in enumerate(targets, start=1):
        run_analysis(target, f"targets/init_populations_protons/{i}_{target}_init_pop.txt", i)
    print(f"\n⏱️ Total time: {time.time() - start_all:.2f} seconds")
