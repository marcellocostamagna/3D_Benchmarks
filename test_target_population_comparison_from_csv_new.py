import os
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import defaultdict, Counter
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc import io
from ccdc.entry import Entry
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
    
def is_biatomic(fragment):
    return len(fragment.atoms) == 2

def interatomic_distance(sdf_string):
    mol = Molecule.from_string(sdf_string, format="sdf")
    coords = [atom.coordinates for atom in mol.atoms]
    if len(coords) != 2:
        raise ValueError("Expected a biatomic molecule.")
    return np.linalg.norm(np.array(coords[0]) - np.array(coords[1]))
    
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

### ---------- CSV-Based Fragment Loading ---------- ###

def load_population_fragments_from_csv(entry_ids, csv_path):
    df = pd.read_csv(csv_path)
    df = df[df["entry_id"].isin(entry_ids)].copy()

    df["fp"] = df["fp"].apply(lambda s: np.array([float(x) for x in s.split(',')]))
    df["formula"] = df["formula"].apply(eval)

    grouped = defaultdict(list)
    for _, row in df.iterrows():
        grouped[row["formula"]].append(row.to_dict())

    return grouped

### ---------- Parallel Similarity Comparison ---------- ###

def compare_group(args):
    target_group, pop_group, threshold = args

    target_status = {
        fingerprint_key(t["fp"]): {
            "data": t,
            "top_matches": [],
            "matched": False
        }
        for t in target_group
    }

    unmatched_keys = set(target_status.keys())

    for pop in pop_group:
        pop_fp = pop["fp"]
        pop_sdf = pop["sdf"]
        pop_formula = pop["formula"]
        pop_atoms = pop["n_atoms"]

        to_remove = []

        for key in unmatched_keys:
            t = target_status[key]["data"]
            is_both_biatomic = t["n_atoms"] == 2 and pop_atoms == 2

            if is_both_biatomic:
                # Ensure we only compare same element pairs
                if t["formula"] != pop_formula:
                    continue
                try:
                    t_dist = interatomic_distance(t["sdf"])
                    p_dist = interatomic_distance(pop_sdf)
                    diff = abs(t_dist - p_dist)

                    if diff <= 0.01:
                        score = 1.0 - diff
                        target_status[key]["top_matches"].append((score, pop_sdf))
                        target_status[key]["matched"] = True
                        to_remove.append(key)
                except Exception as e:
                    print(f"Distance calc error: {e}")
                    continue

            else:
                score = sim.compute_similarity_score(t["fp"], pop_fp)
                if score >= threshold:
                    target_status[key]["top_matches"].append((score, pop_sdf))
                    target_status[key]["matched"] = True
                    to_remove.append(key)

            # Always sort and trim to top 3
            target_status[key]["top_matches"].sort(key=lambda x: -x[0])
            target_status[key]["top_matches"] = target_status[key]["top_matches"][:3]

        for key in to_remove:
            unmatched_keys.discard(key)

        if not unmatched_keys:
            break

    # Final output
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

def run_analysis(entry_id, population_file, idx, csv_path="all_fragments_data_0_5.csv", threshold=0.999):
    start = time.time()
    print(f"\n🔍 Target: {entry_id}")
    
    output_dir = "frag_comparison_results_0_5"
    os.makedirs(output_dir, exist_ok=True)

    t0 = time.time()
    target_frags = process_target(entry_id, threshold)
    t1 = time.time()
    total = sum(len(v) for v in target_frags.values())
    print(f"✅ Unique target fragments: {total} (retrieved in {t1 - t0:.2f}s)")

    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f if line.strip()]
    # # for debugging purposes 
    # pop_ids = pop_ids[:100]
    print(f"⚙️ Population size: {len(pop_ids)}")

    t2 = time.time()
    pop_frags = load_population_fragments_from_csv(pop_ids, csv_path)
    t3 = time.time()
    print(f"📦 Loaded population fragments from CSV in {t3 - t2:.2f}s")

    print("🔗 Running similarity comparison...")
    t4 = time.time()
    threshold = 0.99
    comparisons = compare_fragments_parallel(target_frags, pop_frags, threshold)
    t5 = time.time()
    print(f"✅ Comparison done in {t5 - t4:.2f}s")

    matched = sum(1 for v in comparisons.values() if v["matched"])
    print(f"📊 Matched {matched}/{total} fragments.")

    t6 = time.time()
    with MoleculeWriter(os.path.join(output_dir, f"{idx}_{entry_id}_target_unique_fragments.sdf")) as w:
        for group in target_frags.values():
            for frag in group:
                w.write(Molecule.from_string(frag["sdf"], format="sdf"))
    t7 = time.time()
    print(f"🧪 Wrote target fragments SDF in {t7 - t6:.2f}s")

    for i, (fp_key, comp) in enumerate(comparisons.items(), start=1):
        output_path = os.path.join(output_dir, f"{idx}_{entry_id}_frag{i}_matches.sdf")
        with MoleculeWriter(output_path) as writer:

            # Write the target molecule
            target_mol = Molecule.from_string(comp["target_sdf"], format="sdf")
            writer.write(target_mol)

            # Write the top matches with similarity scores as SD tags
            for sim_score, sdf in comp["top_matches"]:
                match_mol = Molecule.from_string(sdf, format="sdf")
                match_entry = Entry.from_molecule(match_mol)

                # Check if biatomic to store DistanceDifference instead of Similarity
                if len(match_mol.atoms) == 2 and len(target_mol.atoms) == 2:
                    try:
                        dist_target = np.linalg.norm(np.array(target_mol.atoms[0].coordinates) - np.array(target_mol.atoms[1].coordinates))
                        dist_match = np.linalg.norm(np.array(match_mol.atoms[0].coordinates) - np.array(match_mol.atoms[1].coordinates))
                        dist_diff = abs(dist_target - dist_match)
                        match_entry.attributes["DistanceDifference"] = f"{dist_diff:.4f}"
                    except Exception as e:
                        print(f"⚠️ Error calculating distance diff: {e}")
                        match_entry.attributes["DistanceDifference"] = "ERROR"
                else:
                    match_entry.attributes["Similarity"] = f"{sim_score:.4f}"

                writer.write_entry(match_entry)


                
    t8 = time.time()
    print(f"📁 Wrote match SDFs in {t8 - t7:.2f}s")

    total_time = time.time() - start
    print(f"✅ Done in {total_time:.2f}s")


### ---------- Script Entry ---------- ###

if __name__ == "__main__":
    overall_start = time.time()
    targets = ['ABAHIW', 'ABAKIZ', 'ABADOX',
                'ABABIP', 'GASQOK', 'ABEKIE',
                'NIWPUE01', 'ABEKIF', 'APUFEX',
                'ABEHAU', 'TITTUO', 'EGEYOG',
                'ABOBUP', 'XIDTOW', 'ACNCOB10',
                'TACXUQ','ACAZFE', 'NIVHEJ',
                'ADUPAS','DAJLAC', 'OFOWIS',
                'CATSUL','HESMUQ01', 'GUDQOL',
                'ABEVAG', 'AKOQOH', 'ADARUT',
                'AFECIA', 'ACOVUL', 'AFIXEV']
    for i, target in enumerate(targets, start=1):
        run_analysis(target, f"targets/init_populations_protons_2/{i}_{target}_init_pop.txt", i)
    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")
