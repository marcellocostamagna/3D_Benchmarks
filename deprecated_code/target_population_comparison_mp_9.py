import os
import numpy as np
import multiprocessing as mp
from collections import defaultdict, Counter
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from hsr import fingerprint as fp
from hsr import similarity as sim
import time

os.environ["QT_QPA_PLATFORM"] = "xcb"  # For headless systems

### ----------------- UTILITY FUNCTIONS ----------------- ###

def get_array_from_ccdcmol(ccdcmol):
    return np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ]) - np.mean([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ], axis=0)

def fingerprint_key(fp_obj):
    return tuple(fp_obj if isinstance(fp_obj, list) else fp_obj.tolist())

def component_of_interest(molecule):
    components = molecule.components
    if not components:
        return None
    properties = [{
        "component": comp,
        "is_organometallic": comp.is_organometallic,
        "molecular_weight": sum(atom.atomic_weight for atom in comp.atoms),
        "atom_count": len(comp.atoms)
    } for comp in components]
    heaviest = max(properties, key=lambda x: x["molecular_weight"])
    most_atoms = max(properties, key=lambda x: x["atom_count"])
    for prop in properties:
        criteria_met = sum([
            prop["is_organometallic"],
            prop["component"] == heaviest["component"],
            prop["component"] == most_atoms["component"]
        ])
        if criteria_met >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None

def create_fragment(central_atom):
    fragment = Molecule(identifier=f"{central_atom.label}_frag")
    atom_map = {central_atom: fragment.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = fragment.add_atom(neighbor)
    for bond in central_atom.bonds:
        atom1, atom2 = bond.atoms
        if atom1 in atom_map and atom2 in atom_map:
            try:
                fragment.add_bond(bond.bond_type, atom_map[atom1], atom_map[atom2])
            except Exception:
                pass
    return fragment

def get_fragments(molecule):
    return [create_fragment(atom) for atom in molecule.atoms]

def formula_signature(fragment):
    atoms = fragment.atoms
    central_atom = atoms[0].atomic_symbol  # First atom is the central one
    atom_counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{atom_counts[el]}" for el in sorted(atom_counts))
    return (central_atom, formula)

def generate_fp_data(fragment):
    return {
        "sdf": fragment.to_string("sdf"),
        "fp": fp.generate_fingerprint_from_data(get_array_from_ccdcmol(fragment)),
        "n_atoms": len(fragment.atoms),
        "central_atom": fragment.atoms[0].atomic_symbol,
        "formula": formula_signature(fragment)
    }

### ----------------- TARGET PROCESSING ----------------- ###

def process_target(entry_id, threshold=0.999):
    reader = io.EntryReader("CSD")
    entry = reader.entry(entry_id)
    mol = component_of_interest(entry.molecule)
    fragments = get_fragments(mol)

    grouped = defaultdict(list)
    for frag in fragments:
        data = generate_fp_data(frag)
        grouped[data["formula"]].append(data)

    unique_fragments = defaultdict(list)
    for formula, items in grouped.items():
        for item in items:
            if all(sim.compute_similarity_score(item["fp"], other["fp"]) < threshold for other in unique_fragments[formula]):
                unique_fragments[formula].append(item)

    return unique_fragments

### ----------------- POPULATION PROCESSING ----------------- ###

def process_population_entry(entry_id):
    try:
        reader = io.EntryReader("CSD")
        entry = reader.entry(entry_id)
        mol = component_of_interest(entry.molecule)
        return [generate_fp_data(frag) for frag in get_fragments(mol)]
    except Exception:
        return []

def process_population(entry_ids):
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(process_population_entry, entry_ids)
    flat = [item for sublist in results for item in sublist]
    grouped = defaultdict(list)
    for item in flat:
        grouped[item["formula"]].append(item)
    return grouped

### ----------------- PARALLEL COMPARISON ----------------- ###

def wrapped_compare(args):
    target_data, pop_subset, threshold = args
    sims = []
    matched = False
    for pop in pop_subset:
        score = sim.compute_similarity_score(target_data["fp"], pop["fp"])
        sims.append((score, pop["sdf"]))
        if score >= threshold:
            matched = True
            break  # ✅ Early exit if match found
    sims.sort(key=lambda x: -x[0])
    top_3 = sims[:3]
    return fingerprint_key(target_data["fp"]), {
        "target_sdf": target_data["sdf"],
        "top_matches": top_3,
        "matched": matched
    }

def compare_fragments_parallel(target_fragments, population_fragments, threshold=0.999):
    tasks = []
    fp_to_id = {}
    total_targets = 0

    for formula, target_list in target_fragments.items():
        pop_subset = population_fragments.get(formula, [])
        for target in target_list:
            tasks.append((target, pop_subset, threshold))
            fp_to_id[fingerprint_key(target["fp"])] = target
            total_targets += 1

    results = []
    matched_set = set()

    with mp.Pool(mp.cpu_count()) as pool:
        for fp_key, res in pool.imap_unordered(wrapped_compare, tasks):
            results.append((fp_key, res))
            if res["matched"]:
                matched_set.add(fp_key)
            if len(matched_set) == total_targets:
                pool.terminate()
                break

    return {fp_key: res for fp_key, res in results}

### ----------------- MAIN ANALYSIS ----------------- ###

def run_analysis(entry_id, population_file, idx, threshold=0.999):
    start = time.time()
    print(f"🔍 Processing target {entry_id}...")

    target_fragments = process_target(entry_id, threshold)
    n_total = sum(len(v) for v in target_fragments.values())
    print(f"✅ Unique target fragments: {n_total}")

    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f]

    print(f"⚙️ Processing {len(pop_ids)} population entries...")
    population_fragments = process_population(pop_ids)

    print("🔗 Comparing fragments...")
    threshold = 0.99
    comparisons = compare_fragments_parallel(target_fragments, population_fragments, threshold)

    matched_count = sum(1 for comp in comparisons.values() if comp["matched"])
    print(f"📊 Matched {matched_count}/{n_total} unique fragments.")

    # Save target fragments
    with MoleculeWriter(f"{idx}_{entry_id}_target_unique_fragments.sdf") as writer:
        for data in target_fragments.values():
            for d in data:
                writer.write(Molecule.from_string(d["sdf"], format="sdf"))

    # Save match files
    for i, comp in enumerate(comparisons.values(), start=1):
        out_name = f"{idx}_{entry_id}_frag{i}_matches.sdf"
        with MoleculeWriter(out_name) as writer:
            writer.write(Molecule.from_string(comp["target_sdf"], format="sdf"))
            for _, sdf in comp["top_matches"]:
                writer.write(Molecule.from_string(sdf, format="sdf"))

    print(f"✅ SDFs written. ⏱️ Elapsed: {time.time() - start:.2f}s")

### ----------------- MAIN ----------------- ###

if __name__ == "__main__":
    start_global = time.time()
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
        run_analysis(target, f"targets/init_populations_protons/{i}_{target}_init_pop.txt", i)
    print(f"Total Time elapsed: {time.time() - start_global:.2f} seconds")
