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
    
    # Deduplicate per formula
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

### ----------------- COMPARISON ----------------- ###

def compare_fragments(target_fragments, population_fragments, threshold=0.999):
    comparisons = {}
    for formula, targets in target_fragments.items():
        matches = []
        population = population_fragments.get(formula, [])
        for target in targets:
            sims = [(sim.compute_similarity_score(target["fp"], pop["fp"]), pop["sdf"]) for pop in population]
            sims.sort(key=lambda x: -x[0])
            top_3 = sims[:3]
            comparisons[fingerprint_key(target["fp"])] = {
                "target_sdf": target["sdf"],
                "top_matches": top_3,
                "matched": any(s >= threshold for s, _ in top_3)
            }
    return comparisons

### ----------------- MAIN ANALYSIS ----------------- ###

def run_analysis(entry_id, population_file, idx):
    start = time.time()
    print(f"🔍 Processing target {entry_id}...")

    target_fragments = process_target(entry_id)
    n_total = sum(len(v) for v in target_fragments.values())
    print(f"✅ Unique target fragments: {n_total}")

    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f]

    print(f"⚙️ Processing {len(pop_ids)} population entries...")
    population_fragments = process_population(pop_ids)

    print("🔗 Comparing fragments...")
    comparisons = compare_fragments(target_fragments, population_fragments)

    matched_count = sum(1 for v in comparisons.values() if v["matched"])
    print(f"📊 Matched {matched_count}/{n_total} unique fragments.")

    # Save target fragments
    with MoleculeWriter(f"{idx}_{entry_id}_target_unique_fragments.sdf") as writer:
        for data in target_fragments.values():
            for d in data:
                writer.write(Molecule.from_string(d["sdf"], format="sdf"))

    # Save match files
    for i, (fp_key, comp) in enumerate(comparisons.items(), start=1):
        out_name = f"{idx}_{entry_id}_frag{i}_matches.sdf"
        with MoleculeWriter(out_name) as writer:
            writer.write(Molecule.from_string(comp["target_sdf"], format="sdf"))
            for _, sdf in comp["top_matches"]:
                writer.write(Molecule.from_string(sdf, format="sdf"))

    print(f"✅ SDFs written. ⏱️ Elapsed: {time.time() - start:.2f}s")

### ----------------- MAIN ----------------- ###

if __name__ == "__main__":
    start_global = time.time()
    targets = ['ABAHIW']  # Extend as needed
    for i, target in enumerate(targets, start=1):
        run_analysis(target, f"targets/init_populations_protons/{i}_{target}_init_pop.txt", i)
    print(f"Total Time elapsed for {target}: {time.time() - start_global:.2f} seconds")
