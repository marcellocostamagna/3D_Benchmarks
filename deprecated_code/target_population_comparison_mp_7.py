import os
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import defaultdict
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from hsr import fingerprint as fp
from hsr import similarity as sim
import time

# For headless environments using CSD Python API
os.environ["QT_QPA_PLATFORM"] = "xcb"

### ---------- Utility Functions ---------- ###

def get_array_from_ccdcmol(ccdcmol):
    return np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ]) - np.mean([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ], axis=0)

def fingerprint_key(fp_obj):
    """Ensures consistent tuple representation of fingerprints."""
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

def generate_fp_and_sdf(fragment):
    arr = get_array_from_ccdcmol(fragment)
    return fp.generate_fingerprint_from_data(arr), fragment.to_string("sdf")

### ---------- Target Processing ---------- ###

def process_target(entry_id, threshold=0.9990):
    csd_reader = io.EntryReader("CSD")
    entry = csd_reader.entry(entry_id)
    molecule = component_of_interest(entry.molecule)
    fragments = get_fragments(molecule)
    fps_sdf = [generate_fp_and_sdf(frag) for frag in fragments]
    unique_by_atoms = defaultdict(list)
    for i, (fp1, sdf1) in enumerate(fps_sdf):
        is_unique = True
        n_atoms = len(fragments[i].atoms)
        for fp2, _ in unique_by_atoms[n_atoms]:
            if sim.compute_similarity_score(fp1, fp2) >= threshold:
                is_unique = False
                break
        if is_unique:
            unique_by_atoms[n_atoms].append((fp1, sdf1))
    return unique_by_atoms

### ---------- Population Processing ---------- ###

def process_population_entry(entry_id):
    try:
        csd_reader = io.EntryReader("CSD")
        entry = csd_reader.entry(entry_id)
        molecule = component_of_interest(entry.molecule)
        fragments = get_fragments(molecule)
        results = []
        for frag in fragments:
            n_atoms = len(frag.atoms)
            fp_frag = fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag))
            sdf_str = frag.to_string("sdf")
            results.append((n_atoms, fp_frag, sdf_str))
        return results
    except Exception:
        return []

def process_population(entries):
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(process_population_entry, entries)
    flat = [item for sublist in results for item in sublist]
    grouped = defaultdict(list)
    for n_atoms, fp, sdf in flat:
        grouped[n_atoms].append((fp, sdf))
    return grouped

### ---------- Similarity Comparison ---------- ###

def compare_batch(target_batch, pop_batch, threshold=0.9990):
    matched_fragments = set()
    top_matches = defaultdict(list)

    for fp_target, _ in target_batch:
        key = fingerprint_key(fp_target)
        for fp_pop, sdf_pop in pop_batch:
            sim_val = sim.compute_similarity_score(fp_target, fp_pop)

            if sim_val >= threshold:
                matched_fragments.add(key)

            top_matches[key].append((sim_val, sdf_pop))
            top_matches[key].sort(key=lambda x: -x[0])
            top_matches[key] = top_matches[key][:3]  

    return matched_fragments, top_matches

### ---------- Orchestration ---------- ###

def run_analysis(target_entry_id, population_file, threshold=0.9990):
    start = time.time()
    print("🔍 Processing target...")
    target_fragments_by_atoms = process_target(target_entry_id, threshold)
    num_unique_target_fragments = sum(len(v) for v in target_fragments_by_atoms.values())
    print(f"✅ Unique target fragments: {num_unique_target_fragments}")

    with open(population_file) as f:
        population_entries = [line.split()[0] for line in f]

    print("⚙️  Processing population in parallel...")
    pop_by_atoms = process_population(population_entries)

    matched_count = 0
    full_top_matches = []

    print("🔗 Comparing by atom coount...")
    for n_atoms, target_batch in target_fragments_by_atoms.items():
        pop_batch = pop_by_atoms.get(n_atoms, [])
        if not pop_batch:
            continue
        matched_fps, top_matches = compare_batch(target_batch, pop_batch, threshold)
        matched_count += len(matched_fps)
        for fps in top_matches.values():
            full_top_matches.extend(fps)
        if matched_count == num_unique_target_fragments:
            print("\n✅ All target fragments matched.")
            break

    print(f"\n📊 Matched {matched_count}/{num_unique_target_fragments} target fragments.")

    with MoleculeWriter("matched_population_fragments.sdf") as writer:
        for _, sdf in full_top_matches:
            mol = Molecule.from_string(sdf, format="sdf")
            writer.write(mol)

    print(f"\n✅ Saved matched fragments to 'matched_population_fragments.sdf'")
    print(f"⏱️  Time elapsed: {time.time() - start:.2f} seconds")

### ---------- Main Entry Point ---------- ###

if __name__ == "__main__":
    target_entries = ['ABAHIW'] # 'ABAKIZ', 'ABADOX',
    #                   'ABABIP', 'GASQOK', 'ABEKIE',
    #                   'NIWPUE01', 'ABEKIF', 'APUFEX',
    #                   'ABEHAU', 'TITTUO', 'EGEYOG',
    #                   'ABOBUP', 'XIDTOW', 'ACNCOB10',
    #                   'TACXUQ','ACAZFE', 'NIVHEJ',
    #                   'ADUPAS','DAJLAC', 'OFOWIS',
    #                   'CATSUL','HESMUQ01', 'GUDQOL',
    #                   'ABEVAG', 'AKOQOH', 'ADARUT',
    #                   'AFECIA', 'ACOVUL', 'AFIXEV']

    for i, target in enumerate(target_entries, start=1):
        start_global = time.time()
        print(f"Processing target: {target}")
        run_analysis(target, f"targets/init_populations_protons/{i}_{target}_init_pop.txt")
        print(f"Total Time elapsed for {target}: {time.time() - start_global:.2f} seconds")
    