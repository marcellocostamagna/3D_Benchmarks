import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import time
from collections import defaultdict
from ccdc import io
from hsr import fingerprint as fp
from hsr import similarity as sim
from scipy.spatial import KDTree

os.environ["QT_QPA_PLATFORM"] = "xcb"  # Qt fix for X11

### ---------- Fragmentation & Fingerprinting ---------- ###

def get_array_from_ccdcmol(ccdcmol):
    return np.array([
        [a.coordinates[0], a.coordinates[1], a.coordinates[2], np.sqrt(a.atomic_number)]
        for a in ccdcmol.atoms
    ]) - np.mean([
        [a.coordinates[0], a.coordinates[1], a.coordinates[2], np.sqrt(a.atomic_number)]
        for a in ccdcmol.atoms
    ], axis=0)

def component_of_interest(mol):
    components = mol.components
    if not components:
        return None
    props = [{
        "component": c,
        "is_organometallic": c.is_organometallic,
        "mw": sum(atom.atomic_weight for atom in c.atoms),
        "count": len(c.atoms)
    } for c in components]
    heaviest = max(props, key=lambda x: x["mw"])
    largest = max(props, key=lambda x: x["count"])
    for p in props:
        if sum([
            p["is_organometallic"],
            p["component"] == heaviest["component"],
            p["component"] == largest["component"]
        ]) >= 2 and p["count"] >= 5:
            return p["component"]
    return None

def create_fragment(atom):
    return [atom] + list(atom.neighbours)

def get_fragments(mol):
    return [create_fragment(a) for a in mol.atoms]

def fragment_to_fp(fragment_atoms):
    coords = np.array([
        [a.coordinates[0], a.coordinates[1], a.coordinates[2], np.sqrt(a.atomic_number)]
        for a in fragment_atoms
    ])
    coords -= np.mean(coords, axis=0)
    return fp.generate_fingerprint_from_data(coords)

def process_molecule(entry_id):
    try:
        reader = io.EntryReader("CSD")
        entry = reader.entry(entry_id)
        mol = component_of_interest(entry.molecule)
        if mol is None:
            return []
        return [
            (fragment_to_fp(frag_atoms), len(frag_atoms))
            for frag_atoms in get_fragments(mol)
        ]
    except Exception as e:
        print(f"❌ Error with {entry_id}: {e}")
        return []

def process_population(entry_ids, workers=mp.cpu_count()):
    with mp.Pool(workers) as pool:
        results = pool.map(process_molecule, entry_ids)
    return [fp for sub in results for fp in sub]

def group_by_atom_count(fp_list):
    groups = defaultdict(list)
    for fingerprint, n_atoms in fp_list:
        groups[n_atoms].append(fingerprint)
    return groups

def find_unique_fps(fps, threshold=0.999):
    unique = []
    for f in fps:
        if all(sim.compute_similarity_score(f, u) < threshold for u in unique):
            unique.append(f)
    return unique

### ---------- Similarity Matching ---------- ###

def match_population_to_target(target_groups, population_groups, threshold=0.999, top_n=3):
    matched = set()
    top_matches = defaultdict(list)
    unique_population_fps = []

    for atom_count, pop_fps in population_groups.items():
        if atom_count not in target_groups:
            continue

        target_fps = target_groups[atom_count]
        pop_unique = find_unique_fps(pop_fps, threshold)
        unique_population_fps.extend(pop_unique)

        if not pop_unique:
            continue

        tree = KDTree(np.array(pop_unique))
        for t_fp in target_fps:
            dists, indices = tree.query(t_fp, k=min(top_n, len(pop_unique)))
            if not isinstance(indices, (list, np.ndarray)):
                indices = [indices]
                dists = [dists]

            for idx in indices:
                candidate = pop_unique[idx]
                score = sim.compute_similarity_score(t_fp, candidate)
                if score >= threshold:
                    matched.add(tuple(t_fp))
                top_matches[tuple(t_fp)].append((score, candidate))

            # Store only top N
            top_matches[tuple(t_fp)] = sorted(top_matches[tuple(t_fp)], key=lambda x: -x[0])[:top_n]

    return unique_population_fps, top_matches, len(matched)

### ---------- Main Execution ---------- ###

if __name__ == "__main__":
    start = time.time()

    target_id = "ABAHIW"
    with open("targets/init_populations_protons/1_ABAHIW_init_pop.txt") as f:
        pop_ids = [line.strip().split()[0] for line in f][:1000]  # Debug subset

    # pop_ids = ['ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE']
    # print(pop_ids)
    
    # Target processing
    target_fps_raw = process_molecule(target_id)
    target_grouped = group_by_atom_count(target_fps_raw)
    for atom_count in target_grouped:
        target_grouped[atom_count] = find_unique_fps(target_grouped[atom_count])

    # Population processing
    population_fps_raw = process_population(pop_ids)
    population_grouped = group_by_atom_count(population_fps_raw)

    # Matching
    unique_pop_fps, top_matches, n_matched = match_population_to_target(
        target_grouped, population_grouped
    )

    print("\n✅ Matching completed.")
    print(f"Unique target fragments: {sum(len(v) for v in target_grouped.values())}")
    print(f"Unique population fragments: {len(unique_pop_fps)}")
    print(f"Matched target fragments: {n_matched}")
    print(f"⏱ Time: {time.time() - start:.2f} seconds")
