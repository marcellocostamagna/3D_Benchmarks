import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import time
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from hsr import fingerprint as fp
from hsr import similarity as sim

# Ensure proper X11 configuration for Qt
os.environ["QT_QPA_PLATFORM"] = "xcb"

### ------------------ Molecule Processing Functions ------------------ ###

def get_array_from_ccdcmol(ccdcmol):
    """Extracts an array of atomic coordinates and atomic numbers."""
    return np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ]) - np.mean([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ], axis=0)  # Centering data


def component_of_interest(molecule):
    """Identifies the most relevant component of a molecule."""
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
    """Creates a molecular fragment centered on a given atom, including its first neighbors and bonds."""
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
                pass  # Ignore duplicate bonds
    
    return fragment


def get_fragments(molecule):
    """Generates molecular fragments for all atoms."""
    return [create_fragment(atom) for atom in molecule.atoms]


### ------------------ Efficient Fragmentation & Fingerprinting ------------------ ###

def find_unique_fragments_optimized(fragments=None, fps=None, threshold=0.9990):
    """Finds unique fragments efficiently using a hash-based set."""
    unique_fps = []
    unique_fragments = [] if fragments else None  # Only track fragments if provided
    seen = set()

    for i, fp in enumerate(fps):
        fp_tuple = tuple(fp.tolist()) if isinstance(fp, np.ndarray) else tuple(fp)

        if fp_tuple in seen:
            continue

        is_unique = all(sim.compute_similarity_score(fp, stored_fp) < threshold for stored_fp in unique_fps)

        if is_unique:
            unique_fps.append(fp)
            if unique_fragments is not None:
                unique_fragments.append(fragments[i])  # Only store fragments if provided
            seen.add(fp_tuple)

    return (unique_fragments, unique_fps) if unique_fragments is not None else unique_fps


def process_molecule_parallel(entry_id):
    """Parallel function for processing a molecule: returns both fingerprints and fragments."""
    try:
        csd_reader = io.EntryReader('CSD')
        entry = csd_reader.entry(entry_id)
        molecule = component_of_interest(entry.molecule)
        fragments = get_fragments(molecule)

        # Compute fingerprints
        fragment_fps = [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in fragments]

        return list(zip(fragments, fragment_fps))  # ✅ Return (fragment, fingerprint) pairs
    except Exception as e:
        print(f"❌ Error processing {entry_id}: {e}")
        return []


def parallel_population_processing(entry_ids, num_workers=mp.cpu_count()):
    """Processes all molecules in parallel, returning both fragments and fingerprints."""
    with mp.Pool(num_workers) as pool:
        all_fragments_fps = pool.map(process_molecule_parallel, entry_ids)
    
    # Flatten lists while keeping fragment-fingerprint pairs
    all_fragments, all_fps = [], []
    for sublist in all_fragments_fps:
        for frag, fp in sublist:
            all_fragments.append(frag)
            all_fps.append(fp)

    return all_fragments, all_fps


### ------------------ Similarity Analysis ------------------ ###

def compare_target_to_population(target_fps, population_fps, population_fragments, top_n=3):
    """Compares target fragments to population and returns the best-matching fragments."""
    results = []
    matched_fragments = []  # ✅ Store the matched molecular fragments

    for i, target_fp in enumerate(target_fps):
        similarities = np.array([sim.compute_similarity_score(target_fp, pop_fp) for pop_fp in population_fps])
        top_indices = np.argsort(similarities)[-top_n:][::-1]  # Get indices of top N matches
        top_scores = similarities[top_indices]

        # ✅ Store the top-matching fragments
        matched_fragments.extend([population_fragments[idx] for idx in top_indices if top_scores[0] >= 0.9990])

        results.append({
            "Target Fragment": f"Fragment {i+1}",
            "Match Found": any(score >= 0.9990 for score in top_scores),
            "Top 1 Similarity": top_scores[0] if len(top_scores) > 0 else None,
            "Top 2 Similarity": top_scores[1] if len(top_scores) > 1 else None,
            "Top 3 Similarity": top_scores[2] if len(top_scores) > 2 else None,
        })

    return pd.DataFrame(results), matched_fragments  # ✅ Return matched fragments


### ------------------ Main Execution ------------------ ###
start = time.time()

# **Step 1: Process the Target Molecule**
target_entry_id = "ABAHIW"
csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry(target_entry_id)
target_molecule = component_of_interest(entry.molecule)
target_fragments = get_fragments(target_molecule)

# Compute fingerprints
target_fps = [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in target_fragments]

# Find unique target fragments
unique_target_fragments, target_fps = find_unique_fragments_optimized(target_fragments, target_fps)

# Save unique target fragments
with MoleculeWriter("target_unique_fragments.sdf") as writer:
    for frag in unique_target_fragments:
        writer.write(frag)

# **Step 2: Process Population**
population_file = "targets/init_populations_protons/1_ABAHIW_init_pop.txt"
with open(population_file, "r") as f:
    population_entry_ids = [line.split()[0] for line in f]

population_entry_ids = population_entry_ids[:100]  # Debugging
population_fragments, population_fps = parallel_population_processing(population_entry_ids)

# **Step 3: Compare Target to Population**
comparison_df, matched_fragments = compare_target_to_population(target_fps, population_fps, population_fragments)

# **Step 4: Save Matching Population Fragments**
with MoleculeWriter("matched_population_fragments.sdf") as writer:
    for frag in matched_fragments:
        writer.write(frag)

comparison_df.to_csv("fragment_comparison_results.csv", index=False)

print(f"\n✅ Process completed in {time.time() - start:.2f} seconds.")
