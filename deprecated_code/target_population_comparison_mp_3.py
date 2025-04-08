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
    """Parallel function for processing a molecule: returns only fingerprints (no molecule objects)."""
    try:
        csd_reader = io.EntryReader('CSD')
        entry = csd_reader.entry(entry_id)
        molecule = component_of_interest(entry.molecule)
        fragments = get_fragments(molecule)

        # Compute fingerprints only (no molecule objects)
        fragment_fps = [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in fragments]

        return fragment_fps  # ✅ Only return fingerprints
    except Exception as e:
        print(f"❌ Error processing {entry_id}: {e}")
        return []




def parallel_population_processing(entry_ids, num_workers=mp.cpu_count()):
    """Processes all molecules in parallel, returning only fingerprints."""
    with mp.Pool(num_workers) as pool:
        all_fps = pool.map(process_molecule_parallel, entry_ids)
    
    return [fp for sublist in all_fps for fp in sublist]


### ------------------ Similarity Analysis ------------------ ###

def compare_target_to_population(target_fps, population_fps, top_n=3):
    """Compares target fragments to population and returns matches."""
    results = []
    found_target_fps = set()  # Track matched target fragments
    
    for i, target_fp in enumerate(target_fps):
        similarities = np.array([sim.compute_similarity_score(target_fp, pop_fp) for pop_fp in population_fps])
        top_indices = np.argsort(similarities)[-top_n:][::-1]
        top_scores = similarities[top_indices]

        # Convert fingerprint to a tuple before adding to set
        fp_tuple = tuple(target_fp) if isinstance(target_fp, np.ndarray) else tuple(target_fp)

        if any(score >= 0.9990 for score in top_scores):
            found_target_fps.add(fp_tuple)

        results.append({
            "Target Fragment": f"Fragment {i+1}",
            "Match Found": any(score >= 0.9990 for score in top_scores),
            "Top 1 Similarity": top_scores[0] if len(top_scores) > 0 else None,
            "Top 2 Similarity": top_scores[1] if len(top_scores) > 1 else None,
            "Top 3 Similarity": top_scores[2] if len(top_scores) > 2 else None,
        })

    return pd.DataFrame(results), len(found_target_fps)  # Return number of matched fragments

### ------------------ Main Execution ------------------ ###
start = time.time()

# **Step 1: Process the Target Molecule**
target_entry_id = "ABAHIW"

# ✅ Process in main process: Store Molecule objects (to avoid multiprocessing issues)
csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry(target_entry_id)
target_molecule = component_of_interest(entry.molecule)
target_fragments = get_fragments(target_molecule)  # ✅ Store molecules in main process

# ✅ Compute fingerprints in parallel (avoids pickling Molecule objects)
target_fps = [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in target_fragments]

# ✅ Find unique fragments and their corresponding fingerprints
unique_target_fragments, unique_target_fps = find_unique_fragments_optimized(target_fragments, target_fps)

# **Save Unique Target Fragments to SDF**
with MoleculeWriter("target_unique_fragments.sdf") as writer:
    for frag in unique_target_fragments:
        writer.write(frag)
# **Step 2: Process the Initial Population**
population_file = "targets/init_populations_protons/1_ABAHIW_init_pop.txt"

with open(population_file, "r") as f:
    population_entry_ids = [line.split()[0] for line in f]

population_entry_ids = population_entry_ids[:100]  # For debugging
# population_entry_ids.append('ABAHIW')
batch_size = len(population_entry_ids) // 8 + 1

found_fragments = set()
all_population_fps = []

# **Set to track all unique fingerprints in the population**
seen_population_fps = set()
all_population_fps = []

# **Process in Batches**
for i in range(0, len(population_entry_ids), batch_size):
    batch = population_entry_ids[i : i + batch_size]
    print(f"🔹 Processing batch {i//batch_size + 1}/{(len(population_entry_ids) // batch_size) + 1}")
    
    batch_fps = parallel_population_processing(batch)
    
    # ✅ Fix: Ensure uniqueness across ALL population fragments
    batch_fps = find_unique_fragments_optimized(fps=batch_fps)  
    
    # ✅ Now, ensure uniqueness across all batches by checking against seen_population_fps
    filtered_batch_fps = []
    for fingerprint in batch_fps:
        fp_tuple = tuple(fingerprint.tolist()) if isinstance(fingerprint, np.ndarray) else tuple(fingerprint)
        if fp_tuple not in seen_population_fps:
            seen_population_fps.add(fp_tuple)  # Add to the global unique set
            filtered_batch_fps.append(fingerprint)  # Store unique fragment

    # Compare to target fragments
    _, n_matches = compare_target_to_population(target_fps, filtered_batch_fps)

    if n_matches == len(target_fps):
        print("\n✅ Found all unique target fragments before finishing population analysis. Stopping early.")
        break

    # Store globally unique population fingerprints
    all_population_fps.extend(filtered_batch_fps)




# **Step 3: Compare Target to Population**
comparison_df, n_matches = compare_target_to_population(target_fps, all_population_fps)

# **Step 4: Save Results**
comparison_df.to_csv("fragment_comparison_results.csv", index=False)

# **Step 5: Display Results**
print("\n✅ Similarity analysis completed. Results saved to 'fragment_comparison_results.csv'.")
print(comparison_df)

# **Summary Statistics**
num_population_molecules = len(population_entry_ids)
num_unique_fragments = len(target_fps)
num_unique_population_fragments = len(all_population_fps)

print(f"\nPopulation size: {num_population_molecules} molecules.")
print(f"Unique fragments in target molecule: {num_unique_fragments}.")
print(f"Unique fragments in population: {num_unique_population_fragments}.")
print(f"Number of unique fragments in target molecule present in population: {n_matches}/{num_unique_fragments}.")

# **Time taken**
print(f"\nTime elapsed: {time.time() - start:.2f} seconds.")
