import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import time
from ccdc import io
from ccdc.molecule import Molecule
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

def find_unique_fragments_optimized(fps, threshold=0.9990):
    """Finds unique fragments efficiently using a hash-based set."""
    unique_fps = []
    seen = set()

    for fp in fps:
        fp_tuple = tuple(fp.tolist()) if isinstance(fp, np.ndarray) else tuple(fp)

        if fp_tuple in seen:
            continue

        is_unique = all(sim.compute_similarity_score(fp, stored_fp) < threshold for stored_fp in unique_fps)

        if is_unique:
            unique_fps.append(fp)
            seen.add(fp_tuple)

    return unique_fps


def process_molecule_parallel(entry_id):
    """Parallel version of process_molecule for multiprocessing."""
    try:
        csd_reader = io.EntryReader('CSD')
        entry = csd_reader.entry(entry_id)
        molecule = component_of_interest(entry.molecule)
        fragments = get_fragments(molecule)

        return [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in fragments]
    except Exception as e:
        print(f"❌ Error processing {entry_id}: {e}")
        return []


def parallel_population_processing(entry_ids, num_workers=mp.cpu_count()):
    """Processes all molecules in parallel using starmap for efficiency."""
    with mp.Pool(num_workers) as pool:
        all_fps = pool.map(process_molecule_parallel, entry_ids)
    
    return [fp for sublist in all_fps for fp in sublist]  # Flatten results


def process_population_in_chunks(population_entry_ids, batch_size=500):
    """Processes population molecules in batches to avoid memory overload."""
    all_population_fps = []

    for i in range(0, len(population_entry_ids), batch_size):
        batch = population_entry_ids[i : i + batch_size]
        print(f"🔹 Processing batch {i//batch_size + 1}/{(len(population_entry_ids) // batch_size) + 1}")
        
        batch_fps = parallel_population_processing(batch)
        all_population_fps.extend(batch_fps)

    return all_population_fps


### ------------------ Similarity Analysis ------------------ ###

def compare_target_to_population(target_fps, population_fps, top_n=3):
    """Compares target fragments to population and returns matches."""
    results = []
    
    for i, target_fp in enumerate(target_fps):
        similarities = np.array([sim.compute_similarity_score(target_fp, pop_fp) for pop_fp in population_fps])
        top_indices = np.argsort(similarities)[-top_n:][::-1]
        top_scores = similarities[top_indices]

        results.append({
            "Target Fragment": f"Fragment {i+1}",
            "Match Found": any(score >= 0.9990 for score in top_scores),
            "Top 1 Similarity": top_scores[0] if len(top_scores) > 0 else None,
            "Top 2 Similarity": top_scores[1] if len(top_scores) > 1 else None,
            "Top 3 Similarity": top_scores[2] if len(top_scores) > 2 else None,
        })

    return pd.DataFrame(results), sum(r["Match Found"] for r in results)


### ------------------ Main Execution ------------------ ###
start = time.time()

# **Step 1: Process the Target Molecule**
target_entry_id = "ABAHIW"
target_fps = find_unique_fragments_optimized(process_molecule_parallel(target_entry_id))

# **Step 2: Process the Initial Population**
population_file = "targets/init_populations_protons/1_ABAHIW_init_pop.txt"

with open(population_file, "r") as f:
    population_entry_ids = [line.split()[0] for line in f]

population_entry_ids = population_entry_ids[:100]
batch_size = len(population_entry_ids) // 8 + 1

# **Process in Batches**
population_fps = find_unique_fragments_optimized(process_population_in_chunks(population_entry_ids, batch_size))

# **Step 3: Compare Target to Population**
comparison_df, n_matches = compare_target_to_population(target_fps, population_fps)

# **Step 4: Save Results**
comparison_df.to_csv("fragment_comparison_results.csv", index=False)

# **Step 5: Display Results**
print("\n✅ Similarity analysis completed. Results saved to 'fragment_comparison_results.csv'.")
print(comparison_df)

print(f"\nPopulation size: {len(population_entry_ids)} molecules.")
print(f"Unique fragments in target molecule: {len(target_fps)}.")
print(f"Unique fragments in population: {len(population_fps)}.")
print(f"Number of unique fragments in target molecule present in population: {n_matches}/{len(target_fps)}.")
print(f"\nTime elapsed: {time.time() - start:.2f} seconds.")
