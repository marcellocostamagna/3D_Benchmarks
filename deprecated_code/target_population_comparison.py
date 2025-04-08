import os
import numpy as np
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc.diagram import DiagramGenerator
from hsr import fingerprint as fp
from hsr import similarity as sim
import time

# Ensure proper X11 configuration for Qt
os.environ["QT_QPA_PLATFORM"] = "xcb"
matplotlib.use("TkAgg")  # Ensure interactive Matplotlib backend


### ------------------ Molecule Processing Functions ------------------ ###

def get_array_from_ccdcmol(ccdcmol):
    """Extracts an array of atomic coordinates and atomic numbers."""
    atom_array = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return atom_array - np.mean(atom_array, axis=0)  # Centering data


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


def find_unique_fragments(fragments, similarity_matrix, threshold=0.9990):
    """Finds unique fragments by filtering out similar ones based on a similarity threshold."""
    unique_indices = []
    seen = set()

    for i in range(len(fragments)):
        if i in seen:
            continue
        unique_indices.append(i)
        seen.update(j for j in range(i + 1, len(fragments)) if similarity_matrix[i, j] >= threshold)
    
    return unique_indices


### ------------------ Fragmentation & Fingerprinting ------------------ ###

def process_molecule(entry_id):
    """Processes a molecule: fragmentation, uniqueness check, and fingerprinting."""
    csd_reader = io.EntryReader('CSD')
    entry = csd_reader.entry(entry_id)
    molecule = component_of_interest(entry.molecule)

    fragments = get_fragments(molecule)
    fragment_fps = [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in fragments]

    # Compute similarity matrix
    similarity_matrix = np.array([[sim.compute_similarity_score(fp1, fp2) for fp2 in fragment_fps] for fp1 in fragment_fps])

    # Get unique fragments
    unique_indices = find_unique_fragments(fragments, similarity_matrix)
    unique_fragments = [fragments[i] for i in unique_indices]
    unique_fps = [fragment_fps[i] for i in unique_indices]

    return unique_fragments, unique_fps


### ------------------ Similarity Analysis ------------------ ###

def compare_target_to_population(target_fps, population_fps, top_n=3):
    """
    Compares the unique fragments of a target molecule to the unique fragments of the initial population.

    Returns a DataFrame with:
    - Whether an exact match exists.
    - The top 3 most similar fragments from the population.
    """
    results = []
    
    for i, target_fp in enumerate(target_fps):
        similarities = [sim.compute_similarity_score(target_fp, pop_fp) for pop_fp in population_fps]
        top_indices = np.argsort(similarities)[-top_n:][::-1]  # Get top N matches
        top_scores = [similarities[idx] for idx in top_indices]

        results.append({
            "Target Fragment": f"Fragment {i+1}",
            "Match Found": any(score >= 0.9990 for score in top_scores),
            "Top 1 Similarity": top_scores[0] if len(top_scores) > 0 else None,
            "Top 2 Similarity": top_scores[1] if len(top_scores) > 1 else None,
            "Top 3 Similarity": top_scores[2] if len(top_scores) > 2 else None,
        })

    # get number of matches
    num_matches = sum(result["Match Found"] for result in results)
    return pd.DataFrame(results), num_matches


### ------------------ Main Execution ------------------ ###
start = time.time()

# **Step 1: Process the Target Molecule**
target_entry_id = "ABAHIW"  # Example CSD Entry ID
target_fragments, target_fps = process_molecule(target_entry_id)

# **Step 2: Process the Initial Population**

population_file = "targets/init_populations_protons/1_ABAHIW_init_pop.txt"

with open(population_file, "r") as f:
    # Get the entries from the first field of each line separated by a space
    population_entry_ids = [line.split()[0] for line in f]
    
# TESTING and debugging
# population_entry_ids = ["ACOVUL"]#, "AFECIA", "ABAHIW"]  # Example CSD Entry IDs for population
population_entry_ids = population_entry_ids[:100]  # Limit to first 10 entries for testing

all_population_fps = []

for entry_id in population_entry_ids:
    _, fps = process_molecule(entry_id)
    all_population_fps.extend(fps)  # Collect all unique fingerprints from population

# Compute a similarity matrix across ALL population fragments
population_similarity_matrix = np.array([[sim.compute_similarity_score(fp1, fp2) 
                                          for fp2 in all_population_fps] for fp1 in all_population_fps])

# Find globally unique fragments across the entire population
unique_population_indices = find_unique_fragments(all_population_fps, population_similarity_matrix)
population_fps = [all_population_fps[i] for i in unique_population_indices]

# **Step 3: Compare Target to Population**
comparison_df, n_matches = compare_target_to_population(target_fps, population_fps)

# **Step 4: Save Results**
comparison_df.to_csv("fragment_comparison_results.csv", index=False)

# **Step 5: Display Results**
print("\n✅ Similarity analysis completed. Results saved to 'fragment_comparison_results.csv'.")
print(comparison_df)

# Number of molecules in the population
num_population_molecules = len(population_entry_ids)
print(f"\nPopulation size: {num_population_molecules} molecules.")

# Number of unique fragments in the target molecule
num_unique_fragments = len(target_fragments)        
print(f"Unique fragments in target molecule: {num_unique_fragments}.")

# Number of unique fragments in the population
num_unique_population_fragments = num_unique_population_fragments = len(set(tuple(fp) for fp in population_fps))

print(f"Unique fragments in population: {num_unique_population_fragments}.")

# How many unique fragments in the target molecule are present in the population (x/tot unique frags)
print(f"Number of unique fragments in target molecule present in population: {n_matches}/{num_unique_fragments}.")

# Time taken
print(f"\nTime elapsed: {time.time() - start:.2f} seconds.")