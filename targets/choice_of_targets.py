# Script to perform the creation of a initial population 
# form the choice of targets and a similarity cutoff and 
# check if the initial popultaion is viable; meaning all the
# representatitve molecular patterns present in the targets 
# are also present in the initial population

import ccdc.io
import pandas as pd
from multiprocessing import Pool, cpu_count
import hsr
import time
import numpy as np
from pattern_queries import *
import pandas as pd

# functions
# For getting the arrays of the target molecules
def get_array_from_ccdcmol(ccdcmol):
    # Iterate over the atoms in the molecule
    atom_array = []
    for atom in ccdcmol.atoms:
        # Append the coordinates and the atomic number of the atom (on each row)
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], atom.atomic_number])    
    # Centering the data
    atom_array -= np.mean(atom_array, axis=0)
    return atom_array

def get_fingerprint(ccdcmol):
    # Get the array of the molecule
    atom_array = get_array_from_ccdcmol(ccdcmol)
    # Generate the fingerprint
    return hsr.fingerprint.generate_fingerprint_from_data(atom_array)

# Constants
SIMILARITY_THRESHOLD = 0.4

# PIPELINE:

# 1. Definition of patterns of interest (Done)

# 2. Targets, read from a file (only manual input)
targets_file = f'./targets/targets.txt'
target_entries = []
with open(targets_file, 'r') as f:
    for line in f:
        target_entries.append(line.strip())
target_molecules = [ccdc.io.EntryReader('CSD').entry(e).molecule.heaviest_component for e in target_entries]

# 3. Viable molecules, read from a file (Done)
csv_file = f'./targets/viable_structures.csv'
df = pd.read_csv(csv_file)
entries = df['Identifier'].tolist()
viable_molecules = [ccdc.io.EntryReader('CSD').entry(e).molecule.heaviest_component for e in entries]

# 4. Initial population (potential) : filtering of viable molecules

# TODO: Make the simlarity filtering in multiprocessing

# Calculate the similarity between the target molecules and the entries
target_fps = [get_fingerprint(m) for m in target_molecules]
entry_fps = [get_fingerprint(m) for m in viable_molecules]

mols_per_target = {}
for i, target_fp in enumerate(target_fps):
    mols_per_target[f'Target_{i}'] = []
    for entry_fp, entry in zip(entry_fps, entries):
        similarity = hsr.similarity.compute_similarity_score(target_fp, entry_fp)
        if similarity < SIMILARITY_THRESHOLD:
            mols_per_target[f'Target_{i}'].append(entry)
            
# Collect the entries that have a similarity score lower than a threshold for all the targets
potential_initial_population = set.intersection(*map(set, mols_per_target.values()))

# 5. Check if the initial population is viable
#    - All the representative molecular patterns are at least present in one molecule of the initial population
#    - (Optional) Check the frequency of the patterns in the initial population

# List of queries for the patterns of interest (from pattern_queries.py)
# TODO: As the substructure search is computationally expensive, it has to be done in multiprocessing




# 6. Alternatives:
#   If the initial population is not viable:
#   - Highlight the patterns (and the associated targets) that are not present in the initial population
#   - Demand a new choice of targets
#   
#   If the initial population is viable:
#   - Return the initial population
#   - (Optional) Return the frequency of the patterns in the initial population

