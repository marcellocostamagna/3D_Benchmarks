import ccdc.io
import pandas as pd
from multiprocessing import Pool, cpu_count
import hsr
import time
import numpy as np

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

# Load the dataframe of viable structures
csv_file = f'./targets/viable_structures.csv'
targets_file = f'./targets/targets.txt'

# Load the dataframe
df = pd.read_csv(csv_file)
# For testing purposes get only the first n entries
df = df.head(100)

# Load the target molecules
# Get the fingerprints of the target molecules
target_entries = []
with open(targets_file, 'r') as f:
    for line in f:
        target_entries.append(line.strip())
target_molecules = [ccdc.io.EntryReader('CSD').entry(e).molecule.heaviest_component for e in target_entries]
target_fps = [get_fingerprint(m) for m in target_molecules]

# Get the identifiers of the entries
# Get the fingerprints of the entries
entries = df['Identifier'].tolist()
# Extract the heaviest component of each entry as a molecule
viable_molecules = [ccdc.io.EntryReader('CSD').entry(e).molecule.heaviest_component for e in entries]
entry_fps = [get_fingerprint(m) for m in viable_molecules]
# cretate a dictionary, entry: [molecule, fingerprint]
entry_dict = dict(zip(entries, zip(viable_molecules, entry_fps)))

# Calculate the similarity between the target molecules and the entries
mols_per_target = {}
for i, target_fp in enumerate(target_fps):
    mols_per_target[f'Target_{i}'] = []
    for entry_fp, entry in zip(entry_fps, entries):
        similarity = hsr.similarity.compute_similarity_score(target_fp, entry_fp)
        if similarity < 0.4:
            mols_per_target[f'Target_{i}'].append(entry)
            
# Collect the entries that have a similarity score lower than a threshold for all the targets (e.g., 0.3)
# aka, select the entries that appear in all the lists, of each target, in mols_per_target dictionary
final_mols = set.intersection(*map(set, mols_per_target.values()))
print(f'Number of final molecules: {len(final_mols)}')

# Save the entries that have a similarity score lower than a threshold in a txt file
with open('./targets/initial_population_test.txt', 'w') as f:
    for mol in final_mols:
        f.write(f'{mol}\n')




