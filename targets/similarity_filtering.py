# Script to retreive all the identifiers of molecules that show a similarity score 
# lower than a threshold to a given target molecule.

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

global_start_time = time.time() 

# Load the dataframe of viable structures
csv_file = f'./targets/viable_structures.csv'
# target_file = f'./targets/target.txt'

# Load the dataframe
df = pd.read_csv(csv_file)

# Remove the entries with no explicit hydrogen atoms
df = df[df['Explicit Hs'] == True]
print(f'Number of entries with explicit Hs: {len(df)}')

# For testing purposes get only the first n entries
df = df.head(10000)

# Load the target molecules
# Get the fingerprints of the target molecules
target_entry = 'AABHTZ'
# with open(target_file, 'r') as f:
#     for line in f:
#         target_entry.append(line.strip())
target_molecule = ccdc.io.EntryReader('CSD').entry(target_entry).molecule.heaviest_component
target_fp = get_fingerprint(target_molecule)

# Get the identifiers of the entries
# Get the fingerprints of the entries
entries = df['Identifier'].tolist()
# Extract the heaviest component of each entry as a molecule
viable_molecules = [ccdc.io.EntryReader('CSD').entry(e).molecule.heaviest_component for e in entries]
entry_fps = [get_fingerprint(m) for m in viable_molecules]
# cretate a dictionary, entry: [molecule, fingerprint]
entry_dict = dict(zip(entries, zip(viable_molecules, entry_fps)))

# Calculate the similarity between the target molecule and the entries
print(f'Calculating the similarity scores of {len(entries)} entries...')
mols_per_target = {}
for entry_fp, entry in zip(entry_fps, entries):
    similarity = hsr.similarity.compute_similarity_score(target_fp, entry_fp)
    if similarity < 0.3:
        mols_per_target[entry] = similarity
            
final_mols = list(mols_per_target.keys())
print(f'Number of final molecules: {len(final_mols)}')

# Save the entries that have a similarity score lower than a threshold in a txt file
with open('./targets/initial_population_test.txt', 'w') as f:
    for mol in final_mols:
        f.write(f'{mol}\n')
        
print(f'Time elapsed: {time.time() - global_start_time:.2f} seconds')




