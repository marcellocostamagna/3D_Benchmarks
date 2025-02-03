# Script to retrieve all the identifiers of molecules that show a similarity score 
# lower than a threshold to a given target molecule.

import ccdc.io
import pandas as pd
from multiprocessing import Pool, cpu_count
import hsr
import time
import numpy as np

# Function to generate an array from a molecule
def get_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], atom.atomic_number])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

# Function to generate a fingerprint from a molecule
def get_fingerprint(ccdcmol):
    atom_array = get_array_from_ccdcmol(ccdcmol)
    return hsr.fingerprint.generate_fingerprint_from_data(atom_array)

# Function to compute similarity for a chunk of entries
def process_chunk(chunk_args):
    target_fp, entries_chunk = chunk_args
    results = {}
    entry_reader = ccdc.io.EntryReader('CSD')

    for entry in entries_chunk:
        try:
            molecule = entry_reader.entry(entry).molecule.heaviest_component
            entry_fp = get_fingerprint(molecule)
            similarity = hsr.similarity.compute_similarity_score(target_fp, entry_fp)
            if similarity < 0.3:
                results[entry] = similarity
        except Exception as e:
            print(f"Error processing entry {entry}: {e}")
    return results

if __name__ == "__main__":
    global_start_time = time.time() 

    # Load the dataframe of viable structures
    csv_file = './targets/viable_structures_last.csv'
    df = pd.read_csv(csv_file)

    # Remove the entries with no explicit hydrogen atoms
    df = df[df['Explicit Hs'] == True]
    print(f'Number of entries with explicit Hs: {len(df)}')

    # Limit for testing purposes
    # df = df.head(100000)

    # Load the target molecule
    target_entry = 'ABAHIV'
    target_molecule = ccdc.io.EntryReader('CSD').entry(target_entry).molecule.heaviest_component
    target_fp = get_fingerprint(target_molecule)

    # Get the list of entries
    entries = df['Identifier'].tolist()

    # Split entries into chunks for multiprocessing
    num_chunks = 7 #cpu_count()
    chunk_size = len(entries) // num_chunks + 1
    entry_chunks = [entries[i:i + chunk_size] for i in range(0, len(entries), chunk_size)]

    print(f"Calculating similarity scores for {len(entries)} entries using {num_chunks} processes...")
    start_time = time.time()

    # Run multiprocessing pool
    with Pool(processes=num_chunks) as pool:
        chunk_args = [(target_fp, chunk) for chunk in entry_chunks]
        results_list = pool.map(process_chunk, chunk_args)

    # Combine results from all processes
    mols_per_target = {}
    for partial_results in results_list:
        mols_per_target.update(partial_results)

    # Final list of molecules
    final_mols = list(mols_per_target.keys())
    print(f'Number of final molecules: {len(final_mols)}')

    # Save the entries that have a similarity score lower than a threshold in a text file
    output_file = f'./targets/initial_population_test_{target_entry}.txt'
    with open(output_file, 'w') as f:
        for mol in final_mols:
            f.write(f"{mol}\n")
    
    print(f'Results saved to {output_file}')
    print(f'Total time elapsed: {time.time() - global_start_time:.2f} seconds')
