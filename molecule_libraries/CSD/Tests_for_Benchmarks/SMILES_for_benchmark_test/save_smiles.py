from ccdc import io
import pandas as pd
import time
import os
import multiprocessing as mp

def process_entries(entry_subset, entries):
    """
    Worker function to process a subset of CSD entries.
    Returns a list of SMILES corresponding to the input entries.
    """
    smiles = []
    csd_reader = io.EntryReader('CSD')
    for entry in csd_reader.entries():
        if entry.identifier in entry_subset:
            smiles.append(entry.molecule.heaviest_component.smiles)
    return smiles

def main():
    start = time.time()
    print(os.getcwd())
    
    # Read the CSV file
    df = pd.read_csv('./molecule_libraries/CSD/final_structures.csv')

    # Collect all the entries in the DataFrame
    entries = df['Identifier'].tolist()

    # Initialize CSD reader
    csd_reader = io.EntryReader('CSD')

    # Number of CPU cores to use for multiprocessing
    num_workers = mp.cpu_count()

    # Split the workload into chunks (equal parts for each core)
    chunk_size = len(entries) // num_workers
    entry_chunks = [entries[i:i + chunk_size] for i in range(0, len(entries), chunk_size)]

    # Start multiprocessing pool
    with mp.Pool(processes=num_workers) as pool:
        # Map each chunk to the worker function
        results = pool.starmap(process_entries, [(chunk, entries) for chunk in entry_chunks])

    # Flatten the list of results
    smiles = [smile for sublist in results for smile in sublist]

    # Add SMILES to DataFrame
    df['SMILES'] = smiles

    # Save the updated DataFrame to a CSV file
    df.to_csv('final_structures_with_smiles.csv', index=False)
    
    print("Done!")
    print(f'Total time elapsed: {(time.time() - start)/60} minutes')

if __name__ == "__main__":
    main()
