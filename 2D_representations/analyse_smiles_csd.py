import pandas as pd
import time
from ccdc import io
import numpy as np
from multiprocessing import Pool, cpu_count
import json


def analyse_smiles(entry):
    """Analyse the smiles of a molecule."""
    csd_reader = io.EntryReader('CSD')
    molecule = csd_reader.entry(entry).molecule
    smiles = []
    for component in molecule.components:
        
        smiles.append(component.smiles)
    return entry, smiles

if __name__ == '__main__':
    print(f'Starting retrieving entries...')
    start = time.time()
    csd_reader = io.EntryReader('CSD') 
    df = pd.read_csv('CSD_similarity/final_structures_with_fp.csv')
    # Get the n_enties from the dataframe
    entries = df['Identifier'].values
    # entries = [csd_reader.entry(identifier) for identifier in entries]
    print(f'Entries retrieved: {len(entries)} in {(time.time()-start)/60} minutes')
    
    print(f"Extracting smiles for {len(entries)} entries with {cpu_count()} CPUs...")
    n_cpus = 7
    with Pool(processes=n_cpus) as pool:
        results = pool.map(analyse_smiles, entries)
    print(f'Finished extracting smiles for {len(results)} entries')
    
    ### ANALYSE THE RESULTS ###
    
    # create a dictionary with the entry and the smiles
    entry_smiles = {entry: smiles for entry, smiles in results}
    
    # Save the entry_smiles dictionary as a json file
    print(f'Saving the entry_smiles dictionary...')
    with open('entry_smiles.json', 'w') as f:
        json.dump(entry_smiles, f)
    print(f'entry_smiles dictionary saved')
        
    
    # Count the entries that have None in at least one of the smiles and save the entries
    entries_with_none = []
    for entry, smiles in entry_smiles.items():
        if None in smiles:
            entries_with_none.append(entry)
    
    print(f'Entries with None in at least one of the smiles: {len(entries_with_none)}')
    # print the first 10 entries (one on a new line) with None in at least one of the smiles
    print(entries_with_none[:10])
    print(f'Total time: {(time.time()-start)/60} minutes')
    
    
    
    