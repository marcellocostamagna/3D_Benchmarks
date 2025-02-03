# Script to search for a single query across the entries in viable_structures.csv
# and save the matching identifiers to a text file.

import pandas as pd
import ccdc.io
import ccdc.search
from multiprocessing import Pool, cpu_count
from pattern_queries import *
from pattern_queries import METALS, NON_METALS
import time
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def search_molecules(query, identifiers):
    '''
    Searches for a query in a subset of the CSD defined by 
    the identifiers of the entries.
    The search is operated only on the heaviest component of the molecule.

    Returns the results of the search as a list of the hits' identifiers.
    '''
    results_ids = set()
    entry_reader = ccdc.io.EntryReader('CSD')   
    entries_dict = {id: entry_reader.entry(id) for id in identifiers}
    heaviest_molecules = {
        id: entry.molecule.heaviest_component for id, entry in entries_dict.items()
    }

    search = ccdc.search.SubstructureSearch()
    search.add_substructure(query)      
    results = search.search(list(heaviest_molecules.values()))

    for hit in results:
        matched_molecule = hit.molecule
        for entry_id, molecule in heaviest_molecules.items():
            if molecule == matched_molecule:
                results_ids.add(entry_id)
                break            
    return list(results_ids)

def run_query(args):
    '''
    Helper function for multiprocessing to run a query and return the result.
    '''
    
    query_name, identifiers = args
    
    # print(f"Running query: {query_name}")
    start_time = time.time()

    # Retrieve the query function dynamically
    query_function = queries_dict[query_name]
    results = search_molecules(query_function, identifiers)
    print(f"Completed query: {query_name} in {time.time() - start_time:.2f} seconds")
    return results

# Define the queries dictionary
global queries_dict
queries_dict = {
    'chelating_ligands': chelating_ligands(METALS),
    'mono_bridged_complexes': mono_bridged_complexes(METALS, NON_METALS, 6),
    'multi_bridged_complexes': multi_bridged_complexes(METALS, NON_METALS, 6),
    'linear_ligands': linear_ligands(METALS, NON_METALS),
    'metal_allenes': metal_allenes(METALS, NON_METALS),
    'dative_bonds': dative_bonds(METALS, NON_METALS),
    'dative_bonds_in_rings': dative_in_rings(METALS, NON_METALS, 6),
    'dative_bonds_in_conjugated_chains': dative_in_conjugated_chains(METALS, NON_METALS),
    'dative_bonds_in_conjugated_rings': dative_in_conjugated_rings(METALS, NON_METALS),
    'cyclic_multihapto_ligands': cyclic_multihapto_ligands(METALS, NON_METALS),
    'sigma_bonds': sigma_bonds(METALS),
    'sigma_bonds_in_rings': sigma_bonds_in_rings(METALS, NON_METALS, 6),
    'saturated_rings': saturated_rings(NON_METALS, 8),
    'unsaturated_rings': unsaturated_rings(NON_METALS, 12),
    'unsaturated_conjugated_rings': unsaturated_conjugated_rings(NON_METALS, 8),
    'aromatic_rings': aromatic_rings(NON_METALS),
    'saturated_chains': saturated_chains(NON_METALS),
    'unsaturated_chains': unsaturated_chains(NON_METALS),
    'unsaturated_conjugated_chains': unsaturated_conjugated_chains(NON_METALS),
    'allenes': allenes(NON_METALS),
}

if __name__ == "__main__":
    global_start_time = time.time()
    
    # Load the dataframe of viable structures
    csv_file = './targets/viable_structures.csv'
    df = pd.read_csv(csv_file)

    # print(f'Number of entries: {len(df)}')

    # Remove the entries with no explicit hydrogen atoms
    df = df[df['Explicit Hs'] == True]
    # print(f'Number of entries with explicit Hs: {len(df)}')

    # Get the identifiers from the dataframe
    identifiers = df['Identifier'].tolist()

    # For testing purposes, limit the number of identifiers
    # identifiers = identifiers[:100000]

    # DEFINE YOUR QUERY HERE 
    # Get the first query from the queries dictionary
    query_name = 'chelating_ligands'

    # Run the search for the current query
    print(f"Running query: {query_name} on {len(identifiers)} entries")
    start_time = time.time()
    
    # Split the research into chunks and run them in parallel
    num_chunks = cpu_count()
    chunk_size = len(identifiers) // num_chunks + 1
    identifiers_chunks = [identifiers[i:i + chunk_size] for i in range(0, len(identifiers), chunk_size)]
    
    chunk_args = [(query_name, chunk) for chunk in identifiers_chunks]
    
    with Pool(processes=cpu_count()) as pool:
        # Run the search for the query
        results_list = pool.map(run_query, chunk_args)
        
    # Combine results from all chunks
    results = set()
    for partial_results in results_list:
        results.update(partial_results)
    
    # Save results to a text file
    output_file = f'./targets/patterns_search_results/{query_name}.txt'
    with open(output_file, 'w') as f:
        for identifier in results:
            f.write(f"{identifier}\n")
            
    print(f'Found {len(results)} hits in {len(identifiers)} entries for query {query_name}')
    print(f"Results for query '{query_name}' saved to {output_file}")
    print(f"Total time: {time.time() - global_start_time:.2f} seconds")
