from ccdc import io
import pandas as pd
from multiprocessing import Pool, cpu_count
import time
from ccdc import search


def get_csd_fingerprint(molecule):
    sim_search = search.SimilaritySearch(molecule)
    fp_builder = sim_search._fp
    fp = fp_builder.similarity_fingerprint(molecule._molecule)
    return fp

def calculate_similarity(data):
    query_identifier, target_identifier = data
    csd_reader = io.EntryReader('CSD')
    target_fp = get_csd_fingerprint(csd_reader.entry(target_identifier).molecule)
    query_molecule = csd_reader.entry(query_identifier).molecule
    similarities = []

    for component in query_molecule.components:
        component_fp = get_csd_fingerprint(component)
        
        similarity = component_fp.tanimoto(target_fp)
        similarities.append(similarity)

    return query_identifier, similarities

if __name__ == '__main__':
    start = time.time()
    csd_reader = io.EntryReader('CSD')

    # POLYMERIZATION TARGET
    identifier_1 = 'GUDQOL'
    target_molecule_1 = csd_reader.entry('GUDQOL')
    
    # OLEFIN METATHESIS TARGET
    identifier_2 = 'TITTUO'
    target_molecule_2 = csd_reader.entry('TITTUO')

    # PNP LIGAND TARGET
    identifier_3 = 'NIVHEJ'
    target_molecule_3 = csd_reader.entry('NIVHEJ')

    # SALEN TYPE TARGET
    identifier_4 = 'XIDTOW'
    target_molecule_4 = csd_reader.entry('XIDTOW')

    target_molecules = [target_molecule_1, target_molecule_2, target_molecule_3, target_molecule_4]
    # target_fps = [get_csd_fingerprint(x) for x in target_molecules]
    target_identifiers = [identifier_1, identifier_2, identifier_3, identifier_4]
    results_files_smiles = ['similarity_results_polymerization_csd.txt', 
                            'similarity_results_grubbs_csd.txt', 
                            'similarity_results_PNP_csd.txt',
                            'similarity_results_salen_csd.txt']
    df = pd.read_csv('final_structures_with_fp.csv')
    entries = df['Identifier'].tolist()
    # get only the first 100 entries for testing

    for target_identifier, result_file in zip(target_identifiers, results_files_smiles):
        print(f"Calculating similarity for {len(entries)} entries with {cpu_count()} CPUs...")
        with Pool(processes=cpu_count()) as pool:
            data = [(entry, target_identifier) for entry in entries]
            results = pool.map(calculate_similarity, data)
        print(f'Finished calculating similarity for {len(results)} entries')

        # Sort the results by similarity score
        results = sorted(results, key=lambda x: max(x[1]), reverse=True)

        # Save the results in the two files
        print(f'Saving smiles similarity results to {result_file}...')
        with open(result_file, 'w') as file:
            for result in results:
                file.write(f"{result[0]}: {result[1]}\n")

    print(f'Total time elapsed: {(time.time() - start)/60:.2f} minutes')
