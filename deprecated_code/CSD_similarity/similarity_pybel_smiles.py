import os
import time
import json
from multiprocessing import Pool, cpu_count
from openbabel import pybel
from ccdc import io

# Global variable to hold the JSON data
json_data = None

def init_process(json_file):
    global json_data
    with open(json_file, 'r') as f:
        json_data = json.load(f)

def calculate_similarity(data):
    entry, target_smiles = data
    similarities = []

    # Compute the target fingerprint within the worker process
    target_molecule = pybel.readstring('smi', target_smiles)
    target_fp = target_molecule.calcfp("FP2")

    # Retrieve the smiles from the global json_data
    entry_smiles = json_data.get(entry, [])

    for smiles_str in entry_smiles:
        try:
            ob_molecule = pybel.readstring('smi', smiles_str)
        except Exception as e:
            print(f"Error reading smiles string: {smiles_str} in entry: {entry}")
            print(e)
            break 
        component_fp = ob_molecule.calcfp("FP2")
        # Calculate the Tanimoto similarity between the two fingerprints
        similarity = component_fp | target_fp
        similarities.append(similarity)

    return entry, similarities

if __name__ == '__main__':
    start = time.time()
    json_file = 'CSD_similarity/entry_smiles_pybel.json'

    # POLYMERIZATION TARGET
    identifier_1 = 'GUDQOL'
    # OLEFIN METATHESIS TARGET
    identifier_2 = 'TITTUO'
    # PNP LIGAND TARGET
    identifier_3 = 'NIVHEJ'
    # SALEN TYPE TARGET
    identifier_4 = 'XIDTOW'

    # Get the list of smiles for the target molecules from the json file
    with open(json_file, 'r') as f:
        smiles = json.load(f)  
        target_smiles_1 = smiles[identifier_1][0]
        target_smiles_2 = smiles[identifier_2][0]
        target_smiles_3 = smiles[identifier_3][0]
        target_smiles_4 = smiles[identifier_4][0]

    target_smiles_list = [target_smiles_1, target_smiles_2, target_smiles_3, target_smiles_4]
    dir = 'CSD_similarity'
    results_files_smiles = [f'{dir}/similarity_results_polymerization_pybel.txt',
                            f'{dir}/similarity_results_grubbs_pybel.txt',
                            f'{dir}/similarity_results_PNP_pybel.txt',
                            f'{dir}/similarity_results_salen_pybel.txt']

    entries = list(smiles.keys())
    # for testing purposes get the first 10 entries
    # entries = entries[:10]
    # # add specific entries for testing
    # entries.extend([identifier_1, identifier_2, identifier_3, identifier_4])

    for target_smiles, result_file in zip(target_smiles_list, results_files_smiles):
        n_cpus = 7
        print(f"Calculating similarity for {len(entries)} entries with {n_cpus} CPUs...")

        with Pool(processes=n_cpus, initializer=init_process, initargs=(json_file,)) as pool:
            data = [(entry, target_smiles) for entry in entries]
            results = pool.map(calculate_similarity, data)

        print(f'Finished calculating similarity for {len(results)} entries')
        
        # Get entries with empty lists (count the elements in the list)
        empty_results = [result for result in results if not result[1]]
        # Remove empty results: check if the lists are not empty
        results = [result for result in results if result[1]]

        results = sorted(results, key=lambda x: max(x[1]), reverse=True)

        print(f'Saving smiles similarity results to {result_file}...')
        with open(result_file, 'w') as file:
            for result in results:
                file.write(f"{result[0]}: {result[1]}\n")

    print(f'Total time elapsed: {(time.time() - start)/60:.2f} minutes')
