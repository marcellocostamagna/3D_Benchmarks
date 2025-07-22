import os
import time
import json
import re
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count, Manager
from openbabel import pybel
from ccdc import io
from ccdc import search
from hsr import similarity as sim
from hsr import fingerprint as fp

# Functions from similarity.py and similarity_csd.py
def hsr_fingerprint(molecule):
    """Generates a fingerprint from a molecule object."""
    molecule_array = []
    molecule_array_1 = []
    for atom in molecule.atoms:
        if atom.coordinates is None:
            continue
        x, y, z = atom.coordinates
        atom_array = [x, y, z]
        atom_array_1 = np.array([x, y, z, np.sqrt(atom.atomic_number)])
        molecule_array.append(atom_array)
        molecule_array_1.append(atom_array_1)
    # Centering data
    molecule_array = molecule_array - np.mean(molecule_array, axis=0)
    molecule_array_1 = molecule_array_1 - np.mean(molecule_array_1, axis=0)
    return fp.generate_fingerprint_from_data(np.array(molecule_array)), fp.generate_fingerprint_from_data(np.array(molecule_array_1))

def extract_fps(fp_text):
    """Extracts fingerprints from a formatted string."""
    fps = {}
    pattern = r"(\d+):\s(\[.*?\]|None);?"
    matches = re.finditer(pattern, fp_text)
    for match in matches:
        comp_id, values = match.groups()
        if values != 'None':
            numeric_values = np.fromstring(values.strip('[]'), sep=',')
            fps[comp_id] = numeric_values
    return fps

def calculate_similarity(identifier):
    """Calculates similarities of the fingerprints against the target fingerprint."""
    fp_text_3d = fp_dict_3d[identifier]
    fp_text_4d = fp_dict_4d[identifier]
    fps_3d = extract_fps(fp_text_3d)
    fps_4d = extract_fps(fp_text_4d)
    similarities_3d = []
    similarities_4d = []

    for fp_3d, fp_4d in zip(fps_3d.values(), fps_4d.values()):
        similarities_3d.append(sim.compute_similarity_score(trgt_fp_3d, fp_3d))
        similarities_4d.append(sim.compute_similarity_score(trgt_fp_4d, fp_4d))
    
    return identifier, similarities_3d, similarities_4d

def init_process(target_fp_3d, target_fp_4d, fingerprint_dict_3d, fingerprint_dict_4d):
    """Initializer function to setup global variables."""
    global trgt_fp_3d
    global trgt_fp_4d
    global fp_dict_3d
    global fp_dict_4d
    trgt_fp_3d = target_fp_3d
    trgt_fp_4d = target_fp_4d
    fp_dict_3d = fingerprint_dict_3d
    fp_dict_4d = fingerprint_dict_4d

def init_process_pybel(json_file):
    global json_data
    with open(json_file, 'r') as f:
        json_data = json.load(f)

def calculate_similarity_pybel(data):
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
            # print(f"Error reading smiles string: {smiles_str} in entry: {entry}")
            print(f"Error reading smiles string in entry: {entry}. Error: {e}")
            break 
        component_fp = ob_molecule.calcfp("FP2")
        # Calculate the Tanimoto similarity between the two fingerprints
        similarity = component_fp | target_fp
        similarities.append(similarity)

    return entry, similarities

def get_csd_fingerprint(molecule):
    sim_search = search.SimilaritySearch(molecule)
    fp_builder = sim_search._fp
    fp = fp_builder.similarity_fingerprint(molecule._molecule)
    return fp

def calculate_similarity_csd(data):
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
    overall_start = time.time()
    csd_reader = io.EntryReader('CSD')
    json_file = 'CSD_similarity/entry_smiles_pybel.json'

    # Target identifiers
    target_identifiers = ['GUDQOL', 'TITTUO', 'NIVHEJ', 'XIDTOW']
    target_labels = ['polymerization', 'grubbs', 'PNP', 'salen']

    target_molecules = [csd_reader.entry(id).molecule for id in target_identifiers]
    target_fps_3d = []
    target_fps_4d = []
    for molecule in target_molecules:
        fp_3d, fp_4d = hsr_fingerprint(molecule)
        target_fps_3d.append(fp_3d)
        target_fps_4d.append(fp_4d)

    df = pd.read_csv('CSD_similarity/final_structures_with_fp.csv')
    manager = Manager()
    fingerprint_dict_3d = manager.dict(pd.Series(df.Fingerprint_3D.values, index=df.Identifier).to_dict())
    fingerprint_dict_4d = manager.dict(pd.Series(df.Fingerprint_4D.values, index=df.Identifier).to_dict())
    entries = list(fingerprint_dict_3d.keys())
    # for testing purposes get the first 10 entries
    # entries = entries[:10]

    results_files_3d = [f'similarity_results_{label}_3d.txt' for label in target_labels]
    results_files_4d = [f'similarity_results_{label}_4d.txt' for label in target_labels]
    results_files_pybel = [f'similarity_results_{label}_pybel.txt' for label in target_labels]
    results_files_csd = [f'similarity_results_{label}_csd.txt' for label in target_labels]

    # Load JSON data for pybel
    with open(json_file, 'r') as f:
        smiles = json.load(f)
    target_smiles_list = [smiles[id][0] for id in target_identifiers]

    for target_fp_3d, target_fp_4d, target_smiles, target_identifier, result_file_3d, result_file_4d, result_file_pybel, result_file_csd, label in zip(target_fps_3d, target_fps_4d, target_smiles_list, target_identifiers, results_files_3d, results_files_4d, results_files_pybel, results_files_csd, target_labels):
        target_start = time.time()

        # 3D Similarity calculations
        start_3d = time.time()
        print(f"Calculating 3D similarity for {len(entries)} entries with {cpu_count()} CPUs for target {label}...")
        with Pool(processes=cpu_count(), initializer=init_process, initargs=(target_fp_3d, target_fp_4d, fingerprint_dict_3d, fingerprint_dict_4d)) as pool:
            results = pool.map(calculate_similarity, entries)
        print(f'Saving 3D results to {result_file_3d}...')
        results_3d = [(result[0], result[1]) for result in results]
        results_3d = sorted(results_3d, key=lambda x: max(x[1]) if x[1] else 0, reverse=True)
        with open(result_file_3d, 'w') as file:
            for result in results_3d:
                file.write(f"{result[0]}: {result[1]}\n")
        time_3d = time.time() - start_3d
        print(f'Time taken for 3D similarity: {time_3d/60:.2f} minutes\n')

        # 4D Similarity calculations
        start_4d = time.time()
        print(f"Calculating 4D similarity for {len(entries)} entries with {cpu_count()} CPUs for target {label}...")
        with Pool(processes=cpu_count(), initializer=init_process, initargs=(target_fp_3d, target_fp_4d, fingerprint_dict_3d, fingerprint_dict_4d)) as pool:
            results = pool.map(calculate_similarity, entries)
        print(f'Saving 4D results to {result_file_4d}...')
        results_4d = [(result[0], result[2]) for result in results]
        results_4d = sorted(results_4d, key=lambda x: max(x[1]) if x[1] else 0, reverse=True)
        with open(result_file_4d, 'w') as file:
            for result in results_4d:
                file.write(f"{result[0]}: {result[1]}\n")
        time_4d = time.time() - start_4d
        print(f'Time taken for 4D similarity: {time_4d/60:.2f} minutes\n')

        # Pybel Similarity calculations
        start_pybel = time.time()
        print(f"Calculating pybel similarity for {len(entries)} entries with {cpu_count()} CPUs for target {label}...")
        with Pool(processes=cpu_count(), initializer=init_process_pybel, initargs=(json_file,)) as pool:
            data = [(entry, target_smiles) for entry in entries]
            results = pool.map(calculate_similarity_pybel, data)
        print(f'Saving pybel similarity results to {result_file_pybel}...')
        results_pybel = [result for result in results if result[1]]
        results_pybel = sorted(results_pybel, key=lambda x: max(x[1]), reverse=True)
        with open(result_file_pybel, 'w') as file:
            for result in results_pybel:
                file.write(f"{result[0]}: {result[1]}\n")
        time_pybel = time.time() - start_pybel
        print(f'Time taken for pybel similarity: {time_pybel/60:.2f} minutes\n')

        # CSD Similarity calculations
        start_csd = time.time()
        print(f"Calculating CSD similarity for {len(entries)} entries with {cpu_count()} CPUs for target {label}...")
        with Pool(processes=cpu_count()) as pool:
            data = [(entry, target_identifier) for entry in entries]
            results = pool.map(calculate_similarity_csd, data)
        print(f'Saving CSD similarity results to {result_file_csd}...')
        results_csd = sorted(results, key=lambda x: max(x[1]), reverse=True)
        with open(result_file_csd, 'w') as file:
            for result in results_csd:
                file.write(f"{result[0]}: {result[1]}\n")
        time_csd = time.time() - start_csd
        print(f'Time taken for CSD similarity: {time_csd/60:.2f} minutes\n')

        target_time_elapsed = time.time() - target_start
        print(f'Total time taken for target {label}: {target_time_elapsed/60:.2f} minutes\n')

    overall_time_elapsed = time.time() - overall_start
    print(f'Total time elapsed: {overall_time_elapsed/60:.2f} minutes')
