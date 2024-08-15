from ccdc import io
from hsr import similarity as sim
from hsr import fingerprint as fp
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count, Manager
import time
import re
import os

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
    
if __name__ == '__main__':
    start = time.time()
    csd_reader = io.EntryReader('CSD')
    
    # # POLYMERIZATION TARGET
    # target_molecule_1 = csd_reader.entry('GUDQOL').molecule
    # target_fp_3d_1, target_fp_4d_1 = hsr_fingerprint(target_molecule_1)
    
    # # OLEFIN METATHESIS TARGET
    # target_molecule_2 = csd_reader.entry('TITTUO').molecule
    # target_fp_3d_2, target_fp_4d_2 = hsr_fingerprint(target_molecule_2)
    
    # # PNP LIGAND TARGET
    # target_molecule_3 = csd_reader.entry('NIVHEJ').molecule
    # target_fp_3d_3, target_fp_4d_3 = hsr_fingerprint(target_molecule_3)
    
    # # # SALEN TYPE TARGET
    # target_molecule_4 = csd_reader.entry('XIDTOW').molecule
    # target_fp_3d_4, target_fp_4d_4 = hsr_fingerprint(target_molecule_4)
    
    # ION PAIR TARGET
    file = 'targets/ion_pairs/ion_pair_2.sdf'
    mol_reader = io.MoleculeReader(file)
    target_molecule_1 = mol_reader[0]
    target_fp_3d_1, target_fp_4d_1 = hsr_fingerprint(target_molecule_1)
    
    # target_molecules = [target_molecule_1, target_molecule_2, target_molecule_3, target_molecule_4]
    target_molecules = [target_molecule_1]
    # target_fps_3d = [target_fp_3d_1, target_fp_3d_2, target_fp_3d_3, target_fp_3d_4]
    # target_fps_4d = [target_fp_4d_1, target_fp_4d_2, target_fp_4d_3, target_fp_4d_4]
    target_fps_3d = [target_fp_3d_1]
    target_fps_4d = [target_fp_4d_1]
    # results_files_3d = ['similarity_results_polymerization_3d.txt', 'similarity_results_grubbs_3d.txt', 'similarity_results_PNP_3d.txt', 'similarity_results_salen_3d.txt']
    # results_files_4d = ['similarity_results_polymerization_4d.txt', 'similarity_results_grubbs_4d.txt', 'similarity_results_PNP_4d.txt', 'similarity_results_salen_4d.txt']
    results_files_3d = ['similarity_results_ion_pair_3d.txt']
    results_files_4d = ['similarity_results_ion_pair_4d.txt']
    
    
    df = pd.read_csv('CSD_similarity/final_structures_with_fp.csv')
    manager = Manager()
    fingerprint_dict_3d = manager.dict(pd.Series(df.Fingerprint_3D.values, index=df.Identifier).to_dict())
    fingerprint_dict_4d = manager.dict(pd.Series(df.Fingerprint_4D.values, index=df.Identifier).to_dict())
    entries = list(fingerprint_dict_3d.keys())

    
    for target_fp_3d, target_fp_4d, result_file_3d, result_file_4d in zip(target_fps_3d, target_fps_4d, results_files_3d, results_files_4d):
        print(f"Calculating similarity for {len(entries)} entries with {cpu_count()} CPUs...")
        with Pool(processes=cpu_count(), initializer=init_process, initargs=(target_fp_3d, target_fp_4d, fingerprint_dict_3d, fingerprint_dict_4d)) as pool:
            results = pool.map(calculate_similarity, entries)
        print(f'Finished calculating similarity for {len(results)} entries')
        
        # Get results for 3D and 4D fingerprints in separate lists
        results_3d = [(result[0], result[1]) for result in results]
        results_4d = [(result[0], result[2]) for result in results]

        # Sort results by the maximum similarity score in descending order
        results_3d = sorted(results_3d, key=lambda x: max(x[1]) if x[1] else 0, reverse=True)
        results_4d = sorted(results_4d, key=lambda x: max(x[1]) if x[1] else 0, reverse=True)
        
        # Save the results in the two files
        print(f'Saving 3d results to {result_file_3d}...')
        with open(result_file_3d, 'w') as file:
            for result in results_3d:
                file.write(f"{result[0]}: {result[1]}\n")
        
        print(f'Saving 4d results to {result_file_4d}...')
        with open(result_file_4d, 'w') as file:
            for result in results_4d:
                file.write(f"{result[0]}: {result[1]}\n")
                
    print(f'Total time elapsed: {(time.time() - start)/60:.2f} minutes')
