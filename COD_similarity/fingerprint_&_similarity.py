from ccdc import io
from hsr import fingerprint as fp
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count, Manager
import time
from hsr import similarity as sim
import re
import warnings


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
    # print(np.array(molecule_array))
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

def get_fingerprints(file_path):
    fp_dict = {}
    fp_dict_1 = {}
    mol_reader = io.MoleculeReader(file_path)
    molecule = mol_reader[0]
    for component in molecule.components:
        if len(component.atoms) < 4:
            fp, fp1 = None, None
        else:
            try:
                fp, fp1 = hsr_fingerprint(component) 
            except Exception as e:
                print(f"Error generating fingerprint for {molecule.identifier} in component {component.identifier}: {e}")
                fp, fp1 = None, None
        fp_dict[component.identifier] = fp
        fp_dict_1[component.identifier] = fp1
    fp_dict = '; '.join([f"{k}: {v}" for k, v in fp_dict.items()])
    fp_dict_1 = '; '.join([f"{k}: {v}" for k, v in fp_dict_1.items()])
    return file_path, fp_dict, fp_dict_1

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
    pd.options.mode.chained_assignment = None
    
    start = time.time()
    # FIRST BLOCK: Generate fingerprints for the COD
    entries = []
    cod_dir = '../cif'
    df = pd.read_csv('../molecule_libraries/COD/final_structures.csv')
    # df = df_origin[(df_origin['Polymetallic'] == False) & (df_origin['Connectivity'] == 'Complete') & (df_origin['Metal_Bonds'] == True) & (df_origin['Overlapping Atoms'] == False)]
    
    # Extract the entries from the 'File Name' of the csv dataframe
    # and add '../' to the beginning of the entry
    entries = [f"{entry}" for entry in df['File Name']]
   
    # For testing purposes, limit the number of entries to process
    # entries = entries[:10]
    print(f"\nRetrieved {len(entries)} entries from the COD. In {round((time.time() - start)/60, 2)} minutes.")
    print('Starting the fingerprint generation')
    print(f"Generating fingerprints for {len(entries)} entries with {cpu_count()} CPUs...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(get_fingerprints, entries)  
    print(f'Finished generating fingerprints for {len(results)} entries.')
    print(f'\nUpdating the dataframe...')
    # Correctly extracting fingerprints dictionaries
    fingerprints_dict = {}
    fingerprints_dict_1 = {}
    for entry, fp_dict, fp_dict_1 in results:
        fingerprints_dict[entry] = fp_dict
        fingerprints_dict_1[entry] = fp_dict_1
    
    # Add the fingerprints to the dataframe to the corresponding entry
    # df['Fingerprint'] = df['File Name'].map(fingerprints_dict)
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore", category=pd.core.common.SettingWithCopyWarning) 
    df.loc[:, 'Fingerprint_3D'] = df['File Name'].map(fingerprints_dict)
    df.loc[:, 'Fingerprint_4D'] = df['File Name'].map(fingerprints_dict_1)
        
    print('Fingerprints added to the dataframe.')
    print(f'\nSaving the dataframe to a CSV file...')
    # Save the dataframe to a CSV file
    df.to_csv('final_structures_with_fp.csv', index=False)
    print(f'Dataframe saved.\nTotal time elapsed for fingerprints calculation: {(time.time() - start)/60} minutes\n')
    
    # SECOND BLOCK: Compute similarities to different targets
    start2 = time.time()
    print('Starting the similarity calculation...')
    csd_reader = io.EntryReader('CSD')
    
    # POLYMERIZATION TARGET
    target_molecule_1 = csd_reader.entry('GUDQOL').molecule
    target_fp_3d_1, target_fp_4d_1 = hsr_fingerprint(target_molecule_1)
    
    # OLEFIN METATHESIS TARGET
    mol_reader1 = io.MoleculeReader('../targets/grubbs:olef_met/hb7792sup1.cif')
    target_molecule_2 = mol_reader1[0]
    target_fp_3d_2, target_fp_4d_2 = hsr_fingerprint(target_molecule_2) 
    
    # PNP LIGAND TARGET
    mol_reader2 = io.MoleculeReader('../targets/PNP_type/ScienceDirect_files_19Jul2024_13-16-49.811/Complex_2.final.cif')
    target_molecule_3 = mol_reader2[0]
    target_fp_3d_3, target_fp_4d_3 = hsr_fingerprint(target_molecule_3)
    
    # SALEN TYPE TARGET
    target_molecule_4 = csd_reader.entry('XIDTOW').molecule
    target_fp_3d_4, target_fp_4d_4 = hsr_fingerprint(target_molecule_4)
    
    taget_molecules = [target_molecule_1, target_molecule_2, target_molecule_3, target_molecule_4]
    target_fps_3d = [target_fp_3d_1, target_fp_3d_2, target_fp_3d_3, target_fp_3d_4] 
    target_fps_4d = [target_fp_4d_1, target_fp_4d_2, target_fp_4d_3, target_fp_4d_4]    
    results_files_3d = ['similarity_results_polymerization_3d.txt', 'similarity_results_grubbs_3d.txt', 'similarity_results_PNP_3d.txt', 'similarity_results_salen_3d.txt']
    results_files_4d = ['similarity_results_polymerization_4d.txt', 'similarity_results_grubbs_4d.txt', 'similarity_results_PNP_4d.txt', 'similarity_results_salen_4d.txt']
    # # couple fingerprints and results files
    # target_fp_results = list(zip(target_fps, results_files))
    
    
    df = pd.read_csv('final_structures_with_fp.csv', low_memory=False)
    # Get only the first 10 entries for testing purposes
    # df = df.head(10)
    manager = Manager()
    fingerprint_dict_3d = manager.dict(pd.Series(df['Fingerprint_3D'].values, index=df['File Name']).to_dict())
    fingerprint_dict_4d = manager.dict(pd.Series(df['Fingerprint_4D'].values, index=df['File Name']).to_dict())     
    entries = list(fingerprint_dict_3d.keys())
    
    
    for target_fp_3d, target_fp_4d, result_file_3d, result_file_4d in zip(target_fps_3d, target_fps_4d, results_files_3d, results_files_4d):
        print(f"Calculating similarity to target for {len(entries)} entries with {cpu_count()} CPUs...")
        with Pool(processes=cpu_count(), initializer=init_process, initargs=(target_fp_3d, target_fp_4d, fingerprint_dict_3d, fingerprint_dict_4d)) as pool:
            results = pool.map(calculate_similarity, entries)
        print(f'Finished calculating similarity for {len(results)} entries')
        
        # get results for 3D and 4D fingerprints in separate lists
        results_3d = [(result[0], result[1]) for result in results]
        results_4d = [(result[0], result[2]) for result in results]
        
        # Sort results by similarity score
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

    print(f'Total time for similarities calculation: {(time.time() - start2)/60:.2f} minutes')
    # index = 0
    # for target_fp, result_file in target_fp_results:
    #     index += 1
    #     print(f"Calculating similarity to target {index} for {len(entries)} entries with {cpu_count()} CPUs...")
    #     with Pool(processes=cpu_count(), initializer=init_process, initargs=(target_fp, fingerprint_dict)) as pool:
    #         results = pool.map(calculate_similarity, entries)
    #     print(f'Finished calculating similarity for {len(results)} entries\n')
    
    #     results = sorted(results, key=lambda x: max(x[1]) if x[1] else 0, reverse=True)
    
    #     print(f'Saving results to {result_file}...')
    #     with open(result_file, 'w') as file:
    #         for result in results:
    #             file.write(f"{result[0]}: {result[1]}\n")
    
    # print(f'Total time for similarities calculation: {(time.time() - start2)/60:.2f} minutes')

    # print(f'Total time elapsed: {(time.time() - start)/60:.2f} minutes')