from ccdc import io
from hsr import fingerprint as fp
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import time

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


def get_fingerprints(entry):
    fp_dict = {}
    fp_dict_1 = {}
    csd_reader = io.EntryReader('CSD')
    molecule = csd_reader.entry(entry).molecule
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
    return entry, fp_dict, fp_dict_1

if __name__ == '__main__':
    start = time.time()
    print('Starting the fingerprint generation')
    # Load the dataframe
    df = pd.read_csv('../molecule_libraries/CSD/final_structures.csv')
    entries = df['Identifier'].tolist()
    print(f"Generating fingerprints for {len(entries)} entries with {cpu_count()} CPUs...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(get_fingerprints, entries)  
    print(f'Finished generating fingerprints for {len(results)} entries. Updating the dataframe...')
    # Correctly extracting fingerprints dictionaries
    fingerprints_dict = {}
    fingerprints_dict_1 = {}
    for entry, fp_dict, fp_dict_1 in results:
        fingerprints_dict[entry] = fp_dict
        fingerprints_dict_1[entry] = fp_dict_1
    # Add the fingerprints to the dataframe to the corresponding entry
    df['Fingerprint_3D'] = df['Identifier'].map(fingerprints_dict)
    df['Fingerprint_4D'] = df['Identifier'].map(fingerprints_dict_1)
    print('Fingerprints added to the dataframe. Saving the dataframe to a CSV file...')
    # Save the dataframe to a CSV file
    df.to_csv('final_structures_with_fp.csv', index=False)
    print(f'Dataframe saved.\n Total time elapsed: {(time.time() - start)/60} minutes')
        