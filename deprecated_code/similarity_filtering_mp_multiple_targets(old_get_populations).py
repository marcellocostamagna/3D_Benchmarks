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
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

# Function to generate a fingerprint from a molecule
def get_fingerprint(ccdcmol):
    atom_array = get_array_from_ccdcmol(ccdcmol)
    return hsr.fingerprint.generate_fingerprint_from_data(atom_array)

# Function to determine the component of interest
def component_of_interest(molecule):
    components = molecule.components
    properties = []
    
    for comp in components:
        is_organometallic = comp.is_organometallic
        molecular_weight = sum(atom.atomic_weight for atom in comp.atoms)
        atom_count = len(comp.atoms)
        
        properties.append({
            "component": comp,
            "is_organometallic": is_organometallic,
            "molecular_weight": molecular_weight,
            "atom_count": atom_count
        })
    
    heaviest_component = max(properties, key=lambda x: x["molecular_weight"])
    most_atoms_component = max(properties, key=lambda x: x["atom_count"])
    
    for prop in properties:
        criteria_met = sum([
            prop["is_organometallic"],
            prop["component"] == heaviest_component["component"],
            prop["component"] == most_atoms_component["component"]
        ])
        if criteria_met >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None

# Function to compute similarity for a chunk of entries
def process_chunk(chunk_args):
    target_fp, entries_chunk = chunk_args
    results = {}
    entry_reader = ccdc.io.EntryReader('CSD')

    for entry in entries_chunk:
        try:
            molecule = component_of_interest(entry_reader.entry(entry).molecule)
            entry_fp = get_fingerprint(molecule)
            similarity = hsr.similarity.compute_similarity_score(target_fp, entry_fp)
            mol_smiles = molecule.smiles
            if similarity < 0.3:
                results[entry] = similarity, mol_smiles
        except Exception as e:
            print(f"Error processing entry {entry}: {e}")
    return results

if __name__ == "__main__":
    global_start_time = time.time()
    
    # Load the dataframe of viable structures
    csv_file = './viable_structures.csv'
    df = pd.read_csv(csv_file)
    # entries = df['Identifier'].tolist()[:1000] # For testing purposes
    entries = df['Identifier'].tolist()
    # Define target molecules
    target_entries = ['ABAHIW', 'ABAKIZ', 'ABADOX',
                      'ABABIP', 'GASQOK', 'ABEKIE',
                      'NIWPUE01', 'ABEKIF', 'APUFEX',
                      'ABEHAU', 'TITTUO', 'EGEYOG',
                      'ABOBUP', 'XIDTOW', 'ACNCOB10',
                      'TACXUQ','ACAZFE', 'NIVHEJ',
                      'ADUPAS','DAJLAC', 'OFOWIS',
                      'CATSUL','HESMUQ01', 'GUDQOL',
                      'ABEVAG', 'AKOQOH', 'ADARUT',
                      'AFECIA', 'ACOVUL', 'AFIXEV'] 
    
    
    entry_reader = ccdc.io.EntryReader('CSD')
    
    for i, target_entry in enumerate(target_entries, start=1):
        print(f"Processing target: {target_entry}")
        target_molecule = component_of_interest(entry_reader.entry(target_entry).molecule)
        target_fp = get_fingerprint(target_molecule)
        
        num_chunks = 8  # Adjust as needed
        chunk_size = len(entries) // num_chunks + 1
        entry_chunks = [entries[i:i + chunk_size] for i in range(0, len(entries), chunk_size)]

        print(f"Calculating similarity scores for {len(entries)} entries using {num_chunks} processes...")
        start_time = time.time()
        
        with Pool(processes=num_chunks) as pool:
            chunk_args = [(target_fp, chunk) for chunk in entry_chunks]
            results_list = pool.map(process_chunk, chunk_args)
        
        mols_per_target = {}
        for partial_results in results_list:
            mols_per_target.update(partial_results)
        
        final_mols = list(mols_per_target.keys())
        print(f'Number of final molecules for {target_entry}: {len(final_mols)}')
        
        output_file = f'./{i}_{target_entry}_init_pop.txt'
        with open(output_file, 'w') as f:
            for mol in final_mols:
                similarity, smiles = mols_per_target[mol]
                f.write(f"{mol} {smiles} {similarity:.4f}\n") 
        
        print(f'Results saved to {output_file}')
        print(f'Time elapsed for {target_entry}: {time.time() - start_time:.2f} seconds')
    
    print(f'Total time elapsed: {time.time() - global_start_time:.2f} seconds')
