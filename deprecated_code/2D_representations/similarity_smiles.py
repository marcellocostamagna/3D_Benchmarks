from ccdc import io
import pandas as pd
from multiprocessing import Pool, cpu_count, Manager
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity


    
def calculate_similarity(identifier):
    """Calculates similarities of the fingerprints against the target fingerprint."""
   
    csd_reader = io.EntryReader('CSD')
    molecule = csd_reader.entry(identifier).molecule
    similarities = []    
    
    for component in molecule.components:
        try:
            smiles = component.smiles
            rdkit_molecule = Chem.MolFromSmiles(smiles, sanitize=False)
            query_fp = AllChem.RDKFingerprint(rdkit_molecule)
            similarities.append(TanimotoSimilarity(trgt_fp, query_fp))
        except:
            similarities.append(None)
          
    return identifier, similarities

def init_process(target_fp):
    """Initializer function to setup global variables."""
    global trgt_fp
    trgt_fp = target_fp

if __name__ == '__main__':
    start = time.time()
    csd_reader = io.EntryReader('CSD')
    
    # POLYMERIZATION TARGET
    target_molecule_1 = csd_reader.entry('GUDQOL').molecule
    smiles = target_molecule_1.smiles
    rdkit_molecule_1 = Chem.MolFromSmiles(smiles, sanitize=False)
  
    # OLEFIN METATHESIS TARGET
    target_molecule_2 = csd_reader.entry('TITTUO').molecule
    smiles = target_molecule_2.smiles
    rdkit_molecule_2 = Chem.MolFromSmiles(smiles, sanitize=False)

    # PNP LIGAND TARGET
    # mol_reader = io.MoleculeReader('targets/PNP_type/ScienceDirect_files_19Jul2024_13-16-49.811/Complex_2.final.cif')
    # print(os.getcwd())  
    target_molecule_3 = csd_reader.entry('NIVHEJ').molecule
    smiles = target_molecule_3.smiles
    rdkit_molecule_3 = Chem.MolFromSmiles(smiles, sanitize=False)
  
    # # SALEN TYPE TARGET
    target_molecule_4 = csd_reader.entry('XIDTOW').molecule
    smiles = target_molecule_4.smiles
    rdkit_molecule_4 = Chem.MolFromSmiles(smiles, sanitize=False)
   

    target_molecules = [rdkit_molecule_1, rdkit_molecule_2, rdkit_molecule_3, rdkit_molecule_4]
    target_fps = [AllChem.RDKFingerprint(x) for x in target_molecules] 
    results_files_smiles = ['similarity_results_polymerization_smiles.txt', 'similarity_results_grubbs_smiles.txt', 'similarity_results_PNP_smiles.txt', 'similarity_results_salen_smiles.txt']
    df = pd.read_csv('final_structures_with_fp.csv')
    manager = Manager()
    entries = df['Identifier'].tolist()
    
    for target_fp, result_file in zip(target_fps, results_files_smiles):
        print(f"Calculating similarity for {len(entries)} entries with {cpu_count()} CPUs...")
        with Pool(processes=cpu_count(), initializer=init_process, initargs=(target_fp, )) as pool:
            results = pool.map(calculate_similarity, entries)
        print(f'Finished calculating similarity for {len(results)} entries')
        
        # Sort the results by similarity score
        results = sorted(results, key=lambda x: x[1], reverse=True)
       
        # Save the results in the two files
        print(f'Saving smiles similarity results to {result_file}...')
        with open(result_file, 'w') as file:
            for result in results:
                file.write(f"{result[0]}: {result[1]}\n")
                
    print(f'Total time elapsed: {(time.time() - start)/60:.2f} minutes')
