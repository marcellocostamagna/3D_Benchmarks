
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from ccdc import io
from rdkit import DataStructs
import pandas as pd
import time
from ccdc import search


def get_fingerprint(entry):
    fps = []    
    # try:
    molecule = entry.molecule
    for component in molecule.components:
        # try:
            smiles = component.smiles
            rdkit_molecule = Chem.MolFromSmiles(smiles, sanitize=False)
            # fps.append(AllChem.RDKFingerprint(rdkit_molecule))
            # Does not work because it is slow and crashes for complex molecules
            
            # fps.append(AllChem.GetMorganFingerprintAsBitVect(rdkit_molecule, 2, nBits=1024))
            # Does not work becuse it need kekulization, sanitation and valence model
            
            # fps.append(Chem.GetHashedTopologicalTorsionFingerprintAsBitVect(rdkit_molecule, 1024))
            # does not work
            
            fp = MACCSkeys.GenMACCSKeys(rdkit_molecule)
            fps.append(fp)
            # fps.append(Chem.MACCSkeys.GenMACCSKeys(rdkit_molecule))
            
        # except:
        #     # print(f'Error in {entry.identifier} in {component.identifier}')
        #     fps.append(None)
    return entry.identifier , fps
    # except:
    #     return entry.identifier, None          


if __name__ == '__main__':
    print(f'Starting retrieving entries...')
    start = time.time()
    n_entries = 10
    csd_reader = io.EntryReader('CSD') 
    df = pd.read_csv('CSD_similarity/final_structures_with_fp.csv')
    
    # Pick only the entries with 'Hapticity' True
    entries = df[df['Hapticity'] == True]['Identifier'].tolist()[:n_entries]
    
    # entries = ['GUDQOL', 'TITTUO', 'NIVHEJ', 'XIDTOW']
    entries = [csd_reader.entry(identifier) for identifier in entries]
    #print the smiles of the first 5 entries one per line
    for entry in entries[:5]:
        print(f'{entry.identifier}: {entry.molecule.smiles}')
    # print(f'Entries retrieved: {len(entries)} in {(time.time()-start)/60} minutes')
    
    # molecule1 = entries[0].molecule
    # similarity_search1 = search.SimilaritySearch(molecule1)
    # fingerprint_builder1 = similarity_search1._fp
    # fingerprint1 = fingerprint_builder1.similarity_fingerprint(molecule1._molecule)

    # molecule2 = entries[1].molecule
    # similarity_search2 = search.SimilaritySearch(molecule2)
    # fingerprint_builder2 = similarity_search2._fp
    # fingerprint2 = fingerprint_builder2.similarity_fingerprint(molecule2._molecule)

    # # Compute the similarity between the fingerprints
    
    # dice_similarity_score = fingerprint1.dice(fingerprint2)

    # tanimoto_similarity_score = fingerprint1.tanimoto(fingerprint2)

    # print(f"Dice similarity between molecule1 and molecule2: {dice_similarity_score}")
    # print(f"Tanimoto similarity between molecule1 and molecule2: {tanimoto_similarity_score}")
    
    
    # start2 = time.time()
    # # Get the fingerprints for the first 1000 entries
    # results = [get_fingerprint(entry) for entry in entries]
    # # make a dictionary with the entry and the fingerprint
    # entry_fingerprint = {entry: fingerprints for entry, fingerprints in results if fingerprints is not None} 

    # # get the total number of fingerprints
    # n_fingerprints = sum([len(fingerprints) for fingerprints in entry_fingerprint.values()])
    # print(f'Time taken to compute {n_fingerprints} fingerprints: {(time.time()-start2)/60} minutes')

    # # Compute the similarity of the first against all the other fingerprints of every entry
    # start3 = time.time()
    # similarities = {}
    # # get the fingerrint of the first entry
    # fisrt_key = list(entry_fingerprint.keys())[0]
    # target_fingerprint = entry_fingerprint[fisrt_key][0]
    # for entry_id, fingerprints in entry_fingerprint.items():
    #     entry_similarities = []
    #     for fingerprint in fingerprints:
    #         if fingerprint is not None:
    #             similarity = DataStructs.TanimotoSimilarity(target_fingerprint, fingerprint)
    #             entry_similarities.append(similarity)
    #         else:
    #             entry_similarities.append(None)
    #     similarities[entry_id] = entry_similarities

    # # print each entry and its similarity on a new line
    # count = 0
    # for entry_id, entry_similarities in similarities.items():
    #     count += 1
    #     print(f'{entry_id}: {entry_similarities}')
    #     if count == 5:
    #         break
        
    # # get the total number of similarities
    # n_similarities = sum([len(similarity) for similarity in similarities.values()])
    
    # print(f'Time taken to compute {n_similarities} similarities: {(time.time()-start3)/60} minutes')

# molecules = [rdkit_molecule_1, rdkit_molecule_2, rdkit_molecule_3, rdkit_molecule_4]
# # fpgen = AllChem.GetRDKitFPGenerator()
# fps = [AllChem.RDKFingerprint(x) for x in molecules] 

# # print all the similarities

# for i in range(4):
#     for j in range(i+1, 4):
#         similarity = DataStructs.TanimotoSimilarity(fps[i], fps[j])
#         print(f'Similarity between molecule {i+1} and molecule {j+1}: {similarity}')



# # POLYMERIZATION TARGET
# target_molecule_1 = csd_reader.entry('GUDQOL').molecule
# smiles = target_molecule_1.smiles
# rdkit_molecule_1 = Chem.MolFromSmiles(smiles, sanitize=False)

# # OLEFIN METATHESIS TARGET
# target_molecule_2 = csd_reader.entry('TITTUO').molecule
# smiles = target_molecule_2.smiles
# rdkit_molecule_2 = Chem.MolFromSmiles(smiles, sanitize=False)

# # PNP LIGAND TARGET
# target_molecule_3 = csd_reader.entry('NIVHEJ').molecule
# smiles = target_molecule_3.smiles
# rdkit_molecule_3 = Chem.MolFromSmiles(smiles, sanitize=False)

# # # SALEN TYPE TARGET
# target_molecule_4 = csd_reader.entry('XIDTOW').molecule
# smiles = target_molecule_4.smiles
# rdkit_molecule_4 = Chem.MolFromSmiles(smiles, sanitize=False)
