from ccdc import io
from ccdc import conformer
from ccdc.molecule import Molecule
import pandas as pd

df = pd.read_csv('./molecule_libraries/CSD/final_structures_with_smiles.csv')

smiles = df['SMILES'].tolist()
# Remove NaN values
smiles = [smile for smile in smiles if type(smile) == str]

# get first 5 smiles
smiles = smiles[:5]
conformer_generator = conformer.ConformerGenerator()
mols = []
for smile in smiles:
    mol = Molecule.from_string(smile)
    conf = conformer_generator.generate(mol, hydrogens=False)
    mol_3d = conf.hits[0].molecule
    mols.append(mol_3d)
    
# Save each conformer in a sdf file with a field for the corresponding smiles string
for i, smile, mol in zip(range(len(smiles)), smiles, mols):
    filename = f'./molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_to_3D_with_ccdc/conformer_{i}.sdf'
    with io.MoleculeWriter(filename) as w:
        w.write(mol)
   
    



     
    
    