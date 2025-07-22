# Script to convert a SMILES string to a SMARTS string

from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_smarts(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES string: {smiles}')
    AllChem.Compute2DCoords(mol)
    return Chem.MolToSmarts(mol)

if __name__ == '__main__':
    # GRUBBS #
    # persistent carbene
    # smiles = 'Cc1cc(C)c(c(C)c1)N2CCN([C]2)c3c(C)cc(C)cc3C'
    # general NHC
    # smiles = 'N(C)2CCN(C)(C2(C))'
    
    # PNP #
    # Unbonded PNP
    # smiles = 'C1=CC=C2C=C3C=CC=C(CP(C)(C)(C))C3=NC2=C(CP(C)(C)(C))1'
    
    # P Bounded PNP
    smiles = 'C1=CC=C2C=C3C=CC=C4C3=NC2=C1CPBPC4'
    
    # PNP skeleton
    smiles = 'C1C=NC=CCPBPC1'
    
    # POLYMERIZATION #
    # 1. Zr with two Cp ligands and two Cl ligands
    # smiles = '[CH]1C=CC=C1.[CH]1C=CC=C1.[Zr](Cl)(Cl)'
    
    # 2. ligand unbonded
    smiles = 'C1C=CC=C1[Si][Si]C1C=CC=C1'
    
    # SALEN TYPE #
    # 1. Cu with two O and two N ligands
    
    # 2. Salen specific ligand
    smiles = 'c1cccc2c1c3c(O[Cu]45N(CCCN4(C(C)c6c(cccc6)O5))C3)cc2'
    
    # 3. General Salen ligand
    smiles = 'O1[Cu]23N(CCC1)(CCN2(CCCO3))'
  
    print(smiles_to_smarts(smiles))
   
