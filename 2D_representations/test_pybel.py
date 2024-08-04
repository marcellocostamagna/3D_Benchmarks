from ccdc import io
import os 
from openbabel import pybel
import json

csd_reader = io.EntryReader('CSD')
# Get the fist 10 entries from the CSD
entries = ['AACANI10', 'AACRUB', 'ABAHIV', 'ABANID']
smiles = {}
for entry in entries:
    print(f"\nEntry: {entry}")
    csd_reader = io.EntryReader('CSD')
    molecule = csd_reader.entry(entry).molecule
    temp_file = f"{entry}.sdf"
    smiles[entry] = []
    for component in molecule.components:
        molecule = component
        with io.MoleculeWriter(temp_file) as writer:
            writer.write(molecule)
        ob_molecule = next(pybel.readfile('sdf', temp_file))
        generated_smiles = ob_molecule.write("smi").strip().split()[0]
        smiles[entry].append(generated_smiles)
        os.remove(temp_file)
        print(generated_smiles)

# Save the smiles dictionary as a json file
with open('entry_smiles_test_1.json', 'w') as f:
    json.dump(smiles, f)