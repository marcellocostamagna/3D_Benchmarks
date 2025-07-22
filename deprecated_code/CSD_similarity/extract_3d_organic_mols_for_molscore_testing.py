from ccdc import io
import os

csd_reader = io.EntryReader('CSD')
organic_3d_mols = []
for entry in csd_reader.entries():
    molecule = entry.molecule
    crystal = entry.crystal
    if molecule.is_organic and molecule.is_3d and not crystal.has_disorder and len(molecule.components) == 1:
        organic_3d_mols.append(molecule)
    if len(organic_3d_mols) == 20:
        break
    
print(f'Collected {len(organic_3d_mols)} organic molecules with 3D coordinates')
    
# cretate a directory for storing molecules files (sdf format)
output_dir = './CSD_similarity/organic_molecules_for_benchmark_test'
os.makedirs(output_dir, exist_ok=True)

for mol in organic_3d_mols:
    entry = mol.identifier
    # get the molecule with the highest similarity score
    csd_reader = io.EntryReader('CSD')
    molecule = csd_reader.entry(entry).molecule.heaviest_component
    # save the molecule to a file
    with io.MoleculeWriter(f'{output_dir}/{entry}.sdf') as w:
        w.write(molecule)