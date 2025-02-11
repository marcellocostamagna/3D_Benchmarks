from ccdc import io

# ### RETURNS THE BOND TYPES OF A MOLECULE ###
# csd_reader = io.EntryReader('CSD')
# # entry = csd_reader.entry('DOYQOX')
# # molecule = entry.molecule
# # bond_types = {}
# # for bond in molecule.bonds:
# #     bond_type = bond.bond_type  # Get the numerical value of the bond type
# #     text = bond.bond_type._bond_type.text()   # Get the text representation of the bond type
# #     value = str(bond.bond_type._bond_type.value()) # Get the value of the bond type
    
# #     key = f'{text}, ({value})'
    
# #     if key in bond_types:
# #         bond_types[key] += 1
# #     else:
# #         bond_types[key] = 1
# # print(f'The molecule has {len(molecule.bonds)} bonds')
# # print(f'The bond types are: {bond_types}')

# # ITERATE OVER THE MOLECULES AND COLLECT 20 MOLECULES
# # which are organic and have 3d coordinates
# organic_3d_mols = []
# for entry in csd_reader.entries():
#     molecule = entry.molecule
#     crystal = entry.crystal
#     if molecule.is_organic and molecule.is_3d and not crystal.has_disorder and len(molecule.components) == 1:
#         organic_3d_mols.append(molecule)
#     if len(organic_3d_mols) == 20:
#         break
    
# print(f'Collected {len(organic_3d_mols)} organic molecules with 3D coordinates')
# # print the identifiers of the collected molecules
# for molecule in organic_3d_mols:
#     print(molecule.identifier)

# reader_formats = sorted(io.MoleculeReader.known_formats.keys())
# writer_formats = sorted(io.MoleculeWriter.known_formats.keys())

# print('Reader formats:\n')
# print('\n'.join(reader_formats))

# print('\nWriter formats:\n')
# print('\n'.join(writer_formats))

# how to check if an entry has smiles
csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry('ABASAB')
# for component in entry.molecule.components:
#     if component.smiles:
#         print(component.smiles)
#         print(f'The component has smiles: {component.smiles}')
#     else:
#         print(component.smiles)
#         print('The component does not have smiles')
heavy_comp = entry.molecule.heaviest_component
print(heavy_comp.smiles)