# print the elements of the heaviest component

import ccdc.io

csd_reader = ccdc.io.EntryReader('CSD')
entry = csd_reader.entry('ACUKAM')
molecule = entry.molecule
heavy_comp = molecule.heaviest_component

atoms = []
for atom in heavy_comp.atoms:
    atoms.append(atom.atomic_symbol)

print(f'Number of atoms in the heaviest component: {len(atoms)}')
# print MW
print(f'Molecular weight: {heavy_comp.molecular_weight}')

for component in molecule.components:
    print(f'Number of atoms in the component: {len(component.atoms)}')
    print(f'Molecular weight: {component.molecular_weight}')
    for atom in component.atoms:
        print(atom.atomic_symbol)