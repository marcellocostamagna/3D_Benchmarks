# print the elements of the heaviest component

import ccdc.io

csd_reader = ccdc.io.EntryReader('CSD')
entry = csd_reader.entry('ABEYUC')
molecule = entry.molecule
heavy_comp = molecule.heaviest_component

atoms = []
for atom in heavy_comp.atoms:
    atoms.append(atom.atomic_symbol)

bonds = {}
for bond in heavy_comp.bonds:
    # key: the couple of atomic symbols
    key = f'{bond.atoms[0].atomic_symbol}-{bond.atoms[1].atomic_symbol}'
    # value: the bond type
    value = bond.bond_type
    value = bond.is_cyclic
    bonds[key] = value

print(atoms)
print(bonds)