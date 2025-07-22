from ccdc import io

entries = ['ADIYES', 'BIWHEY', 'CAFWEM']

csd_reader = io.EntryReader('CSD')

for entry in entries:
    print(entry)
    csd_entry = csd_reader.entry(entry)
    for i, component in enumerate(csd_entry.molecule.components):
        print(component.identifier)
        filename = f'{entry}_component{i+1}.sdf'
        with io.MoleculeWriter(filename) as writer:
            writer.write(component)
