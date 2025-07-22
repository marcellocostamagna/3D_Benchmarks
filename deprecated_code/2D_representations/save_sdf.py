from ccdc import io


if __name__ == '__main__':
    csd_reader = io.EntryReader('CSD')
    entry = csd_reader.entry('ABOLOU')
    molecule = entry.molecule
    
    file_name = 'ABOLOU.sdf'
    
    with io.MoleculeWriter(file_name) as mol_writer:
        mol_writer.write(molecule)
    