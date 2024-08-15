from ccdc import io


csd_reader = io.EntryReader('CSD')


# POLYMERIZATION TARGET
target_molecule_1 = csd_reader.entry('GUDQOL').molecule

# OLEFIN METATHESIS TARGET
target_molecule_2 = csd_reader.entry('TITTUO').molecule

# PNP LIGAND TARGET
target_molecule_3 = csd_reader.entry('NIVHEJ').molecule

# SALEN TYPE TARGET
target_molecule_4 = csd_reader.entry('XIDTOW').molecule

# ION PAIR TARGET
file = 'targets/ion_pairs/ion_pair_2.sdf'
target_molecule_5 = io.MoleculeReader(file)[0]

# HOST-GUEST TARGET
file = 'targets/host_guest/1/3185287/ja052912c_si_002.cif'
target_molecule_6 = io.MoleculeReader(file)[0]

targets = [target_molecule_1, target_molecule_2, target_molecule_3, target_molecule_4, target_molecule_5, target_molecule_6]
files = ['polymerization.sdf', 'olefin_metathesis.sdf', 'PNP_ligand.sdf', 'salen.sdf', 'ion_pair.sdf', 'host_guest.sdf']

for target, file in zip(targets, files):
    with io.MoleculeWriter(file) as w:
        w.write(target)
    