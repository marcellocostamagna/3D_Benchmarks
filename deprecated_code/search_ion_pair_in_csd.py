# Script to search for ion pairs in the CSD

from ccdc import io
from ccdc import search
from ccdc.search import SimilaritySearch, SubstructureSearch, MoleculeSubstructure


ion_pair_file = 'ion_pair_2.sdf'

# test (targets/grubbs_olef_met/)
ion_pair_file = 'hb7792sup1.cif'

mol_reader = io.MoleculeReader(ion_pair_file)
mol = mol_reader[0]

with io.MoleculeWriter('mol.sdf') as mol_writer:
    mol_writer.write(mol)

mol_reader = io.MoleculeReader('mol.sdf')
mol = mol_reader[0]

print(len(mol.atoms))
print(len(mol.bonds))

# csd_reader = io.EntryReader('CSD')
# mol = csd_reader.entry('TITTUO').molecule

sim_search = SimilaritySearch(mol)
hits = sim_search.search()

print(f'Number of hits: {len(hits)}')

# Get the identifiers of the hits
identifiers = [hit.identifier for hit in hits]

print(identifiers[:5])

# Other method

mol_substructure = MoleculeSubstructure(mol)
substructure_search = SubstructureSearch()
sub_id = substructure_search.add_substructure(mol_substructure)
hits = substructure_search.search()

print(f'Number of hits: {len(hits)}')

# Get the identifiers of the hits
identifiers = [hit.identifier for hit in hits]
print(identifiers[:5])