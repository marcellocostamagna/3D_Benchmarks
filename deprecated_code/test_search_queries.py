# Script to test the correct functionality of all the queries 
# defined in pattern_queries.py

from pattern_queries import *
import ccdc.io
import ccdc.search
from pattern_queries import METALS, NON_METALS
from pattern_queries import *
import time

# Search for one query
def search_molecules(queries, identifiers):
    '''
    Searches for a list of queries in a subset of the CSD defined by 
    the identifiers of the entries.
    The search is operated only on the heaviest component of the molecule.
    
    Returns the results of the search as a list of the hits' identifiers.
    '''
    results_ids = set()
    entry_reader = ccdc.io.EntryReader('CSD')   
    entries_dict = {id: entry_reader.entry(id) for id in identifiers}
    heaviest_molecules = {
        id: entry.molecule.heaviest_component for id, entry in entries_dict.items()
    }
    # if queries is a single query, convert it to a list
    if not isinstance(queries, list):
        queries = [queries]
    for query in queries:
        search = ccdc.search.SubstructureSearch()
        search.add_substructure(query)
        
        results = search.search(list(heaviest_molecules.values()))
        for hit in results:
            matched_molecule = hit.molecule

            for entry_id, molecule in heaviest_molecules.items():
                if molecule == matched_molecule:
                    results_ids.add(entry_id)
                    break
                
    return list(results_ids)

def search_entries(queries, identifiers):
    '''
    Searches for a list of queries in a subset of the CSD defined by 
    the identifiers of the entries.
    The search is operated only on the heaviest component of the molecule.
    
    Returns the results of the search as a list of the hits' identifiers.
    '''
    results_ids = set()
    for query in queries:
        search = ccdc.search.SubstructureSearch()
        search.add_substructure(query)
        results = search.search(identifiers)
        results_ids.update(hit.identifier for hit in results) 
    return list(results_ids)    
start = time.time()   

# TEST BRIDGED COMPLEXES and CHELATING LIGANDS
mono_bridged_query = mono_bridged_complexes(METALS, NON_METALS, 4)
multi_bridged_query = multi_bridged_complexes(METALS, NON_METALS, 4)
chelating_ligands_query = chelating_ligands(METALS)    

identifiers = ['ABABOY', # mono bridged
               'ABAROO', # multi bridged
               'ABEQUV01', # multi bridged
               'ABEVIN', # mono bridged
               'ABANAV', # multi bridged
               'ABATIH', # multi bridged
               'ABAWUV', # multi bridged
               'AQUXUI', # mono, multi bridged
               'COPMOJ10', # mono bridged
               'KAZSIN'] # mono bridged (The only one that has no chelating ligands)
    
results1 = search_molecules(mono_bridged_query, identifiers)
results2 = search_molecules(multi_bridged_query, identifiers)
results3 = search_molecules(chelating_ligands_query, identifiers)

print(f'Mono bridged complexes: {results1}')
print(f'Multi bridged complexes: {results2}')
print(f'Chelating ligands: {results3}')

# TEST LINEAR LIGANDS , METAL ALLENES
linear_ligands_query = linear_ligands(METALS, NON_METALS)
metal_allenes_query = metal_allenes(METALS, NON_METALS)

identifiers = ['KAZSIN', # 
               'ABICET', # metal allene
               'AKEHIH', # metal allene
               'ABABOY',
               'COPMOJ10',
               'XIDTOW' ] # The only one that has no linear ligands
                
results1 = search_molecules(linear_ligands_query, identifiers)
results2 = search_molecules(metal_allenes_query, identifiers)

print(f'Linear ligands: {results1}')
print(f'Metal allenes: {results2}')

# TEST DATIVE BONDS and SIGMA BONDS
dative_bonds_query = dative_bonds(METALS, NON_METALS)
dative_bonds_in_rings_query = dative_in_rings(METALS, NON_METALS, 6)
dative_bonds_in_conjugated_chains_query = dative_in_conjugated_chains(METALS, NON_METALS)
dative_bonds_in_conjugated_rings_query = dative_in_conjugated_rings(METALS, NON_METALS)
cyclic_multihapto_ligands_query = cyclic_multihapto_ligands(METALS, NON_METALS)       
sigma_bonds_query = sigma_bonds(METALS)
sigma_bonds_in_rings_query = sigma_bonds_in_rings(METALS, NON_METALS, 6)

# All contain dative bonds (The matches returned are correct)
identifiers = ['ABABOY', 
               'ABAFIT',
               'ABAZEK',
               'ABESEF',
               'ABETIK',
               'AXUTIZ',
               'AYIYIQ',
               'FOZNAJ',
               'LISRIO',
               'ABETUW',
               'AFAGAP',
               'AFIJEE',
               'EWIYEP',
               'UJACAK',
               'YIDBUJ',
               'DECXOY',
               'ACAZFE',
               'HESMUQ01',
               'YIDBOD',
               ]

results1 = search_molecules(dative_bonds_query, identifiers)
results2 = search_molecules(dative_bonds_in_rings_query, identifiers)
results3 = search_molecules(dative_bonds_in_conjugated_chains_query, identifiers)
results4 = search_molecules(dative_bonds_in_conjugated_rings_query, identifiers)
results5 = search_molecules(cyclic_multihapto_ligands_query, identifiers)
results6 = search_molecules(sigma_bonds_query, identifiers)
results7 = search_molecules(sigma_bonds_in_rings_query, identifiers)

print(f'Dative bonds: {results1}')
print(f'Dative bonds in rings: {results2}')
print(f'Dative bonds in conjugated chains: {results3}')
print(f'Dative bonds in conjugated rings: {results4}')
print(f'Cyclic multihapto ligands: {results5}')
print(f'Sigma bonds: {results6}')
print(f'Sigma bonds in rings: {results7}')
               
# TEST QUERIES FOR GENERAL COMPOUNDS (NON METALS)
allenes_query = allenes(NON_METALS)
saturated_chains_query = saturated_chains(NON_METALS)
unsaturated_chains_query = unsaturated_chains(NON_METALS)
unsaturated_conjugated_chains_query = unsaturated_conjugated_chains(NON_METALS)
saturated_rings_query = saturated_rings(NON_METALS, 8)
unsaturated_rings_query = unsaturated_rings(NON_METALS, 8)
unsaturated_conjugated_rings_query = unsaturated_conjugated_rings(NON_METALS, 5)
aromatic_rings_query = aromatic_rings(NON_METALS)
               
identifiers = ['AACANI10',
               'AACMHX10',
               'AAGGAG10',
               'AAPUNI',
               'ABABEL',
               'ABAZOW',
               'ABOPAL',
               'ABOSES',
               'AABHTZ',
               'ABAKAR',
               'ARESAR',
               'AACFAZ10',
               'ABATAB',
               'AKAPUZ',
               ]

results1 = search_molecules(allenes_query, identifiers)
results2 = search_molecules(saturated_chains_query, identifiers)
results3 = search_molecules(unsaturated_chains_query, identifiers)
results4 = search_molecules(unsaturated_conjugated_chains_query, identifiers)
results5 = search_molecules(saturated_rings_query, identifiers)
results6 = search_molecules(unsaturated_rings_query, identifiers)
results7 = search_molecules(unsaturated_conjugated_rings_query, identifiers)
results8 = search_molecules(aromatic_rings_query, identifiers)

print(f'Allenes: {results1}')
print(f'Saturated chains: {results2}')
print(f'Unsaturated chains: {results3}')
print(f'Unsaturated conjugated chains: {results4}')
print(f'Saturated rings: {results5}')
print(f'Unsaturated rings: {results6}')
print(f'Unsaturated conjugated rings: {results7}')
print(f'Aromatic rings: {results8}')

end = time.time()
print(f"Time elapsed: {end-start:.2f} s")



