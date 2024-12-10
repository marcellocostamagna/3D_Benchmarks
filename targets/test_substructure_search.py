import ccdc.search
import time
import random
import ccdc.io
from pattern_queries import *
from pattern_queries import METALS, NON_METALS




def get_random_identifiers(number, space):
    csd_reader = ccdc.io.EntryReader('CSD')
    count = 0
    ids = []
    for entry in csd_reader.entries():
        count += 1
        ids.append(entry.molecule.identifier)
        if count == space:
            break
    # Randomly select n entries from the database
    ids = random.sample(ids, number)
    return ids

def search(query, subset_identifiers):
    search = ccdc.search.SubstructureSearch()
    search.add_substructure(query)
    results = search.search(subset_identifiers)
    return results
    

start = time.time()

# SIMULATE A SEARCH OF A SUBSET OF THE CSD
subset_identifiers = get_random_identifiers(100, 100)  
# subset_identifiers = ['YIDBUJ', 'UJACAK', 'HANWIC', 'ABEDES']  

# DFEINITION OF THE QUERY SUBSTRUCTURE
# Chooose a query from the pattern_queries.py
query = metal_in_ring(METALS)

# SEARCH 
# Perform the substructure search
results = search(query, subset_identifiers)

end = time.time()

# POST PROCESSING
# print time in seconds with two decimanl places
print(f"Time elapsed: {end-start:.2f} s")
# # Collect only the uniqie set of identifiers
results = set(hit.identifier for hit in results)
print(f"Number of hits: {len(results)}")
non_hits = [id for id in subset_identifiers if id not in results]
print(f"Number of non hits: {len(non_hits)}")

# Save hits identifiers in a txt file
with open('hits.txt', 'w') as f:
    for hit in results:
        f.write(f'{hit}\n')  
# Save non hits identifiers in a txt file
with open('non_hits.txt', 'w') as f:
    for non_hit in non_hits:
        f.write(f'{non_hit}\n')
        
        
# # # METAL IN RING QUERY
# def metal_in_ring(metals):
#     # Create a new query substructure
#     metal_in_ring_query = ccdc.search.QuerySubstructure()
#     # Add a metal query atom with a cyclic constraint
#     metal = ccdc.search.QueryAtom(metals)
#     metal.cyclic = True  # Add cyclic constraint
#     # Add the QueryAtom to the QuerySubstructure
#     metal_in_ring_query.add_atom(metal)
#     return metal_in_ring_query


# # # TWO METALS IN A RING (Iterative search)
# def two_metals_in_a_ring(metals, ligands, num_intermediates):
#     # Create a new query substructure
#     query = ccdc.search.QuerySubstructure()
#     # Add the first metal atom
#     metal1 = query.add_atom(ccdc.search.QueryAtom(metals))
#     metal1.cyclic = True 
#     # Add intermediate ligand atoms
#     previous_atom = metal1
#     for _ in range(num_intermediates):
#         ligand = query.add_atom(ccdc.search.QueryAtom(ligands))
#         ligand.cyclic = True  
#         query.add_bond('Any', previous_atom, ligand) 
#         previous_atom = ligand  
#     # Add the second metal atom
#     metal2 = query.add_atom(ccdc.search.QueryAtom(metals))
#     metal2.cyclic = True 
#     query.add_bond('Any', previous_atom, metal2)  
#     return query

# # DATIVE BONDS (not in rings)
# def dative_bond(metals, non_metals):
#     query = ccdc.search.QuerySubstructure()
#     # Add the metal atom
#     metal = query.add_atom(ccdc.search.QueryAtom(metals))
#     # Add the ligand atom
#     ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
#     # Add pi bond between the metal and the ligand
#     bond = ccdc.search.QueryBond('Pi')
#     query.add_bond(bond, metal, ligand)
#     return query

# # DATIVE BONDS in rings (Iterative search)
# def dative_in_ring(num_intermediates):
#     query = ccdc.search.QuerySubstructure()
#     non_pi_bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
#     # Add the metal atom
#     metal = query.add_atom(ccdc.search.QueryAtom(metals))
#     previous_atom = metal
#     # Add the ligand atom(s) #start with at least one ligand
#     for _ in range(num_intermediates):
#         ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
#         bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
#         query.add_bond(bond, previous_atom, ligand)
#         previous_atom = ligand 
#     # insert the element in the ring
#     non_metal = query.add_atom(ccdc.search.QueryAtom(non_metals)) 
#     query.add_bond(non_pi_bond, previous_atom, non_metal)
#     # Add dative bond between the metal and the non-metal
#     dative_bond = ccdc.search.QueryBond('Pi')
#     query.add_bond(dative_bond , metal, non_metal)
#     return query

# # SIGMA BONDS (AGOSTIC INTERACTIONS) not in rings
# def sigma_bond(metals, non_metals):
#     query = ccdc.search.QuerySubstructure()
#     # Add the metal atom
#     metal = query.add_atom(ccdc.search.QueryAtom(metals))
#     # Add the ligand atom
#     ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
#     # Add sigma bond between the metal and the ligand
#     bond = ccdc.search.QueryBond('Sigma')
#     query.add_bond(bond, metal, ligand)
#     return query

# # SIGMA BONDS (AGOSTIC INTERACTIONS) in rings
# def sigma_bonds_in_ring(num_intermediates):
#     query = ccdc.search.QuerySubstructure()
#     non_pi_bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
    
#     # Add the metal atom
#     metal = query.add_atom(ccdc.search.QueryAtom(metals))
#     previous_atom = metal
#     # Add the ligand atom(s) #start with at least one ligand
#     for _ in range(num_intermediates):
#         ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
#         bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
#         query.add_bond(bond, previous_atom, ligand)
#         previous_atom = ligand 
#     # insert the element in the ring
#     hydrogen = query.add_atom(ccdc.search.QueryAtom('H')) 
#     query.add_bond(non_pi_bond, previous_atom, hydrogen)
#     # Add dative bond between the metal and the non-metal
#     dative_bond = ccdc.search.QueryBond('Pi')
#     query.add_bond(dative_bond , metal, hydrogen)
#     return query
