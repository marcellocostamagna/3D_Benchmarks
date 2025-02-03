from pattern_queries import *
import ccdc.io
import ccdc.search

# Search for one query
def search(queries, identifiers):
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
        molecules = [ccdc.io.EntryReader('CSD').entry(id).molecule.heaviest_component for id in identifiers]
        results = search.search(molecules)
        results_ids = set(hit.identifier for hit in results)
    return results_ids

# Metals in rings
# one metal in a ring (chelating ligands)
# two metals ina aring (bridging ligands)
# more than two metals in a ring (clusters, not searched)
def metals_in_rings(METALS, NON_METALS, identifiers):
    '''
    Searches for the presence of metals in rings:
    - chelating ligands
    - Multi bridged complexes'''
    query1 = chelating_ligands(METALS, NON_METALS)
    query2 = multi_bridged_complexes(METALS, NON_METALS, 4) 
    results1 = search(query1, identifiers)
    results2 = search(query2, identifiers)
    # Subtract bridged complexes from chelating ligands
    results1 = results1 - results2 if results1 and results2 else results1
    # Return the results as a dictionary
    results_dict = {'Chelating ligands': results1, 'Multi bridged complexes': results2}
    return results_dict

def dative_patterns(METALS, NON_METALS, identifiers):
    query1 = dative_bonds(METALS, NON_METALS)
    query2 = cyclic_multihapto_ligands(METALS, NON_METALS)
    query3 = dative_in_rings(METALS, NON_METALS, 6)
    query4 = dative_in_conjugated_chains(METALS, NON_METALS)
    query5 = dative_in_conjugated_rings(METALS, NON_METALS)
    
    results1 = search(query1, identifiers)
    # Since the first search is the most general we can limit the next searches to the non hits of the first search
    results2 = search(query2, results1)
    # Since the second
    
    results3 = search(query3, results1)
    results4 = search(query4, results1)
    results5 = search(query5, results1)
            
    
    
    
    
    
    
    



