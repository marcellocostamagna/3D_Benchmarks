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
# subset_identifiers = get_random_identifiers(100, 1000)  
subset_identifiers = ['ABETUW', 'ACAZFE', 'ABIGIC', 'ABEDES']  

# DFEINITION OF THE QUERY SUBSTRUCTURE
# Chooose a query from the pattern_queries.py
query = dative_in_conjugated_chains(METALS, NON_METALS)

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
        
 