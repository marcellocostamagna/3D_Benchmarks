import ccdc 
import pandas
import numpy

# text smiles with the entry and the similarity score
sim_text = '/Users/marcellocostamagna/3D_Benchmarks/CSD_similarity/similarity_results/similarity_results_grubbs_4d.txt'

# extract the names before ':' 
identifiers = []
with open(sim_text, 'r') as f:
    for line in f:
        identifiers.append(line.split(':')[0].strip())  
        
# retrive 100 smiles skipping the first one
# identifiers = identifiers[1:101]
        
# csv file with the smiles
csv_file = '/Users/marcellocostamagna/3D_Benchmarks/molecule_libraries/CSD/final_structures_with_smiles.csv'

# extract the smiles from the csv file for the identifiers
df = pandas.read_csv(csv_file)
smiles = []
for identifier in identifiers:
    smiles.append(df[df['Identifier'] == identifier]['SMILES'].values[0])
    
# write the smiles to a file
with open('/Users/marcellocostamagna/3D_Benchmarks/molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_for_benchmark_test/smiles_for_grubbs_4d.txt', 'w') as f:
    for smile in smiles:
        f.write(f'{smile}\n')