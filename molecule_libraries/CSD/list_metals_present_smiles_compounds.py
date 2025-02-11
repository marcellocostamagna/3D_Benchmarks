import pandas as pd
import numpy as np

csv_file = 'molecule_libraries/CSD/final_structures_with_smiles.csv'

# get all the entries that have a SMILES string not '' and not NaN
df = pd.read_csv(csv_file)

df = df[df['SMILES'].notnull()]

# Extract the information under 'Metals' and remove the brackets
# to extract the elements symbols. Create a set to remove duplicates

metals = set()

for metal_txt in df['Metals']:
    if pd.isnull(metal_txt):
        continue
    # if there is a , there are multiple metals. Collect them all
    if ',' in metal_txt:
        # Get eveything between the commas as different elements in the list
        metals_to_add = metal_txt.split(',')
    else:
        # If there is only one metal, add it to the list
        metals_to_add = [metal_txt]
    
    # The metals have a () after them with the number of occurences. Remove them
    for metal in metals_to_add:
        m = metal.strip().split('(')[0]
        metals.add(m)
        
# Save the metals to a file
with open('molecule_libraries/CSD/metals_present.txt', 'w') as f:
    for metal in metals:
        f.write(f'{metal}\n')
    
  
    
