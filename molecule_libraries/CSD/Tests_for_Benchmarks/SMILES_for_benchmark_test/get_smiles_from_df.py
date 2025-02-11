import pandas as pd

df = pd.read_csv('./molecule_libraries/CSD/final_structures_with_smiles.csv')

# Get all the smiles strings
smiles = df['SMILES'].tolist()

# Get only the smiles that contain Ru
smiles = [smile for smile in smiles if type(smile) == str]
# smiles = [smile for smile in smiles if 'Ru' in smile]

# Save the smiles to a text file
output_dir = './molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_for_benchmark_test'
with open(f'{output_dir}/smiles.txt', 'w') as file:
    for smile in smiles:
        # do not write Nan values
        if type(smile) == str:
            file.write(f"{smile}\n")
