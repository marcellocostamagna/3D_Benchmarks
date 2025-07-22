import pandas as pd

# Get the minimum values of N_atoms from the csv file
df = pd.read_csv('final_structures_last.csv')

# min_n_atoms = df['N_Atoms'].min()
# print(f'The smallest molecule has {min_n_atoms} atoms.')
# #     entries = df['Identifier'].tolist()

# # reorder the entries by the number of atoms and
# # Print the 10 smallest molecules
# df = df.sort_values(by='N_Atoms')
# print(df.head(10))

# chech if the 'Overlapping Atoms' column has any True values
overlapping_atoms = df['Overlapping Atoms'].unique()
print(overlapping_atoms)
# get the number of entries with overlapping atoms
overlapping_atoms = df[df['Overlapping Atoms'] == True]
print(f'The number of entries with overlapping atoms is {len(overlapping_atoms)}')
# print the first 10 entries with overlapping atoms
print(overlapping_atoms.head(10))

# check if the 'Connectivity' column has any partial values
connectivity = df['Connectivity'].unique()
print(connectivity)
# get the number of entries with partial connectivity
partial_connectivity = df[df['Connectivity'] == 'Partial']
print(f'The number of entries with partial connectivity is {len(partial_connectivity)}')
# print the first 10 entries with partial connectivity
print(partial_connectivity.head(10))

# Check if the 'Explicit Hs' column has any False values
explicit_hs = df['Explicit Hs'].unique()
print(explicit_hs)
# get the number of entries with no explicit hydrogen atoms
no_explicit_hs = df[df['Explicit Hs'] == False]
print(f'The number of entries with no explicit hydrogen atoms is {len(no_explicit_hs)}')
# print the first 10 entries with no explicit hydrogen atoms
print(no_explicit_hs.head(10))

