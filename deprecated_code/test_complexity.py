from ccdc import io
import pandas as pd


csd_reader = io.EntryReader('CSD')

df = pd.read_csv('final_structures_with_fp.csv')

entries = df['Identifier'].values

# with open('batch_18.txt', 'r') as f:
#     # read the lines of the file
#     lines = f.readlines()
#     # remove the newline character from each line
#     entries = [line.strip() for line in lines]

entries = [csd_reader.entry(identifier) for identifier in entries]

ratios = {}
for entry in entries:
    # print(f'Identifier: {entry.identifier}')
    # print(f'N_ATOMS: {len(entry.molecule.atoms)}')
    # print(f'N_BONDS: {len(entry.molecule.bonds)}')
    # print ratio of bonds to atoms
    ratios[entry.identifier] = len(entry.molecule.bonds)/len(entry.molecule.atoms)
    # print(f'Ratio: {len(entry.molecule.bonds)/len(entry.molecule.atoms)}')
    
# print mac, min and average ratio with the corresponding identifier
copy_ratios = ratios.copy()
max_ratio = max(ratios, key=ratios.get)

min_ratio = min(ratios, key=ratios.get)
average_ratio = sum(ratios.values())/len(ratios)
print(f'Max ratio: {max_ratio} with ratio {ratios[max_ratio]}')

# get the other 10 max ratios
for i in range(10):
    copy_ratios.pop(max_ratio)
    max_ratio = max(copy_ratios, key=copy_ratios.get)
    print(f'Max ratio {i+2}: {max_ratio} with ratio {copy_ratios[max_ratio]}')



print(f'Min ratio: {min_ratio} with ratio {ratios[min_ratio]}')
print(f'Average ratio: {average_ratio}')
