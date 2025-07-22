import pandas as pd
from ccdc import io
from ccdc import molecule
import os
import subprocess

# Head of the DataFrame
# Identifier,Metals,Component Count,Components With Metals,Polymolecular,Polymetallic,Organometallic,Explicit Hs,Hapticity,Delocalized,N_Atoms,Has Bonds,Metal_Bonds,Overlapping Atoms,All Bonds
# AACANI10,Ni(1),3,01: Ni(1),True,False,True,False,False,False,42.0,True,True,False,True
# AACANI11,Ni(1),3,01: Ni(1),True,False,True,True,False,False,42.0,True,True,False,True
# AACRUB,Ru(2),2,01: Ru(2),True,True,True,False,False,True,41.0,True,True,False,True
# AACRUB01,Ru(2),2,01: Ru(2),True,True,True,False,False,True,41.0,True,True,False,True

if __name__ == '__main__':
    # Load the CSV file into a DataFrame
    df = pd.read_csv('final_structures.csv')
    
    # Get, list and visulize DataFrame information
    
    # Numer of entries
    print(f'The number of entries is {len(df)}\n')
    
    # Number of columns (features)
    print(f'The number of features is {len(df.columns)}:\n')
    
    # List of features with description
    print(f'First feature: {df.columns[0]} - Identifier: CSD Identifier of the structure.')
    print(f'Second feature: {df.columns[1]} - Metals: Defines the metals in the structure and the number of occurrences of each metal. e.g. "Fe(3), Mn(1), Ti(1)"')
    print(f'Third feature: {df.columns[2]} - Component Count: Number of components in the structure, i.e. number of separarate molecules in the structure.')
    print(f'Fourth feature: {df.columns[3]} - Components With Metals: List of components with metals. e.g. "04: Fe(1); 05: Mn(1)"')
    print(f'Fifth feature: {df.columns[4]} - Polymolecular: True if the structure is polymolecular, False otherwise.')
    print(f'Sixth feature: {df.columns[5]} - Polymetallic: True if the structure has multiple occurrences of metals, False otherwise.')
    print(f'Seventh feature: {df.columns[6]} - Organometallic: True if the structure posses a metal and a Carbon atom, False otherwise.')
    print(f'Eighth feature: {df.columns[7]} - Explicit Hs: True if the structure has explicit Hydrogen atoms, False otherwise.')
    print(f'Tenth feature: {df.columns[9]} - Hapticity: True if the structure has hapto bonds, False otherwise.')
    print(f'Eleventh feature: {df.columns[10]} - N_Atoms: Number of atoms in the structure.')
    print(f'Twelfth feature: {df.columns[11]} - Connectivity: Complete, Partial or None.')
    print(f'Twelfth feature: {df.columns[11]} - Metal_Bonds: True if the structure has bonds involving metal atoms, False otherwise.')
    print(f'Thirteenth feature: {df.columns[12]} - Overlapping Atoms: True if the structure has overlapping atoms, False otherwise.')
    
    # Get the number of entries for different properties
    
    # Get the number of entries with hapticity
    n_hapticity = len(df[df['Hapticity'] == True])
    print(f'The number of entries with hapticity is {n_hapticity}')
    
    # Get the number of entries with polymolecular structures
    n_polymolecular = len(df[df['Polymolecular'] == True])
    print(f'The number of entries with polymolecular structures is {n_polymolecular}')
    
    # Get the number of entries with polymetallic structures
    n_polymetallic = len(df[df['Polymetallic'] == True])
    print(f'The number of entries with polymetallic structures is {n_polymetallic}')
    
    # Get the number of entries with organometallic structures
    n_organometallic = len(df[df['Organometallic'] == True])
    print(f'The number of entries with organometallic structures is {n_organometallic}')
    
    # Get the number of entries with metal bonds
    n_metal_bonds = len(df[df['Metal_Bonds'] == True])
    print(f'The number of entries with metal bonds is {n_metal_bonds}')
    
    # Get the number of entries with overlapping atoms
    n_overlapping_atoms = len(df[df['Overlapping Atoms'] == True])
    print(f'The number of entries with overlapping atoms is {n_overlapping_atoms}')
    
    # Attempt of filtering the structures that should be interesting and polish
    # Filter the dtaframe to get the structures with these properties:
    # - No polymetallic structures
    # - Complete Connectivity
    # - Has Metal bonds
    # - No Overlapping Atoms
    
    filtered_df = df[(df['Polymetallic'] == False) & (df['Connectivity'] == 'Complete') & (df['Metal_Bonds'] == True) & (df['Overlapping Atoms'] == False)]
    
    print(f'The number of entries that could meet the criteria is {len(filtered_df)}')
    
    # # Save the filtered DataFrame to a CSV file
    # filtered_df.to_csv('good_structures_monometallic.csv', index=False)

    # # File paths
    # file_names = filtered_df['File Name'][:1000].tolist() 
    # file_names = [os.path.abspath(file_name) for file_name in file_names]

    # # Mercury command
    # mercury = '/Users/marcellocostamagna/CCDC/ccdc-software/mercury/mercury.app/Contents/MacOS/mercury'

    # # Open the files with mercury
    # subprocess.run([mercury] + file_names, check=True)
        
        
    
    
    
    
    
    # # Get entries that have bonds
    # filtered_df = df[df['Overlapping Atoms'] == True]
    # print(len(filtered_df))
    # # print the first 50 entries
    # print(filtered_df.head(50))
    
    # # get the structure with the minimum number of atoms
    # min_atoms = filtered_df['N_Atoms'].min()
    
    # # get the structures with the minimum number of atoms
    # min_atoms_df = filtered_df[filtered_df['N_Atoms'] > min_atoms]
    
    # # get the maximum number of atoms
    # max_atoms = filtered_df['N_Atoms'].max()
    # max_atoms_entry = filtered_df[filtered_df['N_Atoms'] == max_atoms]
    # print(f'The entry {max_atoms_entry["File Name"].values[0]} has the maximum number of atoms: {max_atoms}')
    
    # # count the number of structures
    # print(f'The number of structures with {min_atoms} atoms is {len(min_atoms_df)}')
    
    # # Get the entry file name
    # min_atoms_entry = filtered_df[filtered_df['N_Atoms'] == min_atoms]
    # print(f'The entry {min_atoms_entry["File Name"].values[0]} has the minimum number of atoms: {min_atoms}')
    
    
    
    # # Iterate over the DataFrame
    # for index, row in df.iterrows():
    #     # Extract the file name
    #     file_path = row['File Name']
    #     # read the molecule
    #     mol_reader = io.MoleculeReader(file_path)
    #     mol = mol_reader[0]
    #     mol.assign_bond_types()
    #     print(f'Processing {file_path}...')
    #     # Iterate through the bonds
    #     for bond in mol.bonds:
    #         # print bond type
    #         print(bond.bond_type)
    #         if bond.bond_type == 9:
    #             # STOP the program
    #             print(f'Found bond type 9 in {file_path}.')
    #             exit()
    #     print(f'Next file...')
            

    