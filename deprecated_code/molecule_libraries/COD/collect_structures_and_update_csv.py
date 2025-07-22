import os
import shutil
import pandas as pd
from multiprocessing import Pool, cpu_count
import time

def copy_file(file_info, destination_directory):
    root, file = file_info
    source_path = os.path.join(root, file)
    destination_path = os.path.join(destination_directory, file)
    shutil.copy2(source_path, destination_path)
    return destination_path  # Return the new path for the file

def collect_and_copy_cif_files(directory, cod_directory, csv_file, new_csv_file):
    all_cif_files = []
    start = time.time()
    print(f'\nListing all files in {directory}...')
    for root, dirs, files in os.walk(directory):
        root_files = [(root, file) for file in files if file.endswith('.cif')]
        all_cif_files.extend(root_files)
    print(f'Listed all files in {time.time() - start:.2f} seconds.')
    print(f'Found {len(all_cif_files)} CIF files in the directory.\n')
    
    # Ensure the destination directory exists
    if not os.path.exists(cod_directory):
        os.makedirs(cod_directory)

    # Copy files using multiprocessing for efficiency
    print(f'Starting copy of CIF files to {cod_directory}...')
    with Pool(processes=cpu_count()) as pool:
        new_paths = pool.starmap(copy_file, [(file_info, cod_directory) for file_info in all_cif_files])

    print(f'Completed copy of CIF files to {cod_directory}.')

    # Load the original DataFrame
    df = pd.read_csv(csv_file)

    # Map old paths to new paths
    path_mapping = {os.path.join(root, file): os.path.join(cod_directory, file) for root, file in all_cif_files}
    
    # Update the DataFrame
    df['File Name'] = df['File Name'].apply(lambda x: path_mapping[x])

    # Save the updated DataFrame
    df.to_csv(new_csv_file, index=False)

    return df

if __name__ == '__main__':
    general_start = time.time()
    
    # Original csv file
    csv_file = 'final_structures_refined.csv'
    
    # Define your original directory containing CIF files
    original_directory = 'cif'
    
    # Create new directory to store the files
    cod_directory = 'COD_structures'
    
    # Create a new csv file to store the updated file paths
    new_csv_file = 'cod_structures.csv'

    # Process directory to collect all files and copy them
    updated_df = collect_and_copy_cif_files(original_directory, cod_directory, csv_file, new_csv_file)
    
    print(f'Total time taken: {time.time() - general_start:.2f} seconds.')
