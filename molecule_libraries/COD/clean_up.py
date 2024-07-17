import os
import pandas as pd
from multiprocessing import Pool, cpu_count
import time

# Load the CSV
csv_file = 'final_structures_refined.csv'  
df = pd.read_csv(csv_file)

# Extract the full paths of valid files from the CSV
valid_files = set(df['File Name'].str.strip())

def remove_file(file_tuple):
    root, file = file_tuple
    file_path = os.path.join(root, file)
    if file_path not in valid_files:
        os.remove(file_path)
        return file_path
    return None

def process_directory(directory):
    all_cif_files = []
    start = time.time()
    print(f'\nListing all files in {directory}...')
    for root, dirs, files in os.walk(directory):
        # Create tuples of root and file for multiprocessing
        root_files = [(root, file) for file in files if file.endswith('.cif')]
        all_cif_files.extend(root_files)
    print(f'Listed all files in {time.time() - start:.2f} seconds.')
    print(f'Found {len(all_cif_files)} CIF files in the directory.\n')
    return all_cif_files

if __name__ == '__main__':
    
    general_start = time.time()
    # Define your original directory containing CIF files
    original_directory = 'cif'  

    # Process directory to collect all files
    all_files = process_directory(original_directory)

    print('Removing unwanted files...')
    # Use multiprocessing to remove files
    with Pool(cpu_count()) as pool:
        removed_files = pool.map(remove_file, all_files)
        
    # Filter None values and count removed files
    removed_files = [file for file in removed_files if file is not None]
 
    print(f"Removed {len(removed_files)} unwanted files from the directory.")
    print(f'Total time taken: {time.time() - general_start:.2f} seconds.')
