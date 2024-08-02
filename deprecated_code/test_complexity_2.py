from ccdc import io
import pandas as pd
from multiprocessing import Pool, cpu_count
import time
import json

def calculate_ratio(identifier):
    try:
        csd_reader = io.EntryReader('CSD')
        entry = csd_reader.entry(identifier)
        num_bonds = len(entry.molecule.bonds)
        num_atoms = len(entry.molecule.atoms)
        ratio = num_bonds / num_atoms
    except Exception as e:
        print(f'Error with identifier {identifier}: {e}')
        return identifier, {"identifier": identifier, "num_bonds": 0, "num_atoms": 0, "ratio": 0}
    
    return identifier, {"identifier": identifier, "num_bonds": num_bonds, "num_atoms": num_atoms, "ratio": ratio}

def get_max_ratios(ratios, n=10):
    sorted_ratios = sorted(ratios.items(), key=lambda x: x[1]["ratio"], reverse=True)
    return sorted_ratios[:n]

def get_closest_ratios(ratios, target_ratio, n=5):
    above_target = {identifier: data for identifier, data in ratios.items() if data["ratio"] > target_ratio}
    below_target = {identifier: data for identifier, data in ratios.items() if data["ratio"] < target_ratio}
    
    above_sorted = sorted(above_target.items(), key=lambda x: abs(x[1]["ratio"] - target_ratio))
    below_sorted = sorted(below_target.items(), key=lambda x: abs(x[1]["ratio"] - target_ratio))
    
    closest_above = above_sorted[:n]
    closest_below = below_sorted[:n]
    
    return closest_below, closest_above

def save_ratios_to_file(ratios, filename):
    with open(filename, 'w') as f:
        json.dump(ratios, f, indent=4)

if __name__ == '__main__':
    start = time.time()
    csd_reader = io.EntryReader('CSD')

    df = pd.read_csv('final_structures_with_fp.csv')
    entries = df['Identifier'].values

    with Pool(processes=cpu_count()) as pool:
        results = pool.map(calculate_ratio, entries)

    ratios = dict(results)
    
    max_ratios = get_max_ratios(ratios)
    min_ratio = min(ratios, key=lambda x: ratios[x]["ratio"])
    average_ratio = sum(data["ratio"] for data in ratios.values()) / len(ratios)
    
    # Save the ratios dictionary to a file
    save_ratios_to_file(ratios, 'ratios_1.json')
    
    # Get the 10 entries that have ratio closest around 1.5
    target_ratio = 1.5
    closest_below, closest_above = get_closest_ratios(ratios, target_ratio) 
    
    print(f'Max ratio: {max_ratios[0][0]} with ratio {max_ratios[0][1]["ratio"]}')
    for i, (identifier, data) in enumerate(max_ratios[1:], start=2):
        print(f'Max ratio {i}: {identifier} with ratio {data["ratio"]}')
        
    print(f'\nClosest ratios below {target_ratio}:')
    for i, (identifier, data) in enumerate(closest_below, start=1):
        print(f'{i}: {identifier} with ratio {data["ratio"]:.5f}')
        
    print(f'\nClosest ratios above {target_ratio}:')
    for i, (identifier, data) in enumerate(closest_above, start=1):
        print(f'{i}: {identifier} with ratio {data["ratio"]:.5f}')

    print(f'\nMin ratio: {min_ratio} with ratio {ratios[min_ratio]["ratio"]}')
    print(f'Average ratio: {average_ratio}\n')
    print(f'Time taken: {(time.time() - start)/60} minutes')
