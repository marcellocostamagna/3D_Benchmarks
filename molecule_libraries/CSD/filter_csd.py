# Script to filter and analyse the CSD to retrieve structures of interest.
#
# Filtering: step to eliminate structures that do not meet these minimum criteria:
#
# 1. Has 3D coordinates
# 2. Contain metals
# 3. Are not polymeric
# 4. Do not contain disordered atoms
#
# Analysis: collects the structures that meet the above criteria and saves them in a
# CSV file with the following information:
#
# 1. Identifier
# 2. Metals present in the structure
# 3. Number of components in the structure
# 4. Components with metals
# 5. Polymolecular (True/False) (these is deduced from the number of components)
# 6. Polymetallic (True/False)
# 7. Organometallic (True/False) (csd definition)
# 8. Explicicit Hydrogens (True/False) 
# 9. Hapticity
# 10. Number of atoms
# 11. Connectivity (Complete/Partial/None)
# 12. Metal Bonds
# 13. Overlapping Atoms (Second check of disordered atoms)

import csv
from ccdc import io
import time
import re
from collections import Counter
import pandas as pd
from multiprocessing import Pool, cpu_count
from ccdc import descriptors

# Define your transition metals and rare earth elements
transition_metals = [
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
]
rare_earths = [
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np', 
    'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
]
elements = set(transition_metals + rare_earths)

def parse_formula(formula):
    """Parses a chemical formula to count occurrences of each element."""
    element_pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    elements = element_pattern.findall(formula)
    element_count = Counter()
    for element, count in elements:
        count = int(count) if count else 1
        element_count[element] += count
    return element_count

def check_multiple_metals(element_count, metal_set):
    """Checks for multiple metals or multiple occurrences of the same metal."""
    metals_in_formula = {element: count for element, count in element_count.items() if element in metal_set}
    multiple_metals = sum(count for count in metals_in_formula.values()) > 1
    return multiple_metals, metals_in_formula

def metals_and_components(molecule):
    "Analyse the components of the structure and checks for metals"
    # analysis of each component
    components_with_metals = {}
    polymetallic_flag = False
    for component in molecule.components:
        component_formula = parse_formula(component.formula)
        multiple_metals, component_metals = check_multiple_metals(component_formula, elements)
        if component_metals:
            metals_list = ', '.join([f"{el}({count})" for el, count in component_metals.items()])
            components_with_metals[component.identifier] = metals_list
            if multiple_metals:
                polymetallic_flag = True
    
    # analysis of the whole structure    
    element_count = parse_formula(molecule.formula)
    _, metals_in_formula = check_multiple_metals(element_count, elements)
    metals_in_structure = ', '.join([f"{el}({count})" for el, count in metals_in_formula.items()])
    
    # summary of metals and components
    metals = metals_in_structure # 2
    components = len(molecule.components) # 3
    components_with_metals = '; '.join([f"{k}: {v}" for k, v in components_with_metals.items()]) # 4
    
    return metals, components, components_with_metals, polymetallic_flag

def check_connectivity(molecule):
    "Check the connectivity of the molecule"
    metal_bonds = False 
    # Check if the molecule has bonds
    if len(molecule.bonds) > 0:
        # Check if the connecticity is partial or complete
        connectivity = 'Partial' if any(len(atom.bonds) == 0 for atom in molecule.atoms) else 'Complete'
                
        # Check if there are bonds involving metals
        for bond in molecule.bonds:
            atom1, atom2 = bond.atoms
            if atom1.atomic_symbol in elements or atom2.atomic_symbol in elements:
                metal_bonds = True
                break   
    else: 
        connectivity = 'None'
        
    return connectivity, metal_bonds
        

def filter_and_analyse(entry_id):
    reader = io.EntryReader('CSD')
    entry = reader.entry(entry_id)
    # Filtering (operates on eache entry in the CSD)
    molecule = entry.molecule
    crystal = entry.crystal
    # Minimum criteria to accept an entry
    
    # 1. Has 3D coordinates
    if not molecule.is_3d:
        return None
    
    # 2. Contain metals
    if not any(element in entry.formula for element in elements):
        return None
    
    # 3. Are not polymeric
    if molecule.is_polymeric:
        return None
    
    # 4. Do not contain disordered atoms
    if crystal.has_disorder:
        return None
    
    # Analysis (operates on each molecule that passes the filtering)
    
    # 1. Identifier
    id = entry.identifier
    
    # 2/3/4. Components and metals
    metals, components, components_with_metals, polymetallic_flag = metals_and_components(molecule)
    
    # 5. Polymolecular (True/False)
    polymolecular = 'True' if components > 1 else 'False'
    
    # 6. Polymetallic (True/False)
    polymetallic = 'True' if polymetallic_flag else 'False'
    
    # 7. Organometallic (True/False)
    organometallic = 'True' if molecule.is_organometallic else 'False'
    
    # 8. Explicicit Hydrogens (True/False)
    explicit_hydrogens = all(atom.atomic_number == 1 and atom.coordinates is not None for atom in molecule.atoms if atom.atomic_number == 1)
    
    # 9. Hapticity (True if the molucule has at least a bond of type 9)
    hapticity = any(bond.bond_type == 9 for bond in molecule.bonds)
    
    # 10. Number of atoms
    n_atoms = len(molecule.atoms)
    
    # 11/12. Connectivity (Complete/Partial/None) and Metal Bonds (True/False)
    connectivity, metal_bonds = check_connectivity(molecule)
    
    # 13. Overlapping Atoms (Second check of disordered atoms)
    desc = descriptors.MolecularDescriptors()
    overlapping_atoms = any((d := desc.atom_distance(atom1, atom2)) is not None and d < 0.4 
                            for atom1 in molecule.atoms 
                            for atom2 in molecule.atoms if atom1 != atom2)   
    
    return id, metals, components, components_with_metals, polymolecular, polymetallic, organometallic, explicit_hydrogens, hapticity, n_atoms, connectivity, metal_bonds, overlapping_atoms
          
def filtering_and_analysis(csd_reader):   
    entries = []
    stop = [False,10000] # Set a limit to the number of entries to process for testing
    print("Retrieving entries from the CSD...")
    for entry in csd_reader.entries():
        entries.append(entry.identifier)
        
        if stop[0] and len(entries) == stop[1]:
            break
    print(f"Retrieved {len(entries)} entries from the CSD. In {(time.time() - start)/60} minutes.")
    print(f"Filtering and analysing the CSD ({len(entries)} entries) with {cpu_count()} CPUs...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(filter_and_analyse, entries)
    
    # Pass on only the results that are not None
    results = [result for result in results if result is not None]
    return results
        

def saving(results):
    # Saves the results to a CSV file
    with open('final_structures.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Identifier', 
                         'Metals', 
                         'Component Count', 
                         'Components With Metals', 
                         'Polymolecular', 
                         'Polymetallic', 
                         'Organometallic', 
                         'Explicit Hs', 
                         'Hapticity', 
                         'N_Atoms', 
                         'Connectivity', 
                         'Metal_Bonds', 
                         'Overlapping Atoms'])
        for result in results:
            writer.writerow(result)
    
def main():
    csd_reader = io.EntryReader('CSD')
    results = filtering_and_analysis(csd_reader)
    print(f"Filtering and analysis complete in {(time.time() - start)/60} minutes.")
    
    print(f" Saving {len(results)} entries to 'final_structures.csv'...")
    saving(results)
    print(f"Saving complete in {(time.time() - start)/60} minutes.")
    
    print(f'Total time elapsed: {(time.time() - start)/60} minutes')

if __name__ == '__main__':
    print('Starting...')
    start= time.time()
    main()
    print(f'Finished in {(time.time() - start)/60} minutes')
