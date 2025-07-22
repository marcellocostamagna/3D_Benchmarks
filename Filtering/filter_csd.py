# Script to filter and analyse the CSD to retrieve structures of interest.
#
# Filtering: step to eliminate structures that do not meet these minimum criteria:
#
#   1. Has 3D coordinates
#   2. Is not polymeric
#   3. Has explicit hydrogens
#   4. Component of interest has more than 5 atoms
#   5. Component of interest has no disordered/not overlapping atoms
#   6. Component of interest has SMILES
#
# Analysis: collects the structures that meet the above criteria and saves them in a
# CSV file with the following information:
#
#   1. Identifier
#   2. Explicit Hydrogens (True/False)
#   3. Number of atoms in the component of interest
#   4. Connectivity (Complete/Partial/None)
#   5. Overlapping Atoms (check of disordered atoms in the component)
#   6. Disorder (True/False) (crystal)
#   7. SMILES of the component of interest

import csv
from ccdc import io
import time
from multiprocessing import Pool, cpu_count
from ccdc import descriptors

def check_connectivity(molecule):
    "Check the connectivity of the molecule"
    # Check if the molecule has bonds
    if len(molecule.bonds) > 0:
        # Check if the connecticity is partial or complete
        connectivity = 'Partial' if any(len(atom.bonds) == 0 for atom in molecule.atoms) else 'Complete'
    else: 
        connectivity = 'None'
    
    return connectivity

def component_of_interest(molecule):
    """
    Returns the component of interest: The component that checks at least two of the following criteria:
    1. Is organometallic
    2. Is the heaviest component
    3. Has the most atoms
    
    If the chosen component has fewer than 5 atoms, returns None.
    """
    components = molecule.components
    
    # Evaluate properties for each component
    properties = []
    for comp in components:
        is_organometallic = comp.is_organometallic
        molecular_weight = sum(atom.atomic_weight for atom in comp.atoms)
        atom_count = len(comp.atoms)
        
        properties.append({
            "component": comp,
            "is_organometallic": is_organometallic,
            "molecular_weight": molecular_weight,
            "atom_count": atom_count
        })
    # Determine heaviest component and component with most atoms
    heaviest_component = max(properties, key=lambda x: x["molecular_weight"])
    most_atoms_component = max(properties, key=lambda x: x["atom_count"])
    
    # Select the component that matches at least two criteria
    for prop in properties:
        criteria_met = 0
        if prop["is_organometallic"]:
            criteria_met += 1
        if prop["component"] == heaviest_component["component"]:
            criteria_met += 1
        if prop["component"] == most_atoms_component["component"]:
            criteria_met += 1
        
        if criteria_met >= 2:
            # Check if the selected component has at least 5 atoms
            if prop["atom_count"] >= 5:
                return prop["component"]
    return None
        
def filter_and_analyse(entry_id):
    reader = io.EntryReader('CSD')
    entry = reader.entry(entry_id)
    # Filtering (operates on each entry in the CSD)
    molecule = entry.molecule
    crystal = entry.crystal
    # Minimum criteria to accept an entry
    
    # 1. Has 3D coordinates
    if not molecule.is_3d:
        return None
    
    # 2. Is not polymeric
    if molecule.is_polymeric:
        return None
    
    # 3. Has explicit hydrogens
    explicit_hydrogens = all(atom.atomic_number == 1 and atom.coordinates is not None for atom in molecule.atoms if atom.atomic_number == 1)
    if not explicit_hydrogens:
        return None
    
    # 4. Component of interest has more than 5 atoms
    comp_of_interest = component_of_interest(molecule)
    if not comp_of_interest:
        return None
    
    # 5. Component of interest has no overlapping atoms
    desc = descriptors.MolecularDescriptors()
    overlapping_atoms = any((d := desc.atom_distance(atom1, atom2)) is not None and d < 0.4
                            for atom1 in comp_of_interest.atoms
                            for atom2 in comp_of_interest.atoms if atom1 != atom2)
    if overlapping_atoms:
        return None
    
    # 6. Component of interest with SMILES
    if not comp_of_interest.smiles:
        return None
    
    # Analysis: collect and return required information only
    
    # 1. Identifier
    id = entry.identifier
    n_atoms = len(comp_of_interest.atoms)
    connectivity = check_connectivity(comp_of_interest)
    disorder = bool(crystal.has_disorder)
    smiles = comp_of_interest.smiles
    
    return id, explicit_hydrogens, n_atoms, connectivity, overlapping_atoms, disorder, smiles
          
def filtering_and_analysis(csd_reader):   
    entries = []
    stop = [False,100000] # Set to True and a limit to the number of entries to process for testing
    print("Retrieving entries from the CSD...")
    for entry in csd_reader.entries():
        entries.append(entry.identifier)
        
        if stop[0] and len(entries) == stop[1]:
            break
    print(f"Retrieved {len(entries)} entries from the CSD. In {(time.time() - start)/60} minutes.")
    print(f"Filtering and analysing the CSD ({len(entries)} entries) with {cpu_count()} CPUs...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(filter_and_analyse, entries)
    results = [result for result in results if result is not None]
    return results
        
def saving(results):
    # Saves the results to a CSV file
    with open('viable_structures.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Identifier', 
                         'Explicit Hs', 
                         'N_Atoms', 
                         'Connectivity', 
                         'Overlapping Atoms',
                         'Disorder',
                         'SMILES'])
        for result in results:
            writer.writerow(result)
    
def main():
    csd_reader = io.EntryReader('CSD')
    results = filtering_and_analysis(csd_reader)
    print(f"Filtering and analysis complete in {(time.time() - start)/60} minutes.")
    
    print(f" Saving {len(results)} entries to 'viable_structures.csv'...")
    saving(results)
    print(f"Saving complete in {(time.time() - start)/60} minutes.")
    print(f'Total time elapsed: {(time.time() - start)/60} minutes')

if __name__ == '__main__':
    print('Starting...')
    start= time.time()
    main()
    print(f'Finished in {(time.time() - start)/60} minutes')
