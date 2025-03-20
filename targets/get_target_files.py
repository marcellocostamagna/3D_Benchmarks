import ccdc.io
import pandas as pd
from multiprocessing import Pool, cpu_count
import hsr
import time
import numpy as np

# Function to determine the component of interest
def component_of_interest(molecule):
    components = molecule.components
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
    
    heaviest_component = max(properties, key=lambda x: x["molecular_weight"])
    most_atoms_component = max(properties, key=lambda x: x["atom_count"])
    
    for prop in properties:
        criteria_met = sum([
            prop["is_organometallic"],
            prop["component"] == heaviest_component["component"],
            prop["component"] == most_atoms_component["component"]
        ])
        if criteria_met >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None


if __name__ == "__main__":

    targets = ['ABAHIW', 'ABAKIZ', 'ABADOX', 
               'ABABIP', 'GASQOK', 'ABEKIE',
               'NIWPUE01','ABEKIF', 'APUFEX', 
               'ABEHAU', 'TITTUO', 'EGEYOG',
               'ABOBUP', 'XIDTOW', 'ACNCOB10',
               'TACXUQ', 'ACAZFE', 'NIVHEJ',
               'ADUPAS', 'DAJLAC', 'OFOWIS',
               'CATSUL', 'HESMUQ01', 'GUDQOL',
               'ABEVAG', 'AKOQOH', 'ADARUT',
               'AFECIA', 'ACOVUL', 'AFIXEV', ]
    
    entry_reader = ccdc.io.EntryReader('CSD')
    
    for index, target_entry in enumerate(targets, start=1):
        print(f"Processing target: {target_entry}")
        target_molecule = component_of_interest(entry_reader.entry(target_entry).molecule)
        # get the SMILES of the component of interest of the target molecule
        target_smiles = target_molecule.smiles
        print(f"Target {target_entry}: {target_smiles}")
        
        # Save the target entry and the SMILES separated by a : in a text file with the target name
        with open(f'./targets/{index}_{target_entry}.txt', 'w') as f:
            f.write(f"{target_entry}:{target_smiles}")
        print(f"Target {target_entry} processed.")
        
        
      