import ccdc
from ccdc import io
from ccdc.io import MoleculeReader, MoleculeWriter

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

csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry('AVIHEU')
molecule = component_of_interest(entry.molecule)

# Save AVIHEU on sdf
with io.MoleculeWriter('AVIHEU.sdf') as mol_writer:
    mol_writer.write(molecule)