import ccdc.io
import numpy as np
import hsr
import os

# --- Utility Functions ---

def get_p_array_from_ccdcmol(ccdcmol):
    """Extract atomic coordinates + sqrt(Z) for HSR."""
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([
                atom.coordinates[0],
                atom.coordinates[1],
                atom.coordinates[2],
                np.sqrt(atom.atomic_number)
            ])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

def get_p_q_array_from_ccdcmol(ccdcmol):
    """Extract atomic coordinates + sqrt(Z) + formal charge."""
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([
                atom.coordinates[0],
                atom.coordinates[1],
                atom.coordinates[2],
                np.sqrt(atom.atomic_number),
                atom.formal_charge
            ])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

def component_of_interest(molecule):
    """Selects the key component from the molecule using 3 criteria."""
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

    targets = [
        'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
        'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
        'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
        'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
        'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV',
        'ABAYAF', 'RULJAM'
    ]

    output_dir = './target_files/'
    os.makedirs(output_dir, exist_ok=True)

    entry_reader = ccdc.io.EntryReader('CSD')

    for index, target_entry in enumerate(targets, start=1):
        print(f"Processing target: {target_entry}")
        target_molecule = component_of_interest(entry_reader.entry(target_entry).molecule)
        if target_molecule is None:
            print(f"  Skipped: no valid component of interest found.")
            continue
        # Get the SMILES of the component of interest
        target_smiles = target_molecule.smiles
        print(f"Target {target_entry}: {target_smiles}")

        # Get the fingerprint of the molecule
        if target_entry in ('ABAYAF', 'RULJAM'):
            fingerprint = hsr.generate_fingerprint_from_data(get_p_q_array_from_ccdcmol(target_molecule))
        else:
            fingerprint = hsr.generate_fingerprint_from_data(get_p_array_from_ccdcmol(target_molecule))

        # Save the target entry, the SMILES and the fingerprint
        filename = os.path.join(output_dir, f"{index}_{target_entry}.txt")
        with open(filename, 'w') as f:
            f.write(f"{target_entry} {target_smiles} {fingerprint}")
        print(f"Target {target_entry} processed and saved to {filename}.")

