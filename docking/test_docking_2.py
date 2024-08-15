from vina import Vina
# from ccdc import io
import oddt.toolkits.extras.rdkit as ob
import rdkit.Chem as Chem
import pymol
from pymol import cmd


def prepare_molecule(input_file, output_file, ligand=False):
    # Load molecule using Pybel
    # mol = toolkit.readfile('mol2', input_file)

    # open molecule with rdKit
    # mol_rdkit = Chem.MolFromMol2File(input_file, sanitize=False, removeHs=False)
    
    # read mol from sdf file
    suppl = Chem.SDMolSupplier(input_file, sanitize=False, removeHs=False )
    mol_rdkit = next(suppl)
    
    mol_pdbqt = ob.MolToPDBQTBlock(mol_rdkit, flexible=False)
    
    if ligand:
        # get the first line of the pdbqt file
        first_line = mol_pdbqt.split('\n')[0]
        # add 'ROOT' in the line after REMARK by the remark line with remark line + \n + ROOT
        mol_pdbqt = mol_pdbqt.replace(first_line, f'{first_line}\nROOT')
        # Add at the end of the file the line 'ENDROOT \n TORSDOF 0
        mol_pdbqt += '\nENDROOT\nTORSDOF 0\n'
    
    with open(output_file, 'w') as f:
        f.write(mol_pdbqt)
        
    return mol_pdbqt
    

def calculate_box_parameters(mol_path):
    
    # Open the molecule file with rdkit
    suppl = Chem.SDMolSupplier(mol_path, sanitize=False, removeHs=False )
    mol = next(suppl)
    
    # Calculate the geometric center of the molecule
    # coordinates_1 = mol.GetConformer().GetPositions()
    coordinates = []
    for atom in mol.GetAtoms():
        position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        coordinates.append(position)
        
    # get the geometric center of the molecule
    x, y, z = zip(*coordinates)
    center_x = sum(x) / len(x)
    center_y = sum(y) / len(y)
    center_z = sum(z) / len(z)
    center = [center_x, center_y, center_z]

    # Define box size as the range in each dimension plus a buffer
    buffer = 50  # buffer in angstroms to extend the box beyond the furthest atoms
    size_x = max(x) - min(x) + buffer
    size_y = max(y) - min(y) + buffer
    size_z = max(z) - min(z) + buffer
    box_size = [size_x, size_y, size_z]

    return center, box_size

def dock_molecules(receptor_path, ligand_path, output_path, center, box_size):
    
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_path)
    v.set_ligand_from_file(ligand_path)


    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=32, n_poses=20)
    v.write_poses(output_path, n_poses=5, overwrite=True)

    print(f"Docking completed. Results are saved in {output_path}")
    
# Covert file to pdbqt format
receptor = prepare_molecule('receptor.sdf', 'receptor.pdbqt', ligand=False)
ligand = prepare_molecule('ligand.sdf', 'ligand.pdbqt', ligand=True)

# Calculate center and box size based on receptor
center, box_size = calculate_box_parameters('receptor.sdf')

receptor_input_path = 'receptor.pdbqt'
ligand_input_path = 'ligand.pdbqt'
final_output_path = 'docked_output.pdbqt'

dock_molecules(receptor_input_path, ligand_input_path, final_output_path, center, box_size)

# Visualize the results with PyMOL
# visualize_with_pymol('receptor.pdbqt', final_output_path)

pymol.finish_launching()
cmd.load('receptor.pdbqt', 'receptor')
cmd.load('docked_output.pdbqt', 'docked')

# cmd.show_as('cartoon', 'receptor')
# cmd.show_as('sticks', 'docked')
cmd.zoom('receptor')