import sys
from openbabel import pybel
import openbabel as ob
import pymol
from pymol import cmd
from pymol.cgo import *

def read_sdf(file_path):
    # Read the SDF file using Pybel
    mols = list(pybel.readfile("sdf", file_path))
    return mols

def calculate_center(points):
    x, y, z = zip(*points)
    center_x = sum(x) / len(x)
    center_y = sum(y) / len(y)
    center_z = sum(z) / len(z)
    return [center_x, center_y, center_z]

def modify_bonds(mol):
    # Remove bonds between atom 0 and specified atoms
    bonds_to_remove = [2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 32]  
    for obatom in mol:
        if obatom.coordidx == 0:
            neighbour_atoms = ob.OBAtomAtomIter(obatom)
            
                # bond = obatom.GetBond(neighbour_atom)
                # if neighbour_atom.GetIdx() in bonds_to_remove:
                #     mol.DeleteBond(bond)
    return mol

def save_sdf(mol, file_path):
    # Save the modified molecule as an SDF file using Pybel
    mol.write("sdf", file_path, overwrite=True)

def visualize_in_pymol(sdf_file, center1, center2):
    # Visualize the molecule in PyMOL and add dotted lines
    pymol.finish_launching()
    cmd.load(sdf_file, 'molecule')
    cmd.show('sticks')
    cmd.zoom()

    # Add dotted lines
    obj = [
        BEGIN, LINES,
        COLOR, 1.0, 0.0, 0.0,
        VERTEX, center1[0], center1[1], center1[2],
        VERTEX, center2[0], center2[1], center2[2],
        END
    ]
    cmd.load_cgo(obj, 'dotted_lines')
    
    # Add dotted lines between atom 0 and the two points defined
    atom_0_pos = cmd.get_atom_coords('molecule and id 1')
    obj1 = [
        BEGIN, LINES,
        COLOR, 0.0, 1.0, 0.0,
        VERTEX, atom_0_pos[0], atom_0_pos[1], atom_0_pos[2],
        VERTEX, center1[0], center1[1], center1[2],
        END
    ]
    obj2 = [
        BEGIN, LINES,
        COLOR, 0.0, 1.0, 0.0,
        VERTEX, atom_0_pos[0], atom_0_pos[1], atom_0_pos[2],
        VERTEX, center2[0], center2[1], center2[2],
        END
    ]
    cmd.load_cgo(obj1, 'line_to_center1')
    cmd.load_cgo(obj2, 'line_to_center2')

# Main code starts here

sdf_file = 'polymerization.sdf'
mols = read_sdf(sdf_file)

if not mols:
    print("No valid molecules found in the SDF file.")
    sys.exit(1)

mol = mols[0]

# Get the center of atoms 5, 6, 7, 8, 34 and the center of atoms 4, 9, 10, 14, 15
points1 = []
points2 = []
for atom in mol:
    idx = atom.idx - 1  # Convert to 0-based index
    pos = atom.coords
    if idx in [4, 5, 6, 7, 33]:
        points1.append(pos)
    elif idx in [3, 8, 9, 13, 14]:
        points2.append(pos)

center1 = calculate_center(points1)
center2 = calculate_center(points2)

# Modify bonds
mol = modify_bonds(mol)

# Save the modified molecule
modified_sdf_file = 'modified_polymerization.sdf'
save_sdf(mol, modified_sdf_file)

# Visualize in PyMOL
visualize_in_pymol(modified_sdf_file, center1, center2)
