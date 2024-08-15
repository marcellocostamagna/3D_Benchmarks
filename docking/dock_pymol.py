from vina import Vina
# from ccdc import io
import oddt.toolkits.extras.rdkit as ob
import rdkit.Chem as Chem
import pymol
from pymol import cmd
import threading

def prepare_molecule(input_file, output_file, ligand=False):
    suppl = Chem.SDMolSupplier(input_file, sanitize=False, removeHs=False)
    mol_rdkit = next(suppl)
    mol_pdbqt = ob.MolToPDBQTBlock(mol_rdkit, flexible=False)

    if ligand:
        first_line = mol_pdbqt.split('\n')[0]
        mol_pdbqt = mol_pdbqt.replace(first_line, f'{first_line}\nROOT')
        mol_pdbqt += '\nENDROOT\nTORSDOF 0\n'

    with open(output_file, 'w') as f:
        f.write(mol_pdbqt)

    return mol_pdbqt

def calculate_box_parameters(mol_path):
    suppl = Chem.SDMolSupplier(mol_path, sanitize=False, removeHs=False)
    mol = next(suppl)

    coordinates = [mol.GetConformer().GetAtomPosition(atom.GetIdx()) for atom in mol.GetAtoms()]
    x, y, z = zip(*coordinates)
    center = [sum(x) / len(x), sum(y) / len(y), sum(z) / len(z)]
    buffer = 10
    box_size = [max(x) - min(x) + buffer, max(y) - min(y) + buffer, max(z) - min(z) + buffer]

    return center, box_size

def dock_molecules(receptor_path, ligand_path, output_path, center, box_size):
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_path)
    v.set_ligand_from_file(ligand_path)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=32, n_poses=20)
    v.write_poses(output_path, n_poses=5, overwrite=True)
    print(f"Docking completed. Results are saved in {output_path}")

def visualize_with_pymol(receptor_file, docked_file):
    pymol.finish_launching()
    cmd.load(receptor_file, 'receptor')
    cmd.load(docked_file, 'docked')
    cmd.align('docked', 'receptor')
    cmd.show_as('cartoon', 'receptor')
    cmd.show_as('sticks', 'docked')
    cmd.zoom('receptor')
    cmd.save('docking_visualization.pse')

def main():
    # file = 'docking/ja052912c_si_002.cif'
    # mol_reader = io.MoleculeReader(file)
    # molecule = mol_reader[0]

    # for component in molecule.components:
    #     print(f'{component.identifier}: {len(component.atoms)} atoms and MW: {component.molecular_weight}')

    # comp1, comp2, comp3 = sorted(molecule.components, key=lambda x: x.molecular_weight, reverse=True)[:3]
    # print(f"Component 1: {comp1.identifier} with {len(comp1.atoms)} atoms")
    # print(f"Component 2: {comp2.identifier} with {len(comp2.atoms)} atoms")
    # print(f"Component 3: {comp3.identifier} with {len(comp3.atoms)} atoms")

    # receptor = comp1
    # ligand = comp3

    # with io.MoleculeWriter('receptor.sdf') as receptor_writer:
    #     receptor_writer.write(receptor)
    # with io.MoleculeWriter('ligand.sdf') as ligand_writer:
    #     ligand_writer.write(ligand)

    receptor = prepare_molecule('receptor.sdf', 'receptor.pdbqt', ligand=False)
    ligand = prepare_molecule('ligand.sdf', 'ligand.pdbqt', ligand=True)

    center, box_size = calculate_box_parameters('receptor.sdf')

    receptor_input_path = 'receptor.pdbqt'
    ligand_input_path = 'ligand.pdbqt'
    final_output_path = 'docked_output.pdbqt'

    dock_molecules(receptor_input_path, ligand_input_path, final_output_path, center, box_size)

    # Run PyMOL visualization on the main thread
    pymol_thread = threading.Thread(target=visualize_with_pymol, args=('receptor.pdbqt', final_output_path))
    pymol_thread.start()
    pymol_thread.join()

if __name__ == '__main__':
    main()
