from ccdc import io
from hsr import similarity as sim
from hsr import fingerprint as fp
from hsr import pre_processing as pp
from hsr import pca_transform as pca
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count, Manager
import time
import re
import os
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSR
from rdkit import Chem


def hsr_fingerprint(molecule):
    """Generates a fingerprint from a molecule object."""
    molecule_array = []
    for atom in molecule.atoms:
        if atom.coordinates is None:
            continue
        x, y, z = atom.coordinates
        # atom_array = np.array([x, y, z, np.sqrt(atom.atomic_number)])
        # atom_array = np.array([x, y, z])
        atom_array = [x, y, z]
        molecule_array.append(atom_array)
    # Centering data
    molecule_array = molecule_array - np.mean(molecule_array, axis=0)
    # print(np.array(molecule_array))
    return fp.generate_fingerprint_from_data(np.array(molecule_array))

# setmolecule1 = pp.load_molecules_from_sdf('target_molecule.sdf')
# fp_1 = fp.generate_fingerprint_from_molecule(molecule1)
molecule1 = io.MoleculeReader('target_molecule.sdf')[0]
fp_1 = hsr_fingerprint(molecule1)


csd_reader = io.EntryReader('CSD')
molecule2 = csd_reader.entry('POZHIW').molecule
fp_2 = hsr_fingerprint(molecule2)

similarity = sim.compute_similarity_score(fp_1, fp_2)



# print(f'Fingerprint 1: {fp_1}') 
# print(f'Fingerprint 2: {fp_2}')
print(f'Similarity: {similarity}')

# rdkit USR similarity
# convert molecules to rdkit molecules
# save the two molecules to sdf files
mol_writer = io.MoleculeWriter('molecule1.sdf')
mol_writer.write(molecule1)
mol_writer.close()
mol_writer = io.MoleculeWriter('molecule2.sdf')
mol_writer.write(molecule2)
mol_writer.close()

suppl = Chem.SDMolSupplier('molecule1.sdf', removeHs=False, sanitize=False)
molecule1 = next(suppl)
# print the coordinates of the molecule
# for i,atom in enumerate(molecule1.GetAtoms()):
#     coords = molecule1.GetConformer().GetPositions()
#     symbol = atom.GetSymbol()
#     print(f'{symbol}: {coords[i]}')
array = pp.molecule_to_ndarray(molecule1, features=None)
# print(f' ARRAY :{array}')
suppl = Chem.SDMolSupplier('molecule2.sdf', removeHs=False, sanitize=False)
molecule2 = next(suppl)
fp_1 = GetUSR(molecule1)
fp_2 = GetUSR(molecule2)
similarity = GetUSRScore(fp_1, fp_2)
# print(f'Fingerprint 1: {fp_1}')
# print(f'Fingerprint 2: {fp_2}')
print(f'Similarity: {similarity}')

# test from sdf files
mol_1 = io.MoleculeReader('molecule1.sdf')[0]
# print the coordinates of the molecule
# for atom in mol_1.atoms:
#     symbol = atom.atomic_symbol
#     x, y, z = atom.coordinates
#     print(f'{symbol}: {x}, {y}, {z}')
mol_2 = io.MoleculeReader('molecule2.sdf')[0]
fp_1 = hsr_fingerprint(mol_1)
fp_2 = hsr_fingerprint(mol_2)
similarity = sim.compute_similarity_score(fp_1, fp_2)
print(f'Similarity: {similarity}')

