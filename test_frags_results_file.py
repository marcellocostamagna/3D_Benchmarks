from ccdc import io
import os
import numpy as np
from hsr import fingerprint as fp
from hsr import similarity as sim

def get_array_from_ccdcmol(ccdcmol):
    """Extracts an array of atomic coordinates and atomic numbers."""
    atom_array = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return atom_array - np.mean(atom_array, axis=0)  # Centering data

def get_array_from_ccdcmol_3d(ccdcmol):
    """Extracts an array of atomic coordinates and atomic numbers."""
    atom_array = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2]]
        for atom in ccdcmol.atoms
    ])
    return atom_array - np.mean(atom_array, axis=0)

src_dir = 'frag_comparison_results/'

#Count the number of files in the directory

file_count = sum([len(files) for r, d, files in os.walk(src_dir)])

for i in range(1, file_count):
    # open files of this format src_dir/*_frag{i}_matches.sdf
    file_name = src_dir + f'25_ABEVAG_frag{i}_matches.sdf'
    print(f"Processing {file_name}")
    mol_reader = io.MoleculeReader(file_name)
    mol_fps= []
    for mol in mol_reader:
        mol_fp = fp.generate_fingerprint_from_data(get_array_from_ccdcmol(mol))
        mol_fps.append(mol_fp)

    # Print similarity between the first molecule and the rest
    for i in range(1, len(mol_fps)):
        similarity = sim.compute_similarity_score(mol_fps[0], mol_fps[i])
        print(f"Similarity between molecule 0 and molecule {i}: {similarity:.4f}")

