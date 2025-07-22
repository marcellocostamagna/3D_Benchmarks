from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeReader, MoleculeWriter


from ccdc import io
# print(io.MoleculeReader.csd_version())

# ### RETURNS THE BOND TYPES OF A MOLECULE ###
# csd_reader = io.EntryReader('CSD')
# entry = csd_reader.entry('DOYQOX')
# molecule = entry.molecule
# bond_types = {}
# for bond in molecule.bonds:
#     bond_type = bond.bond_type  # Get the numerical value of the bond type
#     text = bond.bond_type._bond_type.text()   # Get the text representation of the bond type
#     value = str(bond.bond_type._bond_type.value()) # Get the value of the bond type
    
#     key = f'{text}, ({value})'
    
#     if key in bond_types:
#         bond_types[key] += 1
#     else:
#         bond_types[key] = 1
# print(f'The molecule has {len(molecule.bonds)} bonds')
# print(f'The bond types are: {bond_types}')

# # ITERATE OVER THE MOLECULES AND COLLECT 20 MOLECULES
# # which are organic and have 3d coordinates
# organic_3d_mols = []
# for entry in csd_reader.entries():
#     molecule = entry.molecule
#     crystal = entry.crystal
#     if molecule.is_organic and molecule.is_3d and not crystal.has_disorder and len(molecule.components) == 1:
#         organic_3d_mols.append(molecule)
#     if len(organic_3d_mols) == 20:
#         break
    
# # print(f'Collected {len(organic_3d_mols)} organic molecules with 3D coordinates')
# # # print the identifiers of the collected molecules
# # for molecule in organic_3d_mols:
# #     print(molecule.identifier)

# # reader_formats = sorted(io.MoleculeReader.known_formats.keys())
# # writer_formats = sorted(io.MoleculeWriter.known_formats.keys())

# # print('Reader formats:\n')
# # print('\n'.join(reader_formats))

# # print('\nWriter formats:\n')
# # print('\n'.join(writer_formats))

# # how to check if an entry has smiles
# csd_reader = io.EntryReader('CSD')
# entry = csd_reader.entry('ABASAB')
# # for component in entry.molecule.components:
# #     if component.smiles:
# #         print(component.smiles)
# #         print(f'The component has smiles: {component.smiles}')
# #     else:
# #         print(component.smiles)
# #         print('The component does not have smiles')
# heavy_comp = entry.molecule.heaviest_component
# print(heavy_comp.smiles)

from openbabel import pybel as pb
import numpy as np
from hsr import fingerprint as fp


# def get_array_from_pybelmol(pybelmol):
#     # Iterate over the atoms in the molecule
#     atom_array = []
#     for atom in pybelmol.atoms:
#         # Append the coordinates and the atomic number of the atom (on each row)
#         atom_array.append([atom.coords[0], atom.coords[1], atom.coords[2], atom.atomicnum])    
#     # Centering the data
#     atom_array -= np.mean(atom_array, axis=0)
#     return atom_array

# # smiles = 'CC#N[Co+]1(C(c2ccccc2)#C1c1ccccc1)([P](C)(C)C)([P](C)(C)C)[P](C)(C)C' #15_ACNCOB10
# # smiles = 'CC(=O)C1=CNC2=C3C4=C1[Fe]234(C#O)(C#O)C#O' # 17_ACAZFE
# # smiles = 'C[Si](C)(C)C(=C=[Cr]12345(C#O)(C#O)c6c1c2c3c4c56)[Si](C)(C)C' # 20_DAJLAC
# smiles = 'Br[La](Br)(Br)(N#CC)(N#CC)(N#CC)(N#CC)N#CC' #25_ABEVAG
# ref_mol = pb.readstring("smi", smiles)
# # ref_mol.addh()
# ref_mol.make3D()
# # ref_mol.localopt()
# ref_mol_array = get_array_from_pybelmol(ref_mol)

# ref_mol_fp = fp.generate_fingerprint_from_data(ref_mol_array)
# print(ref_mol_fp)

# import os
# import glob
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from hsr import fingerprint as fp
# from hsr.utils import PROTON_FEATURES

# # Directory containing target files
# files_dir = './targets/'

# # Dictionary to store CSD entry and corresponding SMILES
# target_molecules = {}

# # Retrieve SMILES and CSD entries from all N_*.txt files
# for i in range(1, 31):
#     file_pattern = os.path.join(files_dir, f"{i}_*.txt")
#     matching_files = glob.glob(file_pattern)

#     for file_path in matching_files:
#         with open(file_path, 'r') as f:
#             for line in f:
#                 parts = line.strip().split(':')
#                 if len(parts) == 2:
#                     csd_entry, smiles = parts[0].strip(), parts[1].strip()
#                     target_molecules[csd_entry] = smiles  # Store in dictionary

# # Process molecules and track failures
# failed_entries = []

# for csd_entry, smiles in target_molecules.items():
#     try:
#         mol = Chem.MolFromSmiles(smiles)
#         mol = Chem.AddHs(mol)
#         AllChem.EmbedMolecule(mol)
#         AllChem.MMFFOptimizeMolecule(mol)
#         ref_mol_array = fp.generate_fingerprint_from_molecule(mol, PROTON_FEATURES)
#     except Exception as e:
#         failed_entries.append((csd_entry, smiles, str(e)))

# # Print failures, if any
# if failed_entries:
#     print("\n❌ The following entries failed RDKit processing:")
#     for entry, smiles, error in failed_entries:
#         print(f"CSD Entry: {entry} | SMILES: {smiles}")
# else:
#     print("\n✅ All molecules processed successfully!")

# from ccdc import io
# from ccdc.molecule import Molecule
# from ccdc.io import MoleculeWriter
# from ccdc.diagram import DiagramGenerator
# from hsr import fingerprint as fp
# from hsr import similarity as sim
# from hsr.utils import PROTON_FEATURES
# import numpy as np
# import matplotlib.pyplot as plt
# import math
# from PIL import Image, ImageDraw, ImageFont
# import matplotlib
# import matplotlib.image as mpimg

# # Force Matplotlib to use 'Agg' (prevents Qt issues)
# matplotlib.use("Agg")

# def get_array_from_ccdcmol(ccdcmol):
#     # Iterate over the atoms in the molecule
#     atom_array = []
#     for atom in ccdcmol.atoms:
#         # Append the coordinates and the atomic number of the atom (on each row)
#         atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)])    
#     # Centering the data
#     atom_array -= np.mean(atom_array, axis=0)
#     return atom_array

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

# def visualize_molecule_ccdc(molecule):
#     """
#     Generates and displays a 2D diagram of a CCDC molecule with atom labels.
#     Uses PIL instead of Qt to avoid fontconfig errors.
#     """
#     diagram_generator = DiagramGenerator()

#     # Force PIL rendering to avoid Qt issues
#     diagram_generator.settings.return_type = 'PIL'

#     # Customize diagram settings
#     diagram_generator.settings.font_size = 12
#     diagram_generator.settings.line_width = 1.6
#     diagram_generator.settings.image_width = 500
#     diagram_generator.settings.image_height = 500
#     diagram_generator.settings.shrink_symbols = False  # Ensure labels are readable

#     # Label all atoms
#     atom_labels = list(molecule.atoms)  # Pass all atoms to be labeled

#     # Generate image with atom labels
#     img = diagram_generator.image(molecule, label_atoms=atom_labels)
#     img.show()


# def save_molecule_image(molecule, filename):
#     """
#     Generates and saves a 2D diagram of a molecule with atom labels.

#     :param molecule: CCDC molecule instance
#     :param filename: Filename to save the image (PNG format)
#     """
#     diagram_generator = DiagramGenerator()

#     # Force PIL rendering to avoid Qt dependency
#     diagram_generator.settings.return_type = "PIL"

#     # Customize diagram settings
#     diagram_generator.settings.font_size = 12
#     diagram_generator.settings.line_width = 1.6
#     diagram_generator.settings.image_width = 500
#     diagram_generator.settings.image_height = 500
#     diagram_generator.settings.shrink_symbols = False  # Ensure labels are readable

#     try:
#         img = diagram_generator.image(molecule, label_atoms=list(molecule.atoms))

#         if img is not None:
#             img = img.convert("RGB")  # Convert to a format that supports saving
#             img.save(filename, format="PNG")  # Save as PNG
#             print(f"✅ Image saved: {filename}")
#         else:
#             print(f"❌ Error: Image generation failed for {filename}")

#     except Exception as e:
#         print(f"❌ Exception during image generation: {e}")
    
# def create_molecule_grid(image_filenames, output_filename, columns=4, label_font_size=20, fragment_scale=0.02):
#     """
#     Combines multiple molecule images into a single image arranged in a grid with labels above.

#     :param image_filenames: List of image file paths (first is the original molecule).
#     :param output_filename: Output file name for the combined image.
#     :param columns: Number of columns in the grid (default: 4).
#     :param label_font_size: Font size for molecule labels.
#     :param fragment_scale: Scaling factor for fragment images (default: 0.5 for half-size).
#     """
#     images = [Image.open(f) for f in image_filenames]

#     # Determine original and fragment sizes
#     original_width, original_height = images[0].size
#     fragment_width = int(original_width * fragment_scale)
#     fragment_height = int(original_height * fragment_scale)

#     num_fragments = len(images) - 1  # Exclude the original molecule
#     rows = math.ceil(num_fragments / columns) + 1  # +1 for original molecule row

#     # Space for text above each molecule
#     title_height = label_font_size + 10

#     # Calculate final grid dimensions
#     grid_width = columns * fragment_width
#     grid_height = (rows * (fragment_height + title_height))  # Extra space above images for labels

#     # Create blank grid canvas
#     grid_image = Image.new("RGB", (grid_width, grid_height), "white")

#     # Load font (fallback to default if unavailable)
#     try:
#         font = ImageFont.truetype("arial.ttf", label_font_size)
#     except IOError:
#         font = ImageFont.load_default()

#     draw = ImageDraw.Draw(grid_image)

#     # Place original molecule (centered in first row)
#     x_original = (grid_width - original_width) // 2
#     y_original = title_height  # Move down after title space
#     draw.text((x_original + original_width // 4, y_original - title_height), "Original Molecule", fill="black", font=font)
#     grid_image.paste(images[0], (x_original, y_original))

#     # Resize and place fragments in the grid
#     for i, img in enumerate(images[1:], start=1):
#         img = img.resize((fragment_width, fragment_height), Image.ANTIALIAS)

#         row = (i - 1) // columns + 1  # Start from second row
#         col = (i - 1) % columns
#         x_offset = col * fragment_width
#         y_offset = (row * (fragment_height + title_height))  # Space above for label

#         # Draw title ABOVE the fragment
#         draw.text((x_offset + fragment_width // 4, y_offset - title_height), f"Fragment {i}", fill="black", font=font)

#         # Paste fragment below the title
#         grid_image.paste(img, (x_offset, y_offset))

#     # Save and display the final image
#     grid_image.save(output_filename)
#     grid_image.show()

# def create_molecule_grid_matplotlib(image_filenames, output_filename, columns=4):
#     """
#     Uses Matplotlib to generate a grid of molecule images with labels ABOVE each image.
    
#     :param image_filenames: List of image file paths (first is the original molecule).
#     :param output_filename: Output file name for the combined image.
#     :param columns: Number of columns in the grid (default: 4).
#     """
#     num_images = len(image_filenames)
#     rows = math.ceil(num_images / columns)  # Auto-adjust rows

#     fig, axes = plt.subplots(rows, columns, figsize=(columns * 2.5, rows * 2.5))

#     # If there's only one row, make sure axes is iterable
#     if rows == 1:
#         axes = [axes]

#     # Load and display each image
#     for i, ax in enumerate(axes.flat):
#         if i < num_images:
#             img = mpimg.imread(image_filenames[i])  # Load image
#             ax.imshow(img)
#             ax.set_title(f"Fragment {i}" if i > 0 else "Original Molecule", fontsize=10)
#         ax.axis("off")  # Hide axes

#     plt.tight_layout()
#     plt.savefig(output_filename, dpi=300)  # Save final image
#     plt.show()

# def create_fragment(central_atom):
#     """
#     Creates a molecular fragment centered on the given atom, including its first neighbors.

#     :param central_atom: A ccdc.molecule.Atom instance
#     :return: A ccdc.molecule.Molecule fragment
#     """
#     # Create an empty molecule for the fragment
#     fragment = Molecule(identifier=f"{central_atom.label}_frag")

#     # Add central atom to the fragment
#     fragment.add_atom(central_atom)

#     # Add first neighbors
#     [fragment.add_atom(neighbour) for neighbour in central_atom.neighbours]

#     return fragment

# def get_fragments(molecule):
#     """
#     Generates fragments for all atoms in the given molecule.
    
#     :param molecule: A ccdc.molecule.Molecule instance
#     :return: List of ccdc.molecule.Molecule fragments
#     """
#     fragments = [create_fragment(atom) for atom in molecule.atoms]
#     return fragments

# def print_molecule_atoms(molecule, title="Molecule"):
#     """
#     Prints atomic symbols and coordinates for a given molecule.

#     :param molecule: A ccdc.molecule.Molecule instance
#     :param title: Title to display for clarity
#     """
#     print(f"\n=== {title} ===")
#     for atom in molecule.atoms:
#         coords = atom.coordinates
#         print(f"{atom.label} ({atom.atomic_symbol}): ({coords.x:.3f}, {coords.y:.3f}, {coords.z:.3f})")

# Example Usage
csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry('OFOWIS')
mol = component_of_interest(entry.molecule)

# Save all fragments to an SDF file
with MoleculeWriter("OFOWIS.sdf") as writer:
    writer.write(mol)


# # Load from CSD
# csd_reader = io.EntryReader('CSD')
# entry = csd_reader.entry('ABAHIW')
# molecule = component_of_interest(entry.molecule)

# # Save the original molecule as an image
# save_molecule_image(molecule, "original_molecule.png")

# # Print the selected component's atoms
# print_molecule_atoms(molecule, title="Component of Interest")

# # Generate fragments for the selected component
# fragments = get_fragments(molecule)

# # # Print each fragment from the component of interest
# # for i, frag in enumerate(fragments):
# #     print_molecule_atoms(frag, title=f"Fragment {i+1}")

# # Save each fragment as an image
# fragment_image_files = []
# for i, frag in enumerate(fragments):
#     filename = f"fragment_{i+1}.png"
#     save_molecule_image(frag, filename)
#     fragment_image_files.append(filename)
    

# # Fragments analysis: Get unique fragments
# fragment_fps = []
# for fragment in fragments:
#     # get HSR fingerprint
#     fragment_array = get_array_from_ccdcmol(fragment)
#     fragment_fp = fp.generate_fingerprint_from_data(fragment_array)
#     fragment_fps.append(fragment_fp)
    
# # construct and print the matrix of similarities between fragments
# similarity_matrix = np.zeros((len(fragments), len(fragments)))
# for i, fp1 in enumerate(fragment_fps):
#     for j, fp2 in enumerate(fragment_fps):
#         similarity_matrix[i, j] = sim.compute_similarity_score(fp1, fp2)
        
# # Determine column width (adjust for readability)
# col_width = 5  # Adjust width as needed for spacing

# # Print header row (column numbers) with correct alignment
# print("\nSimilarity Matrix:")

# header = ["{:>{w}}".format(i+1, w=col_width) for i in range(len(fragments))]
# print(" " * col_width + " ".join(header))  # Add leading space for row labels

# # Print rows with formatted values
# for i, row in enumerate(similarity_matrix):
#     formatted_row = ["{:>{w}.2f}".format(value, w=col_width) for value in row]  # Right-align numbers
#     print("{:>{w}} |".format(i+1, w=col_width - 1) + " ".join(formatted_row))  # Add separator for row index


# # Create a grid image with all molecule images
# # create_molecule_grid(["original_molecule.png"] + fragment_image_files, "molecule_grid.png", columns=4)
# create_molecule_grid_matplotlib(["original_molecule.png"] + fragment_image_files, "molecule_grid_matplotlib.png", columns=4)
