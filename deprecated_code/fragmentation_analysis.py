
import os
os.environ["QT_QPA_PLATFORM"] = "xcb"  # Ensure X11 is used, not Wayland

# Set Matplotlib backend BEFORE importing pyplot
import matplotlib
matplotlib.use("TkAgg")  # Change to "TkAgg" if you still get issues

import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.image as mpimg
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc.diagram import DiagramGenerator
from hsr import fingerprint as fp
from hsr import similarity as sim
from hsr.utils import PROTON_FEATURES
from PIL import Image, ImageDraw, ImageFont


def get_array_from_ccdcmol(ccdcmol):
    # Iterate over the atoms in the molecule
    atom_array = []
    for atom in ccdcmol.atoms:
        # Append the coordinates and the atomic number of the atom (on each row)
        atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)])    
    # Centering the data
    atom_array -= np.mean(atom_array, axis=0)
    return atom_array

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

def visualize_molecule_ccdc(molecule):
    """
    Generates and displays a 2D diagram of a CCDC molecule with atom labels.
    Uses PIL instead of Qt to avoid fontconfig errors.
    """
    diagram_generator = DiagramGenerator()

    # Force PIL rendering to avoid Qt issues
    diagram_generator.settings.return_type = 'PIL'

    # Customize diagram settings
    diagram_generator.settings.font_size = 12
    diagram_generator.settings.line_width = 1.6
    diagram_generator.settings.image_width = 500
    diagram_generator.settings.image_height = 500
    diagram_generator.settings.shrink_symbols = False  # Ensure labels are readable

    # Label all atoms
    atom_labels = list(molecule.atoms)  # Pass all atoms to be labeled

    # Generate image with atom labels
    img = diagram_generator.image(molecule, label_atoms=atom_labels)
    img.show()
    

def save_molecule_image(molecule, filename):
    """
    Generates and saves a 2D diagram of a molecule with atom labels.

    :param molecule: CCDC molecule instance
    :param filename: Filename to save the image (PNG format)
    """
    diagram_generator = DiagramGenerator()

    # Force PIL rendering to avoid Qt dependency
    diagram_generator.settings.return_type = "PIL"

    # Customize diagram settings
    diagram_generator.settings.font_size = 17
    diagram_generator.settings.line_width = 1.6
    diagram_generator.settings.image_width = 500
    diagram_generator.settings.image_height = 500
    diagram_generator.settings.shrink_symbols = False  # Ensure labels are readable

    try:
        img = diagram_generator.image(molecule, label_atoms=list(molecule.atoms))

        if img is not None:
            img = img.convert("RGB")  # Convert to a format that supports saving
            img.save(filename, format="PNG")  # Save as PNG
            print(f"✅ Image saved: {filename}")
        else:
            print(f"❌ Error: Image generation failed for {filename}")

    except Exception as e:
        print(f"❌ Exception during image generation: {e}")

def create_molecule_grid(image_filenames, output_filename, columns=4, label_font_size=20):
    """
    Combines multiple molecule images into a single image arranged in a grid.

    :param image_filenames: List of image file paths (first is the original molecule).
    :param output_filename: Output file name for the combined image.
    :param columns: Number of columns in the grid (default: 4).
    :param label_font_size: Font size for molecule labels.
    """
    images = [Image.open(f) for f in image_filenames]

    # Determine grid size
    img_width, img_height = images[0].size
    num_fragments = len(images) - 1  # Excluding the original molecule
    rows = math.ceil(num_fragments / columns) + 1  # +1 for original molecule row

    # Create blank canvas for the grid
    grid_width = columns * img_width
    grid_height = rows * img_height
    grid_image = Image.new("RGB", (grid_width, grid_height), "white")

    # Load a basic font for labeling (adjust if necessary)
    try:
        font = ImageFont.truetype("arial.ttf", label_font_size)
    except IOError:
        font = ImageFont.load_default()

    draw = ImageDraw.Draw(grid_image)

    # Paste original molecule at the top (centered)
    x_original = (grid_width - img_width) // 2
    grid_image.paste(images[0], (x_original, 0))
    draw.text((x_original + img_width // 4, img_height - 40), "Original Molecule", fill="black", font=font)

    # Paste fragments below in a grid
    for i, img in enumerate(images[1:], start=1):  
        row = (i - 1) // columns + 1  # Start from second row
        col = (i - 1) % columns
        x_offset = col * img_width
        y_offset = row * img_height
        grid_image.paste(img, (x_offset, y_offset))

        # Label each fragment with its index
        draw.text((x_offset + img_width // 4, y_offset + img_height - 40), f"Fragment {i}", fill="black", font=font)

    # Save and show the combined image
    grid_image.save(output_filename)
    grid_image.show()
    
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.image as mpimg

# Ensure proper backend for interactive Matplotlib
os.environ["QT_QPA_PLATFORM"] = "xcb"
matplotlib.use("TkAgg")


def create_molecule_grid_matplotlib(image_filenames, output_filename, columns=4, 
                                    original_scale_factor=1.5, fragment_scale_factor=1.2):
    """
    Uses Matplotlib to generate a grid where:
    - The original molecule occupies the first row alone (centered).
    - The fragments are arranged compactly in rows below.
    - The scale factors for the original molecule and fragments are separate.

    :param image_filenames: List of image file paths (first is the original molecule).
    :param output_filename: Output file name for the combined image.
    :param columns: Number of columns for the fragment grid (default: 4).
    :param original_scale_factor: Scale factor for the original molecule (default: 1.5).
    :param fragment_scale_factor: Scale factor for the fragments (default: 1.2).
    """
    num_images = len(image_filenames)
    num_fragments = num_images - 1  # Exclude the original molecule

    # Determine the number of rows for the fragments
    fragment_rows = math.ceil(num_fragments / columns)
    total_rows = fragment_rows + 1  # One extra row for the original molecule

    # Dynamically scale figure size based on image count
    fig_width = columns * 3  # Increase cell width
    fig_height = total_rows * 3  # Increase cell height

    fig, axes = plt.subplots(total_rows, columns, figsize=(fig_width, fig_height))

    # Ensure axes is always iterable
    axes = np.array(axes).reshape(total_rows, columns)

    def load_and_resize_image(image_path, scale_factor):
        """Load an image and scale it while preserving aspect ratio."""
        img = Image.open(image_path)
        width, height = img.size

        # Scale the image size
        new_width = int(width * scale_factor)
        new_height = int(height * scale_factor)
        img = img.resize((new_width, new_height), Image.LANCZOS)

        return np.array(img)

    # Load and display the original molecule centered in the first row
    original_img = load_and_resize_image(image_filenames[0], original_scale_factor)
    for col in range(columns):  # Span across all columns
        if col == 0:
            axes[0, 0].imshow(original_img)
            axes[0, 0].set_title("Original Molecule", fontsize=14, fontweight="bold")
            axes[0, 0].axis("off")
        else:
            axes[0, col].axis("off")  # Hide empty spaces

    # Load and display each fragment in the remaining space
    fragment_index = 0
    for row in range(1, total_rows):  # Start from second row
        for col in range(columns):
            if fragment_index < num_fragments:
                img = load_and_resize_image(image_filenames[fragment_index + 1], fragment_scale_factor)
                axes[row, col].imshow(img)
                axes[row, col].set_title(f"Fragment {fragment_index + 1}", fontsize=12)
                axes[row, col].axis("off")
                fragment_index += 1
            else:
                axes[row, col].axis("off")  # Hide empty subplots

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)  # Save final image
    print(f"✅ Image grid saved as {output_filename}")

    # Display the image properly
    plt.show()

def create_fragment(central_atom):
    """
    Creates a molecular fragment centered on the given atom, including its first neighbors and bonds.

    :param central_atom: A ccdc.molecule.Atom instance
    :return: A ccdc.molecule.Molecule fragment with bond information
    """
    # Create an empty molecule for the fragment
    fragment = Molecule(identifier=f"{central_atom.label}_frag")

    # Mapping of original atoms to new fragment atoms
    atom_map = {}

    # Add central atom to the fragment
    new_central = fragment.add_atom(central_atom)
    atom_map[central_atom] = new_central

    # Add first neighbors and store them in the atom map
    for neighbor in central_atom.neighbours:
        new_atom = fragment.add_atom(neighbor)
        atom_map[neighbor] = new_atom

    # Add bonds from the original molecule
    for bond in central_atom.bonds:
        atom1, atom2 = bond.atoms

        # Ensure both atoms exist in the fragment before adding the bond
        if atom1 in atom_map and atom2 in atom_map:
            try:
                fragment.add_bond(bond.bond_type, atom_map[atom1], atom_map[atom2])
            except Exception as e:
                print(f"❌ Error adding bond between {atom1.label} and {atom2.label}: {e}")

    return fragment


def get_fragments(molecule):
    """
    Generates fragments for all atoms in the given molecule.
    
    :param molecule: A ccdc.molecule.Molecule instance
    :return: List of ccdc.molecule.Molecule fragments
    """
    fragments = [create_fragment(atom) for atom in molecule.atoms]
    return fragments

def print_molecule_atoms(molecule, title="Molecule"):
    """
    Prints atomic symbols and coordinates for a given molecule.

    :param molecule: A ccdc.molecule.Molecule instance
    :param title: Title to display for clarity
    """
    print(f"\n=== {title} ===")
    for atom in molecule.atoms:
        coords = atom.coordinates
        print(f"{atom.label} ({atom.atomic_symbol}): ({coords.x:.3f}, {coords.y:.3f}, {coords.z:.3f})")
        
        
def find_unique_fragments(fragments, similarity_matrix, threshold=0.9990):
    """
    Identifies unique fragments based on the similarity matrix and threshold.

    :param fragments: List of CCDC molecules representing fragments.
    :param similarity_matrix: NxN numpy array with similarity scores.
    :param threshold: Float threshold for uniqueness (default = 0.9990).
    :return: List of unique fragment indices.
    """
    num_fragments = len(fragments)
    unique_indices = []
    seen = set()

    for i in range(num_fragments):
        if i in seen:
            continue  # Skip already grouped fragments

        unique_indices.append(i)
        for j in range(i + 1, num_fragments):
            if similarity_matrix[i, j] >= threshold:
                seen.add(j)  # Mark similar fragments as "seen"

    return unique_indices

# Example Usage
csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry('ABAHIW')
molecule = component_of_interest(entry.molecule)

# Save the original molecule as an image
# save_molecule_image(molecule, "original_molecule.png")


# Print the original molecule
print_molecule_atoms(molecule, title="Original Molecule")

# Generate diagram
# diagram_generator = DiagramGenerator()
# img = diagram_generator.image(mol)
# img.show()

# Visualize the molecule
# visualize_molecule_ccdc(mol)

# Generate fragments
fragments = get_fragments(molecule)

# Print each fragment
for i, frag in enumerate(fragments):
    print_molecule_atoms(frag, title=f"Fragment {i+1}")

# Save all fragments to an SDF file
with MoleculeWriter("fragments.sdf") as writer:
    for frag in fragments:
        writer.write(frag)


# Print the selected component's atoms
print_molecule_atoms(molecule, title="Component of Interest")

# Generate fragments for the selected component
fragments = get_fragments(molecule)

# Save each fragment as an image
fragment_image_files = []
for i, frag in enumerate(fragments):
    filename = f"fragment_{i+1}.png"
    # save_molecule_image(frag, filename)
    fragment_image_files.append(filename)

# Fragment image files from already saved images
# fragment_image_files = [f"fragment_{i+1}.png" for i in range(len(fragments))]
    
# Fragments analysis: Get unique fragments
fragment_fps = []
for fragment in fragments:
    # get HSR fingerprint
    fragment_array = get_array_from_ccdcmol(fragment)
    fragment_fp = fp.generate_fingerprint_from_data(fragment_array)
    fragment_fps.append(fragment_fp)
    
# construct and print the matrix of similarities between fragments
similarity_matrix = np.zeros((len(fragments), len(fragments)))
for i, fp1 in enumerate(fragment_fps):
    for j, fp2 in enumerate(fragment_fps):
        similarity_matrix[i, j] = sim.compute_similarity_score(fp1, fp2)
        
# Determine column width dynamically based on decimal places
decimal_places = 4
col_width = decimal_places + 2  # Adjust width to align properly

# Print header row (column numbers) with correct alignment
print("\nSimilarity Matrix:")

header = ["{:>{w}}".format(i+1, w=col_width) for i in range(len(fragments))]
print(" " * (col_width - 1) + " ".join(header))  # Add leading space for row labels

# Print rows with formatted values
for i, row in enumerate(similarity_matrix):
    formatted_row = ["{:>{w}.{d}f}".format(value, w=col_width, d=decimal_places) for value in row]  # Align numbers
    print("{:>{w}} |".format(i+1, w=col_width - 1) + " ".join(formatted_row))  # Add separator for row index


# Create a grid image with all molecule images
# create_molecule_grid(["original_molecule.png"] + fragment_image_files, "molecule_grid.png", columns=4)
# create_molecule_grid_matplotlib(["original_molecule.png"] + fragment_image_files, "molecule_grid.png", columns=4)


# get unique fragments
print("\nUnique Fragments:")
unique_indices = find_unique_fragments(fragments, similarity_matrix, threshold=0.9990)
unique_fragments = [fragments[i] for i in unique_indices]
for i, frag in enumerate(unique_fragments):
    print_molecule_atoms(frag, title=f"Unique Fragment {i+1}")

# Save png images of unique fragments
unique_fragments_image_files = []
for i, frag in enumerate(unique_fragments):
    filename = f"unique_fragment_{i+1}.png"
    # save_molecule_image(frag, filename)
    unique_fragments_image_files.append(filename)

# Save xml files of unique fragments
# diagram_generator = DiagramGenerator()
# for i, frag in enumerate(unique_fragments):
#     filename = f"unique_fragment_{i+1}.xml"
#     diagram_generator.chemdraw_xml(molecule=frag, file_name=filename)

# Print information: Number of fragments and unique fragments
n_atoms_original = len(molecule.atoms)  
num_fragments = len(fragments)
num_unique_fragments = len(unique_fragments)
print(f"\nOriginal Molecule Atoms: {n_atoms_original}")
print(f"\nTotal Fragments: {num_fragments}")
print(f"Unique Fragments: {num_unique_fragments}")


# Visualize original molecule and unique fragments
# create_molecule_grid_matplotlib(["original_molecule.png"] + unique_fragments_image_files, "molecule_grid.png", columns=4)