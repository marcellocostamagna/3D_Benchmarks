import os
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image, ImageDraw, ImageFont
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc.diagram import DiagramGenerator
from hsr import fingerprint as fp
from hsr import similarity as sim

# Ensure X11 for Qt compatibility
os.environ["QT_QPA_PLATFORM"] = "xcb"
matplotlib.use("TkAgg")  # Ensure interactive Matplotlib backend


### ------------------ Molecule Processing Functions ------------------ ###

def get_array_from_ccdcmol(ccdcmol):
    """Extracts an array of atomic coordinates and atomic numbers."""
    atom_array = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return atom_array - np.mean(atom_array, axis=0)  # Centering data


def component_of_interest(molecule):
    """Identifies the most relevant component of a molecule based on weight and atom count."""
    components = molecule.components
    if not components:
        return None

    properties = [{
        "component": comp,
        "is_organometallic": comp.is_organometallic,
        "molecular_weight": sum(atom.atomic_weight for atom in comp.atoms),
        "atom_count": len(comp.atoms)
    } for comp in components]

    heaviest = max(properties, key=lambda x: x["molecular_weight"])
    most_atoms = max(properties, key=lambda x: x["atom_count"])

    for prop in properties:
        criteria_met = sum([
            prop["is_organometallic"],
            prop["component"] == heaviest["component"],
            prop["component"] == most_atoms["component"]
        ])
        if criteria_met >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    
    return None

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
        

def create_fragment(central_atom):
    """Creates a molecular fragment centered on a given atom, including its first neighbors and bonds."""
    fragment = Molecule(identifier=f"{central_atom.label}_frag")
    atom_map = {central_atom: fragment.add_atom(central_atom)}

    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = fragment.add_atom(neighbor)

    for bond in central_atom.bonds:
        atom1, atom2 = bond.atoms
        if atom1 in atom_map and atom2 in atom_map:
            try:
                fragment.add_bond(bond.bond_type, atom_map[atom1], atom_map[atom2])
            except Exception:
                pass  # Ignore duplicate bonds
    
    return fragment


def get_fragments(molecule):
    """Generates molecular fragments for all atoms."""
    return [create_fragment(atom) for atom in molecule.atoms]


def find_unique_fragments(fragments, similarity_matrix, threshold=0.9990):
    """Finds unique fragments by filtering out similar ones based on a similarity threshold."""
    unique_indices = []
    seen = set()

    for i in range(len(fragments)):
        if i in seen:
            continue
        unique_indices.append(i)
        seen.update(j for j in range(i + 1, len(fragments)) if similarity_matrix[i, j] >= threshold)
    
    return unique_indices


### ------------------ Visualization Functions ------------------ ###

def save_molecule_image(molecule, filename):
    """Generates and saves a 2D diagram of a molecule."""
    diagram_generator = DiagramGenerator()
    diagram_generator.settings.return_type = "PIL"
    diagram_generator.settings.font_size = 17
    diagram_generator.settings.line_width = 1.6
    diagram_generator.settings.image_width = 500
    diagram_generator.settings.image_height = 500
    diagram_generator.settings.shrink_symbols = False

    try:
        img = diagram_generator.image(molecule, label_atoms=list(molecule.atoms))
        if img:
            img.convert("RGB").save(filename, format="PNG")
            print(f"✅ Image saved: {filename}")
    except Exception as e:
        print(f"❌ Error saving image: {e}")


def create_molecule_grid_matplotlib(image_filenames, output_filename, columns=4, 
                                    original_scale_factor=1.5, fragment_scale_factor=1.2):
    """Creates a grid of molecule images with labels above each image."""
    num_images = len(image_filenames)
    num_fragments = num_images - 1
    fragment_rows = math.ceil(num_fragments / columns)
    total_rows = fragment_rows + 1  # One row for the original molecule

    fig, axes = plt.subplots(total_rows, columns, figsize=(columns * 3, total_rows * 3))
    axes = np.array(axes).reshape(total_rows, columns)

    def load_and_resize(image_path, scale_factor):
        img = Image.open(image_path)
        new_size = (int(img.width * scale_factor), int(img.height * scale_factor))
        return np.array(img.resize(new_size, Image.LANCZOS))

    original_img = load_and_resize(image_filenames[0], original_scale_factor)
    axes[0, 0].imshow(original_img)
    axes[0, 0].set_title("Original Molecule", fontsize=14, fontweight="bold")
    axes[0, 0].axis("off")
    for col in range(1, columns):
        axes[0, col].axis("off")  # Hide empty spaces

    fragment_index = 0
    for row in range(1, total_rows):
        for col in range(columns):
            if fragment_index < num_fragments:
                img = load_and_resize(image_filenames[fragment_index + 1], fragment_scale_factor)
                axes[row, col].imshow(img)
                axes[row, col].set_title(f"Fragment {fragment_index + 1}", fontsize=12)
                axes[row, col].axis("off")
                fragment_index += 1
            else:
                axes[row, col].axis("off")

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"✅ Image grid saved as {output_filename}")
    plt.show()


### ------------------ Main Processing ------------------ ###

# Load molecule and select the most relevant component
csd_reader = io.EntryReader('CSD')
entry = csd_reader.entry('ACOVUL')
molecule = component_of_interest(entry.molecule)

# Save the original molecule as an image
save_molecule_image(molecule, "original_molecule.png")
print_molecule_atoms(molecule, title="Original Molecule")

# Generate and save fragments
fragments = get_fragments(molecule)
fragment_image_files = []
# with MoleculeWriter("fragments.sdf") as writer:
#     for i, frag in enumerate(fragments):
#         filename = f"fragment_{i+1}.png"
#         # save_molecule_image(frag, filename)
#         fragment_image_files.append(filename)
#         writer.write(frag)

# Compute similarity matrix
fragment_fps = [fp.generate_fingerprint_from_data(get_array_from_ccdcmol(frag)) for frag in fragments]
similarity_matrix = np.array([[sim.compute_similarity_score(fp1, fp2) for fp2 in fragment_fps] for fp1 in fragment_fps])

# Find and save unique fragments
unique_indices = find_unique_fragments(fragments, similarity_matrix, threshold=0.9990)
unique_fragments = [fragments[i] for i in unique_indices]
unique_fragment_image_files = []
for i, frag in enumerate(unique_fragments):
    filename = f"unique_fragment_{i+1}.png"
    save_molecule_image(frag, filename)
    unique_fragment_image_files.append(filename)
    
unique_fragment_image_files = []    
with MoleculeWriter("unique_fragments.sdf") as writer:
    for i, frag in enumerate(unique_fragments):
        filename = f"unique_fragment_{i+1}.png"
        unique_fragment_image_files.append(filename)
        writer.write(frag)
        
# Print summary
print(f"\nOriginal Molecule Atoms: {len(molecule.atoms)}")
print(f"Total Fragments: {len(fragments)}")
print(f"Unique Fragments: {len(unique_fragments)}")

# Visualize original molecule and unique fragments
create_molecule_grid_matplotlib(["original_molecule.png"] + unique_fragment_image_files, "molecule_grid.png", columns=4)
