from ccdc.io import EntryReader
from ccdc.diagram import DiagramGenerator
import os

def component_of_interest(molecule):
    components = molecule.components
    if not components:
        return None
    props = [{
        "component": c,
        "is_organometallic": c.is_organometallic,
        "mw": sum(atom.atomic_weight for atom in c.atoms),
        "atom_count": len(c.atoms)
    } for c in components]
    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["atom_count"])
    for prop in props:
        if sum([
            prop["is_organometallic"],
            prop["component"] == heaviest["component"],
            prop["component"] == most_atoms["component"]
        ]) >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None

TARGETS = [
    'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE', 'NIWPUE01',
    'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG', 'ABOBUP', 'XIDTOW',
    'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ', 'ADUPAS', 'DAJLAC', 'OFOWIS',
    'CATSUL', 'HESMUQ01', 'GUDQOL', 'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA',
    'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM'
]

# Make sure the folder exists
output_dir = "Targets_diagrams_test"
os.makedirs(output_dir, exist_ok=True)

entry_reader = EntryReader('CSD')
diagram_generator = DiagramGenerator()

## SETTINGS ##
diagram_generator.settings.font_size = 10
diagram_generator.settings.line_width = 1
diagram_generator.settings.image_width = 250
diagram_generator.settings.image_height = 250
diagram_generator.settings.return_type = 'SVG'

for i, target in enumerate(TARGETS, start=1):
    try:
        entry = entry_reader.entry(target)
        mol = entry.molecule
        comp = component_of_interest(mol)
        if comp is None:
            print(f"[{i}] {target}: No suitable component found, skipping.")
            continue
        img = diagram_generator.image(comp)
        filename = os.path.join(output_dir, f"{i}_{target}.svg")
        with open(filename, "w") as f:
            f.write(img)
        print(f"[{i}] {target}: Diagram saved.")
    except Exception as e:
        print(f"[{i}] {target}: Error ({e})")