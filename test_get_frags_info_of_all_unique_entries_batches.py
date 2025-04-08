import os
import csv
import time
import numpy as np
from collections import Counter
from multiprocessing import Pool, cpu_count
from ccdc import io
from ccdc.molecule import Molecule
from hsr import fingerprint as fp

# -------------------------------
# Helper Functions
# -------------------------------

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

def create_fragment(central_atom):
    frag = Molecule(identifier=f"{central_atom.label}_frag")
    atom_map = {central_atom: frag.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = frag.add_atom(neighbor)
    for bond in central_atom.bonds:
        a1, a2 = bond.atoms
        if a1 in atom_map and a2 in atom_map:
            try:
                frag.add_bond(bond.bond_type, atom_map[a1], atom_map[a2])
            except:
                pass
    return frag

def get_fragments(mol):
    return [create_fragment(atom) for atom in mol.atoms]

def get_array_from_ccdcmol(ccdcmol):
    coords = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return coords - coords.mean(axis=0)

def formula_signature(fragment):
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, formula)

def generate_fp_data(fragment):
    fp_data = fp.generate_fingerprint_from_data(get_array_from_ccdcmol(fragment))
    if isinstance(fp_data, np.ndarray):
        fp_data = fp_data.tolist()

    return {
        "sdf": fragment.to_string("sdf"),
        "fp": fp_data,
        "formula": formula_signature(fragment),
        "n_atoms": len(fragment.atoms),
        "central_atom": fragment.atoms[0].atomic_symbol,
    }

def process_population_entry(entry_id):
    try:
        reader = io.EntryReader("CSD")
        entry = reader.entry(entry_id)
        mol = component_of_interest(entry.molecule)
        if mol is None:
            return []

        fragments = get_fragments(mol)
        return [{
            **generate_fp_data(frag),
            "entry_id": entry_id,
            "frag_idx": i
        } for i, frag in enumerate(fragments)]

    except Exception as e:
        print(f"❌ Error processing {entry_id}: {e}")
        return []

def extract_unique_entries_from_populations(folder_path):
    unique_entries = set()
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            print(f"📄 Reading: {filename}")
            with open(os.path.join(folder_path, filename), 'r', encoding='utf-8') as file:
                for line in file:
                    if line.strip():
                        unique_entries.add(line.strip().split()[0])
    print(f"✅ Total unique entries: {len(unique_entries)}")
    return list(unique_entries)

def write_fragments_to_csv(data, output_csv, header=False):
    keys = ["entry_id", "frag_idx", "sdf", "fp", "formula", "n_atoms", "central_atom"]
    with open(output_csv, 'a', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        if header:
            writer.writeheader()
        for row in data:
            writer.writerow(row)

def chunk_list(data, chunk_size):
    for i in range(0, len(data), chunk_size):
        yield data[i:i + chunk_size]

# -------------------------------
# Main
# -------------------------------

if __name__ == "__main__":
    start = time.time()
    folder = "targets/init_population_protons_2/"
    output_csv = "all_fragments_data.csv"
    batch_size = 5000

    if os.path.exists(output_csv):
        os.remove(output_csv)

    entries = extract_unique_entries_from_populations(folder)
    n_proc = min(cpu_count(), 8)
    print(f"⚙️ Using {n_proc} processes")

    for i, batch in enumerate(chunk_list(entries, batch_size), start=1):
        print(f"🔄 Processing batch {i} / {len(entries) // batch_size + 1}")
        results = []
        with Pool(n_proc) as pool:
            for res in pool.imap_unordered(process_population_entry, batch):
                results.extend(res)

        write_fragments_to_csv(results, output_csv, header=(i == 1))
        print(f"✅ Saved batch {i} ({len(results)} fragments)")

    print(f"🏁 Done! Total time: {time.time() - start:.2f} seconds")
