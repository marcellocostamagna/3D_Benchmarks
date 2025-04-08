#!/usr/bin/env python3

import os
import time
import numpy as np
import pandas as pd
import duckdb
from collections import Counter
from multiprocessing import Pool, cpu_count
from ccdc import io
from ccdc.molecule import Molecule
from hsr import fingerprint as fp
from functools import partial

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

# Function to generate an array from a molecule
# EXTRA FEATURE: PROTONS
def get_p_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

# EXTRA FEATURES: PROTONS & FORMAL CHARGES
def get_p_q_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number), atom.formal_charge])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

def formula_signature(fragment):
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula_str = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, formula_str)

def generate_fp_data(fragment, with_charges=False):
    # Choose the right array builder
    if with_charges:
        arr = get_p_q_array_from_ccdcmol(fragment)
    else:
        arr = get_p_array_from_ccdcmol(fragment)

    # Generate fingerprint
    fp_data = fp.generate_fingerprint_from_data(arr)
    if isinstance(fp_data, np.ndarray):
        fp_data = fp_data.tolist()

    # Get formula info
    central, brute_formula = formula_signature(fragment)
    return {
        "sdf": fragment.to_string("sdf"),
        "fp": fp_data,
        "central_atom": central,
        "n_atoms": len(fragment.atoms),
        "formula_str": brute_formula,
    }


def process_population_entry(entry_id, with_charges_ids=None):
    try:
        with_charges = entry_id in with_charges_ids if with_charges_ids else False
        reader = io.EntryReader("CSD")
        entry = reader.entry(entry_id)
        mol = component_of_interest(entry.molecule)
        if mol is None:
            return []

        fragments = get_fragments(mol)
        results = []
        for i, frag in enumerate(fragments):
            frag_data = generate_fp_data(frag, with_charges=with_charges)
            frag_data["entry_id"] = entry_id
            frag_data["frag_idx"] = i
            results.append(frag_data)
        return results

    except Exception as e:
        print(f"❌ Error processing {entry_id}: {e}")
        return []


def extract_entries_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    identifiers = df["Identifier"].dropna().unique().tolist()
    print(f"✅ Total unique entries: {len(identifiers)}")
    return identifiers

def chunk_list(data, chunk_size):
    for i in range(0, len(data), chunk_size):
        yield data[i:i + chunk_size]

def main():
    start = time.time()
    input_csv = "targets/viable_structures.csv"
    duckdb_path = "all_fragments.duckdb"
    batch_size = 5000
    with_charges_ids = ['ABAYAF', 'RULJAM']

    # Remove existing DB if you want a fresh start
    if os.path.exists(duckdb_path):
        os.remove(duckdb_path)

    entries = extract_entries_from_csv(input_csv)
    n_proc = min(cpu_count(), 8)
    print(f"⚙️ Using {n_proc} processes")

    con = duckdb.connect(duckdb_path)
    # Make sure column order in DuckDB: entry_id, frag_idx, sdf, fp, central_atom, n_atoms, formula_str
    con.execute("""
        CREATE TABLE fragments (
            entry_id TEXT,
            frag_idx INT,
            sdf TEXT,
            fp TEXT,
            central_atom TEXT,
            n_atoms INT,
            formula_str TEXT
        )
    """)

    total_fragments_written = 0
    batch_number = 0

    for batch in chunk_list(entries, batch_size):
        batch_number += 1
        print(f"🔄 Processing batch {batch_number} / {len(entries)//batch_size + 1}")

        # Parallel generation
        results = []
        with Pool(n_proc) as pool:
            func = partial(process_population_entry, with_charges_ids=with_charges_ids)
            for frag_list in pool.imap_unordered(func, batch):
                results.extend(frag_list)

        if not results:
            print(f"✅ No fragments in batch", batch_number)
            continue

        # Build DataFrame
        df = pd.DataFrame(results)
        # Force the DataFrame columns to match the EXACT order in the DuckDB table
        # If any column is missing, fill with default or drop columns not needed
        required_cols = [
            "entry_id", 
            "frag_idx", 
            "sdf", 
            "fp", 
            "central_atom", 
            "n_atoms", 
            "formula_str"
        ]
        # ensure all required columns exist (some might be missing if there's an error)
        for col in required_cols:
            if col not in df.columns:
                df[col] = None  # or a default

        df = df[required_cols]
        # Convert 'fp' from list -> string
        df["fp"] = df["fp"].apply(lambda x: str(x))

        if df.empty:
            print(f"✅ DataFrame for batch {batch_number} is empty; skipping insert.")
            continue

        # Register DF, create a temp table by selecting columns in the correct order
        con.register("df_table", df)
        temp_name = f"tmp_fragments_{batch_number}"
        con.execute(f"""
            CREATE TEMPORARY TABLE {temp_name} AS
            SELECT 
                entry_id,
                frag_idx,
                sdf,
                fp,
                central_atom,
                n_atoms,
                formula_str
            FROM df_table
        """)

        # Insert into main table in the exact same order
        con.execute(f"""
            INSERT INTO fragments
            SELECT
                entry_id,
                frag_idx,
                sdf,
                fp,
                central_atom,
                n_atoms,
                formula_str
            FROM {temp_name}
        """)
        con.execute(f"DROP TABLE {temp_name}")

        batch_frag_count = len(df)
        total_fragments_written += batch_frag_count
        print(f"✅ Saved batch {batch_number} ({batch_frag_count} fragments)")

    # Create index if desired
    con.execute("""
        CREATE INDEX idx_formula 
        ON fragments(central_atom, n_atoms, formula_str)
    """)
    con.close()

    elapsed = time.time() - start
    print(f"\n🏁 Done! Inserted {total_fragments_written} fragments into {duckdb_path}.")
    print(f"Total time: {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()
