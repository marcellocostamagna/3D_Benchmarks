"""
This script extracts Connected Atom Environments (CAEs) from molecular entries listed in a CSV file,
computes their fingerprints, and stores their structural and fingerprint data in a DuckDB database ('all_caes.duckdb').
"""

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
from collections import defaultdict

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

def create_cae(central_atom):
    cae = Molecule(identifier=f"{central_atom.label}_cae")
    atom_map = {central_atom: cae.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = cae.add_atom(neighbor)
    for bond in central_atom.bonds:
        a1, a2 = bond.atoms
        if a1 in atom_map and a2 in atom_map:
            try:
                cae.add_bond(bond.bond_type, atom_map[a1], atom_map[a2])
            except:
                pass
    return cae

def get_caes(mol):
    return [create_cae(atom) for atom in mol.atoms]

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

def formula_signature(cae):
    atoms = cae.atoms
    central = atoms[0].atomic_symbol
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula_str = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, formula_str)

def generate_fp_data(cae, with_charges=False):
    arr_p = get_p_array_from_ccdcmol(cae)
    fp_p = fp.generate_fingerprint_from_data(arr_p)
    if isinstance(fp_p, np.ndarray):
        fp_p = fp_p.tolist()

    fp_pq = None
    if with_charges:
        arr_pq = get_p_q_array_from_ccdcmol(cae)
        fp_pq = fp.generate_fingerprint_from_data(arr_pq)
        if isinstance(fp_pq, np.ndarray):
            fp_pq = fp_pq.tolist()

    central, brute_formula = formula_signature(cae)
    return {
        "sdf": cae.to_string("sdf"),
        "fp_p": fp_p,
        "fp_pq": fp_pq,
        "central_atom": central,
        "n_atoms": len(cae.atoms),
        "formula_str": brute_formula,
    }

def process_population_entry(entry_id, with_charges_ids=None):
    try:
        with_charges = True
        reader = io.EntryReader("CSD")
        entry = reader.entry(entry_id)
        mol = component_of_interest(entry.molecule)
        if mol is None:
            return []

        caes = get_caes(mol)
        results = []
        for i, cae in enumerate(caes):
            cae_data = generate_fp_data(cae, with_charges=with_charges)
            cae_data["entry_id"] = entry_id
            cae_data["cae_idx"] = i
            results.append(cae_data)
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
    input_csv = "viable_structures.csv"
    duckdb_path = "all_caes.duckdb"
    batch_size = 5000
    cae_size_counter = defaultdict(int)
    
    # Remove existing DB if you want a fresh start
    if os.path.exists(duckdb_path):
        os.remove(duckdb_path)

    entries = extract_entries_from_csv(input_csv)
    n_proc = min(cpu_count(), 8)
    print(f"⚙️ Using {n_proc} processes")

    con = duckdb.connect(duckdb_path)

    # 💡 Create table ONCE at the start
    con.execute("""
        CREATE TABLE caes (
            entry_id TEXT,
            cae_idx INT,
            sdf TEXT,
            fp_p TEXT,
            fp_pq TEXT,
            central_atom TEXT,
            n_atoms INT,
            formula_str TEXT,
            PRIMARY KEY (entry_id, cae_idx)
        )
    """)

    total_caes_written = 0
    batch_number = 0

    for batch in chunk_list(entries, batch_size):
        batch_number += 1
        print(f"🔄 Processing batch {batch_number} / {len(entries)//batch_size + 1}")
        start_batch = time.time()
        
        results = []
        with Pool(n_proc) as pool:
            for cae_list in pool.imap_unordered(process_population_entry, batch):
                results.extend(cae_list)

        if not results:
            print(f"✅ No CAEs in batch {batch_number}; skipping.")
            continue

        df = pd.DataFrame(results)
        
        # 🔢 Update CAE size distribution
        for size, count in df["n_atoms"].value_counts().items():
            cae_size_counter[size] += count

        # 💡 Ensure all expected columns exist
        required_cols = [
            "entry_id", 
            "cae_idx", 
            "sdf", 
            "fp_p", 
            "fp_pq", 
            "central_atom", 
            "n_atoms", 
            "formula_str"
        ]
        for col in required_cols:
            if col not in df.columns:
                df[col] = None

        df = df[required_cols]

        # 💡 Convert fingerprint lists to strings
        df["fp_p"] = df["fp_p"].apply(lambda x: str(x) if x is not None else None)
        df["fp_pq"] = df["fp_pq"].apply(lambda x: str(x) if x is not None else None)

        if df.empty:
            print(f"✅ DataFrame for batch {batch_number} is empty; skipping insert.")
            continue


        con.register("df", df)
        con.execute("""
            INSERT INTO caes AS c
            SELECT * FROM df
            ON CONFLICT (entry_id, cae_idx) DO UPDATE SET
                fp_p = COALESCE(c.fp_p, EXCLUDED.fp_p),
                fp_pq = COALESCE(c.fp_pq, EXCLUDED.fp_pq),
                sdf = EXCLUDED.sdf,
                central_atom = EXCLUDED.central_atom,
                n_atoms = EXCLUDED.n_atoms,
                formula_str = EXCLUDED.formula_str
        """)


        batch_cae_count = len(df)
        total_caes_written += batch_cae_count
        print(f"✅ Saved batch {batch_number} ({batch_cae_count} CAEs)")
        batch_time = time.time() - start_batch
        print(f"⏱️ Batch {batch_number} completed in {batch_time:.2f} seconds")

    con.execute("""
        CREATE INDEX idx_formula 
        ON caes(central_atom, n_atoms, formula_str)
    """)
    con.close()
    
    # 📊 Analyze size distribution
    print(f"\n📊 Total CAEs processed: {total_caes_written}")
    print("📈 CAE size distribution:")
    for size in sorted(cae_size_counter):
        count = cae_size_counter[size]
        pct = 100 * count / total_caes_written
        print(f"  Size {size}: {count} CAEs ({pct:.2f}%)")
        
    elapsed = time.time() - start
    print(f"\n✅ All data inserted into '{duckdb_path}' table 'caes'.")
    print(f"🔄 Total batches processed: {batch_number}")
    print(f"⏱️ Total time elapsed: {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()
