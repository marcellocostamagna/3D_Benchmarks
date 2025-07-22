import duckdb
import os
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter

# === CONFIGURATION ===
DUCKDB_PATH = "all_fragments_2.duckdb"
ATOM_COUNT = 21
SAVE_EXAMPLE = True
OUTPUT_SDF_PATH = f"example_n{ATOM_COUNT}_fragment.sdf"

# === CONNECT TO DB ===
con = duckdb.connect(DUCKDB_PATH)

# === FETCH FRAGMENTS WITH DESIRED SIZE ===
query = f"""
    SELECT entry_id, frag_idx, sdf
    FROM fragments
    WHERE n_atoms = {ATOM_COUNT}
"""
df = con.execute(query).fetchdf()

# === DISPLAY RESULTS ===
n_total = len(df)
unique_entries = df["entry_id"].nunique()

print(f"🔍 Found {n_total} fragments with {ATOM_COUNT} atoms")
print(f"🧬 Unique CSD entries: {unique_entries}")
# print(f"🆔 Entry IDs: {sorted(df['entry_id'].unique())}")
print(f"🆔 Entry IDs: {sorted(df['entry_id'].unique())[:10]}{'...' if unique_entries > 10 else ''}")


# === SAVE EXAMPLE FRAGMENT ===
if SAVE_EXAMPLE and not df.empty:
    example_sdf = df.iloc[0]["sdf"]
    mol = Molecule.from_string(example_sdf, format="sdf")
    with MoleculeWriter(OUTPUT_SDF_PATH) as writer:
        writer.write(mol)
    print(f"💾 Saved example fragment to: {OUTPUT_SDF_PATH}")

# === FRAGMENTS WITH CENTRAL ATOM = 'La' ===
print("\n📊 Fragment size distribution for central atom 'La':")

la_query = """
    SELECT n_atoms, COUNT(*) AS count
    FROM fragments
    WHERE central_atom = 'La'
    GROUP BY n_atoms
    ORDER BY n_atoms
"""
la_dist = con.execute(la_query).fetchdf()

# === DISPLAY La STATISTICS ===
total_la = la_dist["count"].sum()
print(f"🔬 Total La-centered fragments: {total_la}")
for _, row in la_dist.iterrows():
    size = row["n_atoms"]
    count = row["count"]
    pct = (count / total_la) * 100
    print(f"  - Size {size}: {count} fragments ({pct:.2f}%)")

con.close()



