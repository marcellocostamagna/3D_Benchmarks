import duckdb

# Paths to your databases
old_db = "all_fragments_3.duckdb"
new_db = "all_caes.duckdb"

# Connect to both (DuckDB can only have one connection at a time, so use two conns in sequence)
con_old = duckdb.connect(old_db)
con_new = duckdb.connect(new_db)

# Create the new table in the new DB
con_new.execute("""
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

# Read all data from the old table and insert into the new table with column rename
# (DuckDB supports reading from another DB via 'read_csv_auto', but here we do it in Python)

df = con_old.execute("SELECT * FROM fragments").fetchdf()
df = df.rename(columns={"frag_idx": "cae_idx"})

# Insert into new table
con_new.register("df", df)
con_new.execute("""
    INSERT INTO caes
    SELECT * FROM df
""")

con_old.close()
con_new.close()
print("✅ Migration complete: all_caes.duckdb created with 'caes' table.")
