#!/usr/bin/env python3

import time
import ast
import duckdb
import pandas as pd
from collections import Counter

def build_duckdb_with_all_data(
    csv_file="all_fragments_data_all.csv",
    db_path="all_fragments.duckdb",
    chunk_size=100_000
):
    """
    Reads `all_fragments_data_all.csv` in chunks, parses the 'formula' column,
    and writes *all* fragment data into a single DuckDB table named `fragments`,
    including:
      - chunk_id
      - row_in_chunk
      - entry_id
      - (parsed) central_atom, n_atoms, formula_str
      - plus original CSV columns: sdf, fp, etc.

    Finally creates an index on (central_atom, n_atoms, formula_str, entry_id)
    so we can quickly filter for a formula and a set of entry_ids.
    """

    start_time = time.time()
    total_fragments = 0
    fragment_size_counts = Counter()

    # 1) Connect to or create the DuckDB database.
    con = duckdb.connect(db_path)

    # 2) Create a table to store all the columns we want.
    #    We'll store chunk_id, row_in_chunk, plus the original columns from the CSV.
    #    (Adjust column definitions as needed. We'll store them as TEXT or INT, etc.)
    #    If it doesn't exist yet, create it:
    con.execute("""
        CREATE TABLE IF NOT EXISTS fragments (
            chunk_id INT,
            row_in_chunk INT,
            entry_id TEXT,
            frag_idx INT,
            sdf TEXT,
            fp TEXT,
            -- the "raw" formula column from CSV, if you want to keep it
            formula TEXT,

            -- these are the derived columns from parsing the formula
            central_atom TEXT,
            n_atoms INT,
            formula_str TEXT
        )
    """)

    # We'll read the CSV in chunks, parse each row, and append to the DB table.
    reader = pd.read_csv(csv_file, chunksize=chunk_size)
    chunk_count = 0

    for chunk_id, chunk in enumerate(reader):
        # Update counters
        total_fragments += len(chunk)
        fragment_size_counts.update(chunk["n_atoms"].value_counts().to_dict())

        # We'll build a DataFrame of rows to insert. We combine columns from the CSV
        # with our derived columns (central_atom, n_atoms, formula_str).
        # Note: chunk["n_atoms"] is already numeric; chunk["entry_id"] is text, etc.
        # We'll store the row_in_chunk, which is basically the local index in `chunk`.
        index_data = {
            "chunk_id": [],
            "row_in_chunk": [],
            "entry_id": [],
            "frag_idx": [],
            "sdf": [],
            "fp": [],
            "formula": [],        # store the raw formula string from the CSV
            "central_atom": [],
            "n_atoms": [],
            "formula_str": [],
        }

        for row_idx, row in chunk.iterrows():
            try:
                raw = ast.literal_eval(row["formula"])  
                # raw is something like ("C", "C4H1")
                central_atom = str(raw[0])        # "C"
                formula_str = str(raw[1])         # "C4H1"
                
                # The integer number of atoms is from the separate column
                n_atoms = int(row["n_atoms"])

                # Now store everything
                index_data["chunk_id"].append(chunk_id)
                index_data["row_in_chunk"].append(row_idx)
                index_data["entry_id"].append(str(row["entry_id"]))
                index_data["frag_idx"].append(int(row["frag_idx"]))
                index_data["sdf"].append(str(row["sdf"]))
                index_data["fp"].append(str(row["fp"]))
                index_data["formula"].append(str(row["formula"]))  # the raw string ("('C','C4H1')")

                # Our derived columns
                index_data["central_atom"].append(central_atom)
                index_data["n_atoms"].append(n_atoms)
                index_data["formula_str"].append(formula_str)

            except Exception as e:
                # Print or debug if you want
                continue


        df_to_insert = pd.DataFrame(index_data)
        if not df_to_insert.empty:
            # Register the DataFrame as a DuckDB table
            con.register("df_to_insert", df_to_insert)
            con.execute("CREATE TEMPORARY TABLE tmp_chunk AS SELECT * FROM df_to_insert")
            con.execute("INSERT INTO fragments SELECT * FROM tmp_chunk")
            con.execute("DROP TABLE tmp_chunk")
            con.unregister("df_to_insert")  # optional but clean

        chunk_count += 1
        print(f"Processed chunk {chunk_id} with {len(df_to_insert)} inserted rows.")

    # 3) Create an index so we can quickly filter by (central_atom, n_atoms, formula_str, entry_id).
    #    Adjust columns to match your typical queries:
    con.execute("""
        CREATE INDEX IF NOT EXISTS idx_formula_entry 
            ON fragments(central_atom, n_atoms, formula_str, entry_id)
    """)

    # Summaries
    print(f"\nTotal fragments processed: {total_fragments}")
    print("Fragment size distribution:")
    for size, count in sorted(fragment_size_counts.items()):
        percentage = (count / total_fragments) * 100
        print(f"  Size {size}: {count} fragments ({percentage:.2f}%)")

    con.close()
    elapsed = time.time() - start_time
    print(f"\nAll data inserted into '{db_path}' table 'fragments'.")
    print(f"Total chunks read: {chunk_count}")
    print(f"Total time elapsed: {elapsed:.2f} seconds")


if __name__ == "__main__":
    build_duckdb_with_all_data()
