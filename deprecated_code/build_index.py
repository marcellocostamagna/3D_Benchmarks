import pandas as pd
import time
import ast  # safer than eval
from collections import Counter

start = time.time()

csv_file = "all_fragments_data_all.csv"
chunk_size = 100_000
index_output = "fragment_index.csv"

total_fragments = 0
fragment_size_counts = Counter()
index_rows = []

reader = pd.read_csv(csv_file, chunksize=chunk_size)
for chunk_id, chunk in enumerate(reader):
    total_fragments += len(chunk)
    fragment_size_counts.update(chunk['n_atoms'].value_counts().to_dict())

    for row_idx, row in chunk.iterrows():
        try:
            raw_formula = ast.literal_eval(row["formula"])  # should be like ('C', 'C2O1')
            central_atom = raw_formula[0]
            formula_str = raw_formula[1]
            n_atoms = int(row["n_atoms"])

            formula_key = (central_atom, n_atoms, formula_str)

            index_rows.append({
                "chunk_id": chunk_id,
                "row_in_chunk": row_idx,
                "entry_id": row["entry_id"],
                "formula": formula_key
            })

        except Exception as e:
            continue  # Skip malformed rows

# Save index
pd.DataFrame(index_rows).to_csv(index_output, index=False)

# Report counts
print(f"Total fragments: {total_fragments}")
print("Fragment size distribution:")
for size, count in sorted(fragment_size_counts.items()):
    percentage = (count / total_fragments) * 100
    print(f"  Size {size}: {count} fragments ({percentage:.2f}%)")

print(f"Index saved to {index_output}")
print(f"Total time elapsed: {time.time() - start:.2f} seconds")
