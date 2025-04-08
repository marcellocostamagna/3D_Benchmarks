import pandas as pd
import time

# Start timer
start = time.time()

csv_file = "all_fragments_data_0_5.csv"
chunk_size = 100_000  # Tune based on your RAM

total_fragments = 0
num_biatomic = 0

# Read CSV in chunks
for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
    total_fragments += len(chunk)
    num_biatomic += (chunk['n_atoms'] == 2).sum()

# Compute percentage
percentage = (num_biatomic / total_fragments) * 100

# Output results
print(f"Total fragments: {total_fragments}")
print(f"Biatomic fragments: {num_biatomic}")
print(f"Percentage: {percentage:.2f}%")

# End timer
end = time.time()
print(f"Total time elapsed: {end - start:.2f} seconds")
