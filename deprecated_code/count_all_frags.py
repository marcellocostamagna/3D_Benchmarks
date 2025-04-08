import pandas as pd
import time
from collections import Counter

# Start timer
start = time.time()

csv_file = "all_fragments_data_all.csv"
chunk_size = 100_000  # Adjust based on available RAM

total_fragments = 0
fragment_size_counts = Counter()

# Read CSV in chunks
for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
    total_fragments += len(chunk)
    fragment_size_counts.update(chunk['n_atoms'].value_counts().to_dict())

# Compute and print results
print(f"Total fragments: {total_fragments}")
print("Fragment size distribution:")
for size, count in sorted(fragment_size_counts.items()):
    percentage = (count / total_fragments) * 100
    print(f"  Size {size}: {count} fragments ({percentage:.2f}%)")

# End timer
end = time.time()
print(f"Total time elapsed: {end - start:.2f} seconds")
