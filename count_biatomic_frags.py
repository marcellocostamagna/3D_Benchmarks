import pandas as pd
import time

start = time.time()
# Load the CSV
csv_file = "all_fragments_data_0_5.csv"
df = pd.read_csv(csv_file)

# Count total number of fragments
total_fragments = len(df)

# Filter biatomic fragments (n_atoms == 2)
biatomic_fragments = df[df['n_atoms'] == 2]
num_biatomic = len(biatomic_fragments)

# Compute percentage
percentage = (num_biatomic / total_fragments) * 100

# Output results
print(f"Total fragments: {total_fragments}")
print(f"Biatomic fragments: {num_biatomic}")
print(f"Percentage: {percentage:.2f}%")
print(f'Total time elapsed: {time.time()-start} seconds')