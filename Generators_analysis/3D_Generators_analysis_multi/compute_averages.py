import pandas as pd

# Load your CSV
df = pd.read_csv("average_scores_and_times.csv")

# Compute average similarities for each method
ccdc_avg = df['ccdc_sim_mean'].mean()
rdkit_avg = df['rdkit_sim_mean'].mean()
obabel_avg = df['obabel_sim_mean'].mean()

# Compute average similarities for each method
ccdc_avg_t = df['ccdc_time_mean'].mean()
rdkit_avg_t = df['rdkit_time_mean'].mean()
obabel_avg_t = df['obabel_time_mean'].mean()

print(f"Average CCDC similarity: {ccdc_avg:.3f}")
print(f"Average RDKit similarity: {rdkit_avg:.3f}")
print(f"Average OBabel similarity: {obabel_avg:.3f}")
print('\n')
print(f"Average CCDC time: {ccdc_avg_t:.3f}")
print(f"Average RDKit time: {rdkit_avg_t:.3f}")
print(f"Average OBabel time: {obabel_avg_t:.3f}")