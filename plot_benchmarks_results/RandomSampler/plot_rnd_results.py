import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Load data from CSV files
file_obabel = "plot_benchmarks_results/RandomSampler/average_results_rnd_obabel.csv"
file_ccdc = "plot_benchmarks_results/RandomSampler/average_results_rnd_ccdc.csv"
file_rdkit = "plot_benchmarks_results/RandomSampler/average_results_rnd_rdkit.csv"


df_obabel = pd.read_csv(file_obabel, sep=",")  # Ensure correct separator
df_ccdc = pd.read_csv(file_ccdc, sep=",")
df_rdkit = pd.read_csv(file_rdkit, sep=",")

# Ensure column names are stripped of spaces
df_obabel.columns = df_obabel.columns.str.strip()
df_ccdc.columns = df_ccdc.columns.str.strip()
df_rdkit.columns = df_rdkit.columns.str.strip()

# Extract tasks and sort by their numeric prefix
def extract_task_number(task_name):
    match = re.match(r"(\d+)_", task_name)  # Extract leading number
    return int(match.group(1)) if match else float('inf')  # Handle non-matching cases

# Sorting tasks based on their numeric prefix
df_ccdc["task_order"] = df_ccdc["task"].apply(extract_task_number)
df_obabel["task_order"] = df_obabel["task"].apply(extract_task_number)
df_rdkit["task_order"] = df_rdkit["task"].apply(extract_task_number)

df_ccdc = df_ccdc.sort_values(by="task_order")
df_obabel = df_obabel.sort_values(by="task_order")
df_rdkit = df_rdkit.sort_values(by="task_order")

# Get sorted task names
task_order = df_ccdc["task"].tolist()

# Create a dictionary of scores for lookup
ccdc_dict = dict(zip(df_ccdc["task"], df_ccdc["3D_Benchmark_Score"]))
obabel_dict = dict(zip(df_obabel["task"], df_obabel["3D_Benchmark_Score"]))
rdkit_dict = dict(zip(df_rdkit["task"], df_rdkit["3D_Benchmark_Score"]))

# Fill in missing values with 0
obabel_scores = [obabel_dict.get(task, 0) for task in task_order]
ccdc_scores = [ccdc_dict.get(task, 0) for task in task_order]
rdkit_scores = [rdkit_dict.get(task, 0) for task in task_order]

# Compute sums and averages
total_ccdc = sum(ccdc_scores)
total_obabel = sum(obabel_scores)
total_rdkit = sum(rdkit_scores)
avg_ccdc = total_ccdc / len(task_order)
avg_obabel = total_obabel / len(task_order)
avg_rdkit = total_rdkit / len(task_order)
summary_text = (f"CCDC: {total_ccdc:.2f}/30 (Avg: {avg_ccdc:.2f})\n"
                f"Obabel: {total_obabel:.2f}/30 (Avg: {avg_obabel:.2f})\n"
                f"RDKit: {total_rdkit:.2f}/30 (Avg: {avg_rdkit:.2f})")

# Define X-axis positions
x = np.arange(len(task_order))
width = 0.3  # Adjust width for three bars

# Plot histogram
fig, ax = plt.subplots(figsize=(12, 6))
# ax.bar(x - width, ccdc_scores, width, label='CCDC', color='dodgerblue', alpha=0.7)  # Left
# ax.bar(x, obabel_scores, width, label='Obabel', color='orangered', alpha=0.7)  # Center
# ax.bar(x + width, rdkit_scores, width, label='RDKit', color='mediumseagreen', alpha=0.7)  # Right

ax.bar(x - width, ccdc_scores, width, label='CCDC', color='maroon', alpha=0.7)  # Left
ax.bar(x, obabel_scores, width, label='Obabel', color='firebrick', alpha=0.7)  # Center
ax.bar(x + width, rdkit_scores, width, label='RDKit', color='indianred', alpha=0.7)  # Right


# Labels and formatting
ax.set_xlabel("Tasks")
ax.set_ylabel("3D Benchmark Score")
ax.set_title("3D Benchmarks Results with Random Sampler")
ax.set_xticks(x)
ax.set_xticklabels(task_order, rotation=45, ha="right", fontsize=8)  # Align labels to the right
ax.set_ylim(0, 1)
ax.legend(loc="upper left")  # Move legend to the top left

# Remove border lines (only show axes)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)

# Add sum and average text to the figure in the top right corner
ax.text(0.95, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

# Show plot
plt.tight_layout()
plt.show()
