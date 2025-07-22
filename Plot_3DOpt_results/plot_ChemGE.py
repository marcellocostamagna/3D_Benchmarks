#!/usr/bin/env python3
"""
plotting_ChemGE.py
Plot grouped ChemGE 3Opt scores with consistent styles and layout.

Usage:
    python plotting_ChemGE.py [csv1.csv csv2.csv csv3.csv]
    (defaults to ChemGE CCDC, OBabel, and RDKit)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re
import textwrap

# ---- Set global font to Arial size 10 ----
import matplotlib
matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

# ------------------------------------------------------------
# 1) Define files and labels (or get from CLI)
# ------------------------------------------------------------
file_labels = [
    ('ChemGE_ccdc', 'CCDC'),
    ('ChemGE_obabel', 'OBabel'),
    ('ChemGE_rdkit', 'RDKit'),
]

if len(sys.argv) > 1:
    file_labels = [(Path(f).stem, Path(f).stem.replace("_", " ")) for f in sys.argv[1:]]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # CCDC, OBabel, RDKit

# ------------------------------------------------------------
# 2) Load and align data
# ------------------------------------------------------------
dfs = []
labels = []

for fname, label in file_labels:
    df = pd.read_csv(fname + ".csv")
    if 'task' not in df.columns or '3DOpt_Score' not in df.columns:
        raise ValueError(f"{fname}.csv must contain 'task' and '3DOpt_Score' columns.")
    dfs.append(df)
    labels.append(label)

tasks = dfs[0]['task'].values
n_tasks = len(tasks)
n_methods = len(dfs)

# Build score matrix
all_scores = np.zeros((n_methods, n_tasks))
for i, df in enumerate(dfs):
    all_scores[i] = df['3DOpt_Score'].values

# ------------------------------------------------------------
# 3) Plotting setup
# ------------------------------------------------------------
bar_width = 0.7 / n_methods
x = np.arange(n_tasks)

fig, ax = plt.subplots(figsize=(6.93, 2.75))
fig.subplots_adjust(bottom=0.30)

for i in range(n_methods):
    bars = ax.bar(
        [p + i * bar_width for p in x],
        all_scores[i],
        width=bar_width,
        label=labels[i],
        color=colors[i % len(colors)],
        edgecolor='black',
        linewidth=0.4
    )

# ------------------------------------------------------------
# 4) Axis styling and labels
# ------------------------------------------------------------
def clean_label(t):
    return re.sub(r'^\d+_', '', t)

ax.set_ylabel("3DOpt Score", fontsize=10, fontname='Arial')
ax.set_xlabel("Task", fontsize=10, fontname='Arial')
# ax.set_title("3DOpt Scores for ChemGE", fontsize=10, fontname='Arial', fontweight='bold')

# Center ticks
ax.set_xticks([p + bar_width * (n_methods - 1) / 2 for p in x])
ax.set_xticklabels(
    [textwrap.fill(clean_label(t), 12) for t in tasks],
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=8,
    fontname='Arial'
)

ax.tick_params(axis='both', labelsize=8)

# Format Y axis as 0–1
ax.set_ylim(0, 1)
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: f"{x:.1f}".rstrip('0').rstrip('.') if '.' in f"{x:.1f}" else f"{int(x)}")
)
for label in ax.get_yticklabels():
    label.set_fontsize(10)
    label.set_fontname('Arial')

# Legend
legend = ax.legend(
    loc='lower right',
    fontsize=7,
    title_fontsize=7,
    prop={'family': 'Arial', 'size': 7},
    framealpha=1.0,
    frameon=True,
    edgecolor='black'
)

legend.get_frame().set_linewidth(0.4)

# Hide top and right spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Tight X limits with small buffer
group_width = bar_width * n_methods
buffer = bar_width * 1.2
left = 0 - buffer
right = n_tasks - 1 + group_width + buffer - bar_width
ax.set_xlim(left, right)

# ------------------------------------------------------------
# 5) Save the figure
# ------------------------------------------------------------
fig.tight_layout(pad=0.5)
output_path = Path("3DOpt_ChemGE_Scores.svg")
fig.savefig(output_path, format='svg')
plt.show()
print(f"✅ Plot saved to {output_path.resolve()}")
