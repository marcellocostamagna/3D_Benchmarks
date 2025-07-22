#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Load the summary CSV
df = pd.read_csv("fragment_summary_0_5.csv")

# 2. Compute totals and percentages
df['total'] = df['matched'] + df['fallback'] + df['no_match']
df['pct_matched'] = df['matched'] / df['total'] * 100
df['pct_fallback'] = df['fallback'] / df['total'] * 100
df['pct_no_match'] = df['no_match'] / df['total'] * 100

# 3. Prepare stacked bar positions
targets = df['target']
ind = np.arange(len(df))    # the x locations for the groups
width = 0.8                 # the width of the bars

# 4. Plot
fig, ax = plt.subplots(figsize=(max(8, len(df)*0.3), 6))
p1 = ax.bar(ind, df['pct_matched'], width, label='Matched', color='green')
p2 = ax.bar(ind, df['pct_fallback'], width, bottom=df['pct_matched'], label='Fallback', color='gold')
p3 = ax.bar(ind, df['pct_no_match'], width, bottom=df['pct_matched'] + df['pct_fallback'], label='No Match', color='red')

# 5. Labels and aesthetics
ax.set_ylabel('Percentage of Fragments')
ax.set_xlabel('Target')
ax.set_title('Fragment Matching Summary by Target')
ax.set_xticks(ind)
ax.set_xticklabels(targets, rotation=90, ha='center')
ax.set_ylim(0, 100)
ax.legend(loc='upper right')

plt.tight_layout()

# 6. Save to PNG (do NOT show)
plt.savefig("fragment_summary_0_5.png", dpi=150)
