import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import textwrap
import re
import matplotlib

matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

csv_path = Path("Starting_populations.csv")
output_fig = Path("Starting_populations")

df = pd.read_csv(csv_path)
df.set_index("target", inplace=True)

numeric_cols = [c for c in df.columns if df[c].dtype.kind in "if"]
df = df[numeric_cols]

scale_factor = 1e5
df_scaled = df / scale_factor

num_targets = len(df_scaled)
num_cols = len(df_scaled.columns)
bar_width = 0.8 / num_cols         
x = range(num_targets)              

fig, ax = plt.subplots(figsize=(6.93, 3.3))
fig.subplots_adjust(bottom=0.30)

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for i, col in enumerate(df_scaled.columns):
    bars = ax.bar(
        [p + i * bar_width for p in x],
        df_scaled[col],
        width=bar_width,
        label=str(col),
        color=colors[i % len(colors)],
        edgecolor='black',
        linewidth=0.4
    )

ax.set_ylabel(r"Number of molecules ($\times 10^5$)", fontsize=10, fontname='Arial')
ax.set_xlabel("Task", fontsize=10, fontname='Arial')

def clean_label(t):
    return re.sub(r'^\d+_', '', t)

ax.set_xticks([p + bar_width * (num_cols - 1) / 2 for p in x])
ax.set_xticklabels(
    [textwrap.fill(clean_label(t), 12) for t in df_scaled.index],
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=8,
    fontname='Arial'
)
ax.tick_params(axis='both', labelsize=8)

ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: f"{x:.1f}".rstrip('0').rstrip('.') if '.' in f"{x:.1f}" else f"{int(x)}")
)
for label in ax.get_yticklabels():
    label.set_fontsize(10)
    label.set_fontname('Arial')

ax.legend(
    title="Similarity threshold",
    loc='upper right', 
    fontsize=7,
    title_fontsize=7,
    prop={'family': 'Arial', 'size': 7},
)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

group_width = bar_width * num_cols
buffer = bar_width * 1.2 
left = 0 - buffer
right = num_targets - 1 + group_width + buffer - bar_width
ax.set_xlim(left, right)

fig.tight_layout(pad=0.5)
svg_path = output_fig.with_suffix('.svg')
fig.savefig(svg_path, format='svg')
print(f"✅  Plot saved to {svg_path.resolve()}")
