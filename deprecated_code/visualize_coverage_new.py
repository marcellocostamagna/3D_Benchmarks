#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# list your input CSVs here
csv_files = [
    "fragment_summary_0_3.csv",
    "fragment_summary_0_4.csv",
    "fragment_summary_0_5.csv",
]

for csv_path in csv_files:
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        print(f"⚠️  Skipping {csv_path.name}: file not found")
        continue

    # extract threshold suffix, e.g. "0_3"
    thr = csv_path.stem.split("_")[-1]

    # 1. Load
    df = pd.read_csv(csv_path)

    # 2. Compute totals and percentages
    df['total']       = df['matched'] + df['fallback'] + df['no_match']
    df['pct_perfect'] = df['matched'] / df['total'] * 100
    df['pct_coverage'] = (df['matched'] + df['fallback']) / df['total'] * 100
    df['pct_no_match'] = df['no_match'] / df['total'] * 100

    # overall metrics
    overall_perfect  = df['matched'].sum() / df['total'].sum() * 100
    overall_coverage = (df['matched'] + df['fallback']).sum() / df['total'].sum() * 100

    # prepare x axis
    targets = df['target']
    n       = len(df)
    ind     = np.arange(n)

    # ────────────────────────────────────────────────────────────────────────
    # A) Stacked‐bar breakdown
    # ────────────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(8, n*0.3), 6))
    ax.bar(ind, df['pct_perfect'], 0.8,
           label='Perfect coverage', color='green')
    ax.bar(ind, df['pct_coverage'] - df['pct_perfect'], 0.8,
           bottom=df['pct_perfect'],  label='Fallback part', color='gold')
    ax.bar(ind, df['pct_no_match'], 0.8,
           bottom=df['pct_coverage'],  label='No match',      color='red')

    ax.set_ylabel('Percentage of Fragments')
    ax.set_xlabel('Target')
    ax.set_title(f'Fragment Matching Breakdown (threshold {thr.replace("_",".")})')
    ax.set_xticks(ind)
    ax.set_xticklabels(
        targets,
        rotation=45,
        ha='right',
        rotation_mode='anchor'
    )
    ax.set_ylim(0, 100)

    # opaque legend
    leg1 = ax.legend(
        loc='lower right',
        bbox_to_anchor=(0.04, 0.00),
        framealpha=1)
    leg1.get_frame().set_facecolor('white')

    out1 = f"fragment_summary_{thr}_breakdown.png"
    plt.tight_layout()
    plt.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"✅ Wrote {out1}")

    # ────────────────────────────────────────────────────────────────────────
    # B) Grouped‐bar “histogram” of perfect vs. coverage
    # ────────────────────────────────────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(max(8, n*0.3), 6))
    bar_w = 0.35

    ax2.bar(ind - bar_w/2, df['pct_perfect'],  bar_w,
            label='Perfect coverage', align='center')
    ax2.bar(ind + bar_w/2, df['pct_coverage'], bar_w,
            label='Coverage (incl. fallbacks)', align='center')

    ax2.set_ylabel('Coverage (%)')
    ax2.set_xlabel('Target')
    ax2.set_title(f'Perfect vs. Total Coverage (threshold {thr.replace("_",".")})')
    ax2.set_xticks(ind)
    ax2.set_xticklabels(
        targets,
        rotation=45,
        ha='right',
        rotation_mode='anchor'
    )
    ax2.set_ylim(0, 100)

    # opaque legend
    leg2 = ax2.legend(loc='lower right', framealpha=1)
    leg2.get_frame().set_facecolor('white')

    # annotate overall metrics inside the plot
    textstr = (
        f"Overall perfect: {overall_perfect:.1f}%\n"
        f"Overall total  : {overall_coverage:.1f}%"
    )
    ax2.text(
        0.75, 1.0, textstr,
        transform=ax2.transAxes,
        fontsize=10,
        va='center', ha='left',
        bbox=dict(boxstyle="round,pad=0.3",
                  facecolor="white", alpha=1)
    )

    out2 = f"fragment_summary_{thr}_grouped_coverage.png"
    plt.tight_layout()
    plt.savefig(out2, dpi=150)
    plt.close(fig2)
    print(f"✅ Wrote {out2}\n")
