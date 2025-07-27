#!/usr/bin/env python3

"""
aggregate_3DOpt_results.py
Recursively collects and averages results.csv from 3DOpt runs (RandomSampler or ChemGE),
and writes per-method CSVs (one row per task, averaged across runs) for easy plotting.
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np

# -- List of all results folders (edit as needed) --
RESULT_FOLDERS = [
    'Results_RndSam_CCDC_Entries',
    'Results_RndSam_CCDC_Smiles',
    'Results_RndSam_OBabel',
    'Results_RndSam_RDKit',
    'Results_ChemGE_RDKit',
    'Results_ChemGE_OBabel',
    'Results_ChemGE_CCDC'
]

# def find_results_csvs(folder):
#     """Recursively yield all results.csv files up to two levels deep."""
#     folder = Path(folder)
#     for root, dirs, files in os.walk(folder):
#         rootpath = Path(root)
#         if rootpath == folder:
#             # only go 2 levels deep
#             for d in dirs:
#                 for subd in (rootpath/d).iterdir():
#                     if subd.is_dir():
#                         for f in subd.glob("results.csv"):
#                             yield f
#         else:
#             for f in rootpath.glob("results.csv"):
#                 yield f
                
def find_results_csvs(folder):
    """Yield all results.csv files up to two levels deep (not recursive)."""
    folder = Path(folder)
    # Level 1: immediate subdirectories
    for d1 in folder.iterdir():
        if d1.is_dir():
            # Level 2: subdirectories of subdirectories
            for d2 in d1.iterdir():
                if d2.is_dir():
                    f = d2 / "results.csv"
                    if f.exists():
                        yield f

def aggregate_method(folder, outname):
    """Aggregate results.csv from all runs for one method folder, average by task."""
    csv_files = list(find_results_csvs(folder))
    if not csv_files:
        print(f"⚠️  No results.csv found in {folder}")
        return None
    # Each CSV: at least columns 'task', '3DOpt_Score'
    dfs = []
    for f in csv_files:
        df = pd.read_csv(f)
        if "task" not in df.columns or "3DOpt_Score" not in df.columns:
            print(f"⚠️  {f} missing expected columns.")
            continue
        dfs.append(df[["task", "3DOpt_Score"]])
    if not dfs:
        print(f"❌ No valid results in {folder}")
        return None
    bigdf = pd.concat(dfs, ignore_index=True)
    agg = bigdf.groupby("task")["3DOpt_Score"].mean().reset_index()
    agg = agg.sort_values("task")
    agg.to_csv(outname, index=False)
    print(f"✅ Saved average results for {folder} to {outname}")

if __name__ == "__main__":
    for folder in RESULT_FOLDERS:
        label = (
            "Rnd_ccdc_entries" if folder == "Results_RndSam_CCDC_Entries" else
            "Rnd_ccdc" if folder == "Results_RndSam_CCDC_Smiles" else
            "Rnd_obabel" if folder == "Results_RndSam_OBabel" else
            "Rnd_rdkit" if folder == "Results_RndSam_RDKit" else
            "ChemGE_ccdc" if folder == "Results_ChemGE_CCDC" else
            "ChemGE_obabel" if folder == "Results_ChemGE_OBabel" else
            "ChemGE_rdkit" if folder == "Results_ChemGE_RDKit" else
            folder.replace("Results_", "")
        )
        aggregate_method(folder, f"{label}.csv")
