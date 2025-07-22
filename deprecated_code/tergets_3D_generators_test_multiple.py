#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-run benchmarking for CSD to 3D structure generation.
Organizes results as run_{i}/ for each run, and Average_results/ for stats.
"""

import os, time, numpy as np, multiprocessing, pandas as pd
from collections import defaultdict
import random

from ccdc import io
from ccdc.conformer import ConformerGenerator
from ccdc.molecule import Molecule
from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt


from hsr import pre_processing as pp
from hsr import fingerprint as fp
from hsr import similarity as sim
from hsr.utils import PROTON_FEATURES

# =========================== Parameters ============================
TARGETS = [
    'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE', 'NIWPUE01',
    'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG', 'ABOBUP', 'XIDTOW',
    'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ', 'ADUPAS', 'DAJLAC', 'OFOWIS',
    'CATSUL', 'HESMUQ01', 'GUDQOL', 'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA',
    'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM'
]
N_RUNS = 100

CCDC_TIMEOUT = 30  # seconds per molecule
BASE_OUT = os.path.abspath("3D_Generators_analysis_last_100")
AVG_OUT = os.path.join(BASE_OUT, "Average_results")
os.makedirs(AVG_OUT, exist_ok=True)

NA = "N.A."

def fmt(v):
    return NA if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))) else f"{v:.2f}"

# =========================== Utility functions ============================

def component_of_interest(molecule):
    comps = molecule.components
    if not comps: return None
    props = [{
        "c": c,
        "is_organomet": c.is_organometallic,
        "mw": sum(a.atomic_weight for a in c.atoms),
        "nat": len(c.atoms)
    } for c in comps]
    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["nat"])
    for p in props:
        score = sum([p["is_organomet"], p["c"] is heaviest["c"], p["c"] is most_atoms["c"]])
        if score >= 2 and p["nat"] >= 5: return p["c"]
    return None

def get_array_from_ccdcmol(ccdcmol):
    arr = np.array([[a.coordinates[0], a.coordinates[1], a.coordinates[2], np.sqrt(a.atomic_number)] for a in ccdcmol.atoms])
    arr -= arr.mean(axis=0)
    return arr

def get_array_from_pybelmol(pybelmol):
    arr = np.array([[a.coords[0], a.coords[1], a.coords[2], np.sqrt(a.atomicnum)] for a in pybelmol.atoms])
    arr -= arr.mean(axis=0)
    return arr

# =========================== Timeout-protected CCDC =======================

def _ccdc_worker(smiles, q, run_seed):
    try:
        random.seed(run_seed)
        cg = ConformerGenerator()
        nconf = random.randint(2, 100)     
        cg.settings.max_conformers = nconf
        mol0 = Molecule.from_string(smiles)
        hits = cg.generate(mol0).hits
        if not hits:
            q.put(None)
        else:
            # idx = random.randint(0, len(hits)-1)
            q.put(hits[0].molecule.to_string("sdf"))
    except Exception:
        q.put(None)

def generate_ccdc_sdf(smiles, timeout=CCDC_TIMEOUT, run_seed=None):
    ctx = multiprocessing.get_context("fork")
    q = ctx.Queue()
    p = ctx.Process(target=_ccdc_worker, args=(smiles, q, run_seed))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.kill()
        p.join()
    return q.get() if not q.empty() else None

def gen_with_ccdc(smiles, run_seed):
    sdf = generate_ccdc_sdf(smiles, timeout=CCDC_TIMEOUT, run_seed=run_seed)
    if sdf is None:
        raise RuntimeError("CCDC generation failed or timed-out")
    mol = Molecule.from_string(sdf)
    return mol, sdf

def gen_with_rdkit(smiles, run_seed):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("SMILES cannot be parsed by RDKit")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDG()
    params.randomSeed = run_seed
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol, Chem.MolToMolBlock(mol)

def gen_with_obabel(smiles, run_seed):
    os.environ["OB_RANDOM_SEED"] = str(run_seed)
    mol = pb.readstring("smi", smiles)
    mol.addh()
    mol.make3D()
    mol.localopt()
    return mol, mol.write("sdf")


# =========================== Main Run Loop ================================

reader = io.EntryReader("CSD")
methods = ["ccdc", "rdkit", "obabel"]

# --- store all runs results for averaging ---
all_run_tables = []

for run in range(1, N_RUNS+1):
    run_seed = np.random.randint(0, 2**32-1)
    print(f"\n\n{'='*15} Starting run {run}/{N_RUNS} {'='*15}")
    OUT_DIR = os.path.join(BASE_OUT, f"run_{run}")
    os.makedirs(OUT_DIR, exist_ok=True)
    sdf_dirs = {m: os.path.join(OUT_DIR, m + "_sdfs") for m in methods}
    for d in sdf_dirs.values(): os.makedirs(d, exist_ok=True)
    stats = {m: [] for m in methods}
    run_rows = []
    for i, refcode in enumerate(TARGETS, 1):
        print(f"\nTarget {i}/{len(TARGETS)}  {refcode}")

        # --- Get reference structure ---
        try:
            mol_cryst = component_of_interest(reader.entry(refcode).molecule)
            smiles    = mol_cryst.smiles
            if smiles is None:
                print(f"Could not obtain SMILES from target {refcode}")
            arr_orig  = get_array_from_ccdcmol(mol_cryst)
            fp_orig   = fp.generate_fingerprint_from_data(arr_orig)
        except Exception as e:
            print(f"  ❌  Could not load crystal structure: {e}")
            row = {"Target": refcode}
            for m in methods:
                row[f"{m}_sim"] = NA
                row[f"{m}_time"] = NA
                row[f"{m}_ok"] = "–"
            run_rows.append(row)
            continue

        # --- Try each generator ---
        res = {}
        for m, func in zip(methods, [gen_with_ccdc, gen_with_rdkit, gen_with_obabel]):
            t0 = time.time()
            try:
                mol3d, sdf = func(smiles, run_seed)
                elapsed = time.time() - t0
                if m == "obabel":
                    arr = get_array_from_pybelmol(mol3d)
                elif m == "rdkit":
                    arr = pp.molecule_to_ndarray(mol3d, features=PROTON_FEATURES, removeHs=True)
                else:  # ccdc
                    arr = get_array_from_ccdcmol(mol3d)
                fp_gen  = fp.generate_fingerprint_from_data(arr)
                sim_val = sim.compute_similarity_score(fp_orig, fp_gen)
                stats[m].append((elapsed, sim_val))
                res[m] = {"time": elapsed, "sim": sim_val, "sdf": sdf}
                # Save SDF in dedicated subfolder for this method
                sdf_filename = os.path.join(sdf_dirs[m], f"{i}_{refcode}_{m}.sdf")
                with open(sdf_filename, "w") as fh:
                    fh.write(sdf)
            except Exception as e:
                print(f"  ⚠️  {m.upper()} failed: {e}")
                res[m] = None

        # --- Report ---
        row = {"Target": refcode}
        for m in methods:
            row[f"{m}_sim"] = fmt(res[m]["sim"]) if res[m] else NA
            row[f"{m}_time"] = fmt(res[m]["time"]) if res[m] else NA
            row[f"{m}_ok"] = "✓" if res[m] else "–"
            print(f"  {m.upper():6} time: {row[f'{m}_time']} s   similarity: {row[f'{m}_sim']}")
        run_rows.append(row)

    # --- Save CSV with all per-target scores and times for this run ---
    df = pd.DataFrame(run_rows)
    run_csv = os.path.join(OUT_DIR, f"scores_and_times_run{run}.csv")
    df.to_csv(run_csv, index=False)
    print(f"\n[Run {run}] Results saved to {run_csv}")
    all_run_tables.append(df)

    # --- Save summary table for this run (✓/–) ---
    summary_table = df[["Target"] + [f"{m}_ok" for m in methods]]
    summary_txt = os.path.join(OUT_DIR, f"success_summary_run{run}.txt")
    summary_table.to_csv(summary_txt, sep="\t", index=False)
    print(f"[Run {run}] Success table saved to {summary_txt}")

# =================== Average Across Runs and Save ======================

# Stack all runs into a single DataFrame
all_scores = defaultdict(lambda: defaultdict(list))  # {method: {target: [scores,...]}}

for df in all_run_tables:
    for _, row in df.iterrows():
        tgt = row["Target"]
        for m in methods:
            try:
                sim = float(row[f"{m}_sim"]) if row[f"{m}_sim"] != NA else np.nan
                t   = float(row[f"{m}_time"]) if row[f"{m}_time"] != NA else np.nan
            except Exception:
                sim = np.nan
                t = np.nan
            all_scores[m][tgt].append({"sim": sim, "time": t})

# Now compute averages, stddevs, N for each method & target
avg_rows = []
for tgt in TARGETS:
    row = {"Target": tgt}
    for m in methods:
        vals = all_scores[m][tgt]
        sims = [v["sim"] for v in vals if not np.isnan(v["sim"])]
        times = [v["time"] for v in vals if not np.isnan(v["time"])]
        row[f"{m}_sim_mean"] = fmt(np.mean(sims)) if sims else NA
        row[f"{m}_sim_std"]  = fmt(np.std(sims))  if sims else NA
        row[f"{m}_time_mean"] = fmt(np.mean(times)) if times else NA
        row[f"{m}_time_std"]  = fmt(np.std(times))  if times else NA
        row[f"{m}_n_succ"] = len(sims)
    avg_rows.append(row)

avg_df = pd.DataFrame(avg_rows)
avg_csv = os.path.join(AVG_OUT, "average_scores_and_times.csv")
avg_df.to_csv(avg_csv, index=False)
print(f"\nAveraged results across {N_RUNS} runs saved to: {avg_csv}")

# =====================  Plotting and Success Rate Reporting =======================

method_markers = {
    "ccdc": "o",     
    "rdkit": "s",     
    "obabel": "D",    
}

# 1. Plot: Average times per target per method
plt.figure(figsize=(14,6))
for m in methods:
    plt.scatter(
        avg_df['Target'], 
        [float(x) if x != NA else np.nan for x in avg_df[f"{m}_time_mean"]],
        label=m.upper(), marker=method_markers[m], s=60, alpha=0.8
    )
plt.xticks(rotation=45, ha='right', fontsize=7)
plt.ylabel('Average time per target (s)')
plt.title('Average 3D Generation Time per Target')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(AVG_OUT, "average_times_per_target.png"))
plt.close()

# 2. Plot: Average similarities per target per method
plt.figure(figsize=(14,6))
for m in methods:
    plt.scatter(
        avg_df['Target'], 
        [float(x) if x != NA else np.nan for x in avg_df[f"{m}_sim_mean"]],
        label=m.upper(), marker=method_markers[m], s=60, alpha=0.8,
    )
plt.xticks(rotation=45, ha='right', fontsize=7)
plt.ylabel('Average similarity per target')
plt.title('Average 3D Similarity per Target')
plt.legend()
plt.ylim(0, 1.0)
plt.tight_layout()
plt.savefig(os.path.join(AVG_OUT, "average_similarities_per_target.png"))
plt.close()

# 3. Report & save success rates per method
success_pct = {}
for m in methods:
    num_success = sum([int(x) for x in avg_df[f"{m}_n_succ"]])
    pct = 100.0 * num_success / (len(TARGETS)*N_RUNS)
    success_pct[m] = pct
    print(f"{m.upper():6}:  {num_success}/{len(TARGETS)*N_RUNS} ({pct:.1f}%) successful generations.")

# Save to txt for reference
with open(os.path.join(AVG_OUT, "success_rates.txt"), "w") as f:
    for m in methods:
        f.write(f"{m.upper():6}:  {success_pct[m]:.1f}% successful generations "
                f"({sum([int(x) for x in avg_df[f'{m}_n_succ']])}/{len(TARGETS)*N_RUNS})\n")

print(f"\nPlots and success summary saved in {AVG_OUT}/\n")


# --- Also print simple per-method summary across ALL runs ---
print("\n=== Method averages over all targets ===")
for m in methods:
    all_sims = [float(r[f"{m}_sim_mean"]) for r in avg_rows if r[f"{m}_sim_mean"] != NA]
    all_times = [float(r[f"{m}_time_mean"]) for r in avg_rows if r[f"{m}_time_mean"] != NA]
    print(f"{m.upper():6}:  N={len(all_sims):2d}    Avg similarity={np.mean(all_sims):.3f}    Avg time={np.mean(all_times):.3f}")

print("\n✅  Finished.")

