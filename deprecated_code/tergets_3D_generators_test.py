#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the crystal structure of CSD entries with 3-D models generated
by CCDC, RDKit and Open Babel.

If any generator fails for a particular target, its time and similarity
are reported as "N.A." and the script continues with the next generator.
"""

from ccdc import io
from ccdc.conformer import ConformerGenerator
from ccdc.molecule import Molecule

from openbabel import pybel as pb

from rdkit import Chem
from rdkit.Chem import AllChem

from hsr import pre_processing as pp
from hsr import fingerprint as fp
from hsr import similarity as sim
from hsr.utils import PROTON_FEATURES

import numpy as np
import time, os, sys, multiprocessing   

##############################################################################
# parameters & helpers
##############################################################################

TARGETS = [
    'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE', 'NIWPUE01',
    'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG', 'ABOBUP', 'XIDTOW',
    'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ', 'ADUPAS', 'DAJLAC', 'OFOWIS',
    'CATSUL', 'HESMUQ01', 'GUDQOL', 'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA',
    'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM'
]

OUT_DIR = "output"
os.makedirs(OUT_DIR, exist_ok=True)

NA = "N.A."                    


def fmt(v):
    """Pretty-format a number or return NA on None/str."""
    if v is None or isinstance(v, str):
        return NA
    return f"{v:.2f}"

##############################################################################
# timeout-protected CCDC wrapper (NEW)
##############################################################################

def _ccdc_worker(smiles, q):
    """Child process: build one conformer and put SDF (or None) in queue."""
    try:
        cg   = ConformerGenerator()
        mol0 = Molecule.from_string(smiles)
        hit  = cg.generate(mol0).hits[0]
        q.put(hit.molecule.to_string("sdf"))
    except Exception:
        q.put(None)

def generate_ccdc_sdf(smiles: str, timeout: int = 50):
    """
    Run CCDC conformer generation in a separate process.
    Returns an SDF string, or None if it errors or exceeds *timeout* seconds.
    """
    ctx = multiprocessing.get_context("fork")   # use fork on Unix; change to "spawn" on Windows
    q   = ctx.Queue()
    p   = ctx.Process(target=_ccdc_worker, args=(smiles, q))
    p.start()
    p.join(timeout)

    if p.is_alive():        # timeout hit
        p.kill()
        p.join()

    return q.get() if not q.empty() else None


##############################################################################
# utility functions
##############################################################################

def component_of_interest(molecule):
    comps = molecule.components
    if not comps:
        return None

    props = [{
        "c": c,
        "is_organomet": c.is_organometallic,
        "mw": sum(a.atomic_weight for a in c.atoms),
        "nat": len(c.atoms)
    } for c in comps]

    heaviest   = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["nat"])

    for p in props:
        score = sum([
            p["is_organomet"],
            p["c"] is heaviest["c"],
            p["c"] is most_atoms["c"]
        ])
        if score >= 2 and p["nat"] >= 5:
            return p["c"]
    return None


def get_array_from_ccdcmol(ccdcmol):
    arr = np.array([[a.coordinates[0],
                     a.coordinates[1],
                     a.coordinates[2],
                     np.sqrt(a.atomic_number)] for a in ccdcmol.atoms])
    arr -= arr.mean(axis=0)
    return arr


def get_array_from_pybelmol(pybelmol):
    arr = np.array([[a.coords[0],
                     a.coords[1],
                     a.coords[2],
                     np.sqrt(a.atomicnum)] for a in pybelmol.atoms])
    arr -= arr.mean(axis=0)
    return arr


##############################################################################
# 3-D generators (each one raises on failure)
##############################################################################

def gen_with_ccdc(smiles):
    """Return (ccdc_molecule, sdf) or raise if timeout/failure."""
    sdf = generate_ccdc_sdf(smiles, timeout=30)
    if sdf is None:
        raise RuntimeError("CCDC generation failed or timed-out")
    mol = Molecule.from_string(sdf)   # parse back to Molecule for array calc
    return mol, sdf


def gen_with_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES cannot be parsed by RDKit")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    return mol, Chem.MolToMolBlock(mol)


def gen_with_obabel(smiles):
    mol = pb.readstring("smi", smiles)
    mol.addh()
    mol.make3D()
    mol.localopt()
    return mol, mol.write("sdf")


##############################################################################
# driver
##############################################################################

reader = io.EntryReader("CSD")
summary_rows = []                                
stats = {"ccdc": [], "rdkit": [], "obabel": []}

for i, refcode in enumerate(TARGETS, 1):
    print(f"\nTarget {i}/{len(TARGETS)}  {refcode}")

    # ------------------------------------------------------------------ crystal
    try:
        mol_cryst = component_of_interest(reader.entry(refcode).molecule)
        smiles    = mol_cryst.smiles
        if smiles is None:
            print(f"Could not obtain SMILES from target {refcode}")
        arr_orig  = get_array_from_ccdcmol(mol_cryst)
        fp_orig   = fp.generate_fingerprint_from_data(arr_orig)
    except Exception as e:
        print(f"  ❌  Could not load crystal structure: {e}")
        continue        # nothing to compare, go to next refcode

    # ------------------------------------------------------------------ helpers
    def attempt(label, func):
        """Run a generator + downstream analysis, return dict or None."""
        t0 = time.time()
        try:
            mol3d, sdf = func(smiles)
            elapsed    = time.time() - t0

            # downstream
            if label == "obabel":
                arr = get_array_from_pybelmol(mol3d)
            elif label == "rdkit":
                arr = pp.molecule_to_ndarray(mol3d, features=PROTON_FEATURES,
                                             removeHs=True)
            else:                # ccdc
                arr = get_array_from_ccdcmol(mol3d)

            fp_gen  = fp.generate_fingerprint_from_data(arr)
            sim_val = sim.compute_similarity_score(fp_orig, fp_gen)
            stats[label].append((elapsed, sim_val))

            return {"time": elapsed, "sim": sim_val, "sdf": sdf}
        except Exception as e:
            print(f"  ⚠️  {label.upper()} failed: {e}")
            return None

    # ------------------------------------------------------------------ runs
    res_ccdc  = attempt("ccdc",  gen_with_ccdc)
    res_rdkit = attempt("rdkit", gen_with_rdkit)
    res_obab  = attempt("obabel", gen_with_obabel)

    # ------------------------------------------------------------------ report
    print(f"  CCDC  time: {fmt(res_ccdc  and res_ccdc['time'])} s   "
          f"similarity: {fmt(res_ccdc  and res_ccdc['sim'])}")
    print(f"  RDKit time: {fmt(res_rdkit and res_rdkit['time'])} s   "
          f"similarity: {fmt(res_rdkit and res_rdkit['sim'])}")
    print(f"  OBabel time: {fmt(res_obab  and res_obab['time'])} s   "
          f"similarity: {fmt(res_obab  and res_obab['sim'])}")

    # ------------------------------------------------------------------ save SDFs
    def maybe_write(label, res):
        if res is None:                 # generation failed
            return
        path = os.path.join(OUT_DIR, f"{i}_{refcode}_{label}.sdf")
        with open(path, "w") as fh:
            fh.write(res["sdf"])

    maybe_write("ccdc",  res_ccdc)
    maybe_write("rdkit", res_rdkit)
    maybe_write("obabel", res_obab)
    
    summary_rows.append({                             
    "Target": refcode,    
    "CCDC":   "✓" if res_ccdc  else "–",
    "RDKit":  "✓" if res_rdkit else "–",
    "OBabel": "✓" if res_obab  else "–"
})


##############################################################################
# summary table
##############################################################################

print("\n=== Generation summary ===")

try:
    import pandas as pd                 
    print(pd.DataFrame(summary_rows).to_string(index=False))
except ImportError:                 
    head=f"{'Target':<10} {'CCDC':^5} {'RDKit':^6} {'OBabel':^6}"
    print(head); print("-"*len(head))
    for r in summary_rows:
        print(f"{r['Target']:<10} {r['CCDC']:^5} {r['RDKit']:^6} {r['OBabel']:^6}")


print("\n✅  Finished.")
