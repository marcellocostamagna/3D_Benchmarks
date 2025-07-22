"""
This script compares target Connected Atom Environments (CAEs) against a population database of CAEs
('all_caes.duckdb'), identifies best matches, and outputs similarity data and SDF files.
"""

import os
import time
import sys
import re
import numpy as np
import multiprocessing as mp
import logging
from collections import defaultdict, Counter

from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc import io
from ccdc.entry import Entry

from hsr import fingerprint as fp
from hsr import similarity as sim

import duckdb
from tqdm import tqdm

os.environ["QT_QPA_PLATFORM"] = "xcb"  # For headless systems

# --------- Helper to redirect print to logging ---------
class StreamToLogger:
    def __init__(self, logger, log_level):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass

# --------- Get log file name based on folder ---------
def extract_threshold_from_folder(folder):
    match = re.search(r"Starting_populations_(\d+)_(\d+)", folder)
    if match:
        return f"{match.group(1)}_{match.group(2)}"
    match_simple = re.search(r"Starting_populations_(\d+)", folder)
    if match_simple:
        return match_simple.group(1)
    return "unknown"

# --------- Setup Logging ---------
population_folder = "../Starting_populations/Starting_populations_0_5"  # <--- CHANGE THIS IF NEEDED
threshold_str = extract_threshold_from_folder(population_folder)
log_filename = f"caes_analysis_{threshold_str}.log"
logging.basicConfig(
    filename=log_filename,
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logging.captureWarnings(True)
sys.stdout = StreamToLogger(logging.getLogger("STDOUT"), logging.INFO)
sys.stderr = StreamToLogger(logging.getLogger("STDERR"), logging.ERROR)

### ---------- Utility Functions ---------- ###

def get_p_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)
    return atom_array

def get_p_q_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number), atom.formal_charge])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)
    return atom_array

def fingerprint_key(fp_obj):
    return tuple(fp_obj.tolist() if isinstance(fp_obj, np.ndarray) else fp_obj)

def formula_signature(cae):
    atoms = cae.atoms
    central = atoms[0].atomic_symbol
    n_atoms = len(atoms)
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, n_atoms, formula)

def generate_fp_data(cae, with_charges=False):
    arr = get_p_q_array_from_ccdcmol(cae) if with_charges else get_p_array_from_ccdcmol(cae)
    fp_data = fp.generate_fingerprint_from_data(arr)
    if isinstance(fp_data, np.ndarray):
        fp_data = fp_data.tolist()

    central, n_atoms, brute_formula = formula_signature(cae)
    return {
        "sdf": cae.to_string("sdf"),
        "fp": fp_data, 
        "formula": (central, n_atoms, brute_formula),
        "central_atom": central,
        "n_atoms": n_atoms,
        "formula_str": brute_formula,
    }

def interatomic_distance(sdf_string):
    mol = Molecule.from_string(sdf_string, format="sdf")
    coords = [atom.coordinates for atom in mol.atoms]
    return np.linalg.norm(np.array(coords[0]) - np.array(coords[1]))

def component_of_interest(molecule):
    components = molecule.components
    if not components:
        return None
    props = [{
        "component": c,
        "is_organometallic": c.is_organometallic,
        "mw": sum(atom.atomic_weight for atom in c.atoms),
        "atom_count": len(c.atoms)
    } for c in components]
    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["atom_count"])
    for prop in props:
        if sum([
            prop["is_organometallic"],
            prop["component"] == heaviest["component"],
            prop["component"] == most_atoms["component"]
        ]) >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None

def create_cae(central_atom):
    cae = Molecule(identifier=f"{central_atom.label}_cae")
    atom_map = {central_atom: cae.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = cae.add_atom(neighbor)
    for bond in central_atom.bonds:
        a1, a2 = bond.atoms
        if a1 in atom_map and a2 in atom_map:
            try:
                cae.add_bond(bond.bond_type, atom_map[a1], atom_map[a2])
            except:
                pass
    return cae

def get_caes(mol):
    return [create_cae(atom) for atom in mol.atoms]

def process_target(entry_id, threshold=0.999, with_charges=False):
    reader = io.EntryReader("CSD")
    mol = component_of_interest(reader.entry(entry_id).molecule)
    caes = get_caes(mol)
    grouped = defaultdict(list)
    for cae in caes:
        data = generate_fp_data(cae, with_charges=with_charges)
        grouped[data["formula"]].append(data)
    unique_by_formula = defaultdict(list)
    for formula, caes_group in grouped.items():
        for cae in caes_group:
            if all(sim.compute_similarity_score(cae["fp"], other["fp"]) < threshold
                   for other in unique_by_formula[formula]):
                unique_by_formula[formula].append(cae)
    return unique_by_formula

def _parse_fp_chunk(rows, with_sdf: bool):
    parsed = []
    if with_sdf:
        for cae_id, fp_str, sdf in rows:
            parsed.append(
                {"cae_id": cae_id,
                 "fp"     : _safe_fp_list(fp_str),
                 "sdf"    : sdf}
            )
    else:
        for cae_id, fp_str in rows:
            parsed.append(
                {"cae_id": cae_id,
                 "fp"     : _safe_fp_list(fp_str)}
            )
    return parsed

def _safe_fp_list(fp_str):
    fp_str = str(fp_str).strip()
    if fp_str.startswith("[") and fp_str.endswith("]"):
        try:
            return [float(x) for x in fp_str.strip("[] ").split(",")]
        except ValueError:
            pass
    return []

def fetch_sdf_for_cae_ids(db_path, cae_ids):
    if not cae_ids:
        return {}
    con = duckdb.connect(db_path)
    con.execute("CREATE TEMP TABLE tmp_ids(cae_id INT)")
    con.executemany("INSERT INTO tmp_ids VALUES (?)", [(i,) for i in cae_ids])
    query = """
        SELECT caes.rowid AS cae_id, caes.sdf
        FROM caes
        JOIN tmp_ids ON caes.rowid = tmp_ids.cae_id
    """
    df = con.execute(query).fetchdf()
    con.close()
    return {row["cae_id"]: row["sdf"] for _, row in df.iterrows()}

def compare_formula_streaming(
        target_caes, db_path, pop_ids, formula,
        threshold=0.999, n_processes=4, chunk_limit=600_000,
        with_charges=False):

    (central_atom, n_atoms, formula_str) = formula
    con = duckdb.connect(db_path)
    con.execute("CREATE TEMP TABLE tmp_pop_ids(entry_id TEXT)")
    con.executemany("INSERT INTO tmp_pop_ids VALUES (?)", [(eid,) for eid in pop_ids])
    need_sdf = target_caes[0]["n_atoms"] == 2
    fp_column = "fp_pq" if with_charges else "fp_p"
    select_cols = f"c.rowid, c.{fp_column}" + (", c.sdf" if need_sdf else "")

    base_query = f"""
        SELECT {select_cols}
        FROM caes c
        JOIN tmp_pop_ids p ON c.entry_id = p.entry_id
        WHERE c.central_atom = '{central_atom}'
          AND c.n_atoms     = {n_atoms}
          AND c.formula_str = '{formula_str}'
    """

    offset = 0
    target_keys = {fingerprint_key(t["fp"]) for t in target_caes}
    unmatched_keys = set(target_keys)
    best_matches = defaultdict(list)
    distorted_match_scores = {}

    while True:
        chunk_query = base_query + f" LIMIT {chunk_limit} OFFSET {offset}"
        t0 = time.time()
        rows = con.execute(chunk_query).fetchall()

        if offset == 0 and not rows:
            print(f"⚠️ No CAEs found in DB for formula {formula}")
            for tdata in target_caes:
                tkey = fingerprint_key(tdata["fp"])
                distorted_match_scores[tkey] = (-1, "NO_MATCH_IN_DB")
            break

        if not rows:
            break

        print(f"📦 Retrieved {len(rows):,} rows at offset {offset} (in {time.time()-t0:.2f}s)")

        chunk_size = max(1, len(rows) // n_processes)
        chunks = [rows[i:i+chunk_size] for i in range(0, len(rows), chunk_size)]
        with mp.Pool(n_processes) as pool:
            parsed_lists = pool.starmap(_parse_fp_chunk,
                                        [(c, need_sdf) for c in chunks])
        parsed_chunk = [item for sublist in parsed_lists for item in sublist]

        args = [(parsed_chunk, target_caes, threshold)]
        with mp.Pool(n_processes) as pool:
            batch_updates = pool.map(compare_one_formula, args)[0]
        for tkey, matches in batch_updates:
            best_matches[tkey].extend(matches)
            if matches and tkey in unmatched_keys:
                unmatched_keys.remove(tkey)

        for tdata in target_caes:
            tkey = fingerprint_key(tdata["fp"])
            current_best = distorted_match_scores.get(tkey, (-1, None))
            best_score, best_id = current_best
            for pop_row in parsed_chunk:
                pop_fp = pop_row["fp"]
                if not pop_fp:
                    continue
                score = sim.compute_similarity_score(tdata["fp"], pop_fp)
                if score > best_score:
                    distorted_match_scores[tkey] = (score, pop_row["cae_id"])

        print(f"🔎 Matched {len(target_keys) - len(unmatched_keys)}/{len(target_keys)}")
        if not unmatched_keys:
            print("✅ All CAEs matched (above threshold).")
            break

        offset += chunk_limit

    con.close()

    for tdata in target_caes:
        tkey = fingerprint_key(tdata["fp"])
        matches = best_matches.get(tkey, [])
        if matches:
            matches.sort(key=lambda x: -x[0])
            best_matches[tkey] = matches[:3]
        else:
            distorted_score, distorted_id = distorted_match_scores.get(tkey, (-1, None))
            if distorted_id is not None:
                best_matches[tkey] = [(distorted_score, distorted_id)]

    return best_matches

def compare_one_formula(args):
    population_chunk, target_data_list, threshold = args

    is_biatomic = target_data_list[0]["n_atoms"] == 2
    best = defaultdict(lambda: (-1.0, None))

    for pop_row in population_chunk:
        pop_fp  = pop_row["fp"]
        pop_id  = pop_row["cae_id"]
        pop_sdf = pop_row.get("sdf")   

        if not pop_fp:
            continue

        for tdata in target_data_list:
            tkey = fingerprint_key(tdata["fp"])

            if is_biatomic:
                try:
                    d_diff = abs(
                        interatomic_distance(tdata["sdf"]) -
                        interatomic_distance(pop_sdf)
                    )
                    score = 1.0 - d_diff if d_diff <= 0.1 else -1.0
                except Exception:
                    score = -1.0
            else:
                score = sim.compute_similarity_score(tdata["fp"], pop_fp)

            if score > best[tkey][0]:
                best[tkey] = (score, pop_id)

    return [(tkey, [best[tkey]]) for tkey in best]

def run_analysis(entry_id, population_file, idx,
                 db_path="all_caes_test.duckdb",
                 target_threshold=0.999,
                 comparison_threshold=0.99,
                 threshold_str=str,):

    with_charges_targets = {'ABAYAF', 'RULJAM'}

    start = time.time()
    print(f"\n🔍 Target: {entry_id}")
    output_dir = f"cae_comparison_results_{threshold_str}"
    os.makedirs(output_dir, exist_ok=True)  

    t0 = time.time()
    target_caes_dict = process_target(
        entry_id,
        threshold=target_threshold,
        with_charges=(entry_id in with_charges_targets)
    )
    t1 = time.time()
    total_caes = sum(len(v) for v in target_caes_dict.values())
    print(f"✅ Unique target CAEs: {total_caes} (retrieved in {t1 - t0:.2f}s)")
    print(f"✅ {len(target_caes_dict)} formulas: {list(target_caes_dict.keys())}")

    t2 = time.time()
    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f if line.strip()]
    t3 = time.time()
    print(f"⚙️ Starting population: {len(pop_ids)} molecules (loaded in {t3 - t2:.2f}s)")

    all_formulas_best = {}

    for formula, tcaes in target_caes_dict.items():
        if not tcaes:
            continue
        print(f"⚙️ Comparing formula: {formula} ({len(tcaes)} targets)")
        t4 = time.time()
        best_matches = compare_formula_streaming(
            tcaes, db_path, pop_ids, formula,
            threshold=comparison_threshold, n_processes=6, chunk_limit=600_000, 
            with_charges=(entry_id in with_charges_targets)
        )
        all_formulas_best[formula] = best_matches
        t5 = time.time()
        print(f"✅ Finished comparisons for {formula} (in {t5 - t4:.2f}s)")

    t6 = time.time()
    matched_caes = 0
    distorted_matches = 0
    no_match_in_db = 0
    total_caes = 0

    for formula, tcaes in target_caes_dict.items():
        bests = all_formulas_best.get(formula, {})
        for tdata in tcaes:
            total_caes += 1
            tkey = fingerprint_key(tdata["fp"])
            matches = bests.get(tkey, [])
            if not matches:
                distorted_matches += 1  
                continue

            top_score, cae_id = matches[0]
            is_biatomic = (tdata["n_atoms"] == 2)

            if cae_id == "NO_MATCH_IN_DB":
                print(f"🟥 No population CAE with formula {formula} found for CAE {tkey[:6]}")
                no_match_in_db += 1
            elif is_biatomic and (1.0 - top_score) <= 0.1:
                matched_caes += 1
            elif not is_biatomic and top_score >= comparison_threshold:
                matched_caes += 1
            else:
                distorted_matches += 1
    t7 = time.time()

    t9 = time.time()
    needed_cae_ids = set()
    for formula, tcaes in target_caes_dict.items():
        for tdata in tcaes:
            tkey = fingerprint_key(tdata["fp"])
            top_matches = all_formulas_best[formula].get(tkey, [])
            for (score, fid) in top_matches:
                if fid != "NO_MATCH_IN_DB":
                    needed_cae_ids.add(fid)

    sdf_map = fetch_sdf_for_cae_ids(db_path, list(needed_cae_ids))
    out_index = 0
    missing_sdf_ids = []

    for formula, tcaes in target_caes_dict.items():
        for tdata in tcaes:
            out_index += 1
            tkey = fingerprint_key(tdata["fp"])
            matches = all_formulas_best[formula].get(tkey, [])
            output_path = os.path.join(output_dir, f"{idx}_{entry_id}_cae{out_index}_matches.sdf")

            with MoleculeWriter(output_path) as writer:
                target_mol = Molecule.from_string(tdata["sdf"], format="sdf")
                writer.write(target_mol)

                for (score, fid) in matches:
                    if fid == "NO_MATCH_IN_DB":
                        continue
                    pop_sdf = sdf_map.get(fid)
                    if not pop_sdf:
                        missing_sdf_ids.append(fid)
                        continue

                    match_mol = Molecule.from_string(pop_sdf, format="sdf")
                    match_entry = Entry.from_molecule(match_mol)

                    if len(match_mol.atoms) == 2 and len(target_mol.atoms) == 2:
                        try:
                            d1 = interatomic_distance(target_mol.to_string("sdf"))
                            d2 = interatomic_distance(match_mol.to_string("sdf"))
                            dist_diff = abs(d1 - d2)
                            match_entry.attributes["DistanceDifference"] = f"{dist_diff:.4f}"
                            if dist_diff > 0.1:
                                match_entry.attributes["DistortedMatch"] = "True"
                        except:
                            match_entry.attributes["DistanceDifference"] = "ERROR"
                    else:
                        match_entry.attributes["Similarity"] = f"{score:.4f}"
                        if score < comparison_threshold:
                            match_entry.attributes["DistortedMatch"] = "True"

                    writer.write_entry(match_entry)

    if missing_sdf_ids:
        print(f"⚠️ Warning: {len(missing_sdf_ids)} matched CAE IDs had no SDF and were skipped.")
    print(f"📁 Wrote match SDFs in {time.time() - t9:.2f}s")

    print(f"📊 Matched {matched_caes}/{total_caes} CAEs (above threshold).")
    print(f"🟡 {distorted_matches} matched only as Distorted Matches (below threshold).")
    print(f"🟥 {no_match_in_db} had no population CAE with the same formula.")
    print(f"✅ Done with target {entry_id} in {time.time() - start:.2f}s")

if __name__ == "__main__":
    overall_start = time.time()

    targets = [
        'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
        'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
        'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
        'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
        'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV',
        'ABAYAF', 'RULJAM' 
    ]

    for i, target in enumerate(targets, start=1):
        run_analysis(
            entry_id=target,
            population_file=f"{population_folder}/{i}_{target}_init_pop.txt",
            idx=i,
            db_path="all_caes.duckdb",
            target_threshold=0.999,
            comparison_threshold=0.98,
            threshold_str=threshold_str,
        )

    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")
