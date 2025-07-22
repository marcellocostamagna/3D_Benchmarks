import os
import time
import numpy as np
import multiprocessing as mp
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

### ---------- Utility Functions ---------- ###

# Function to generate an array from a molecule
# EXTRA FEATURE: PROTONS
def get_p_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

# EXTRA FEATURES: PROTONS & FORMAL CHARGES
def get_p_q_array_from_ccdcmol(ccdcmol):
    atom_array = []
    for atom in ccdcmol.atoms:
        if atom.coordinates is not None:
            atom_array.append([atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number), atom.formal_charge])    
    atom_array = np.array(atom_array)
    atom_array -= np.mean(atom_array, axis=0)  # Center the data
    return atom_array

def fingerprint_key(fp_obj):
    return tuple(fp_obj.tolist() if isinstance(fp_obj, np.ndarray) else fp_obj)

def formula_signature(fragment):
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    n_atoms = len(atoms)
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, n_atoms, formula)

def generate_fp_data(fragment, with_charges=False):
    arr = get_p_q_array_from_ccdcmol(fragment) if with_charges else get_p_array_from_ccdcmol(fragment)
    fp_data = fp.generate_fingerprint_from_data(arr)
    if isinstance(fp_data, np.ndarray):
        fp_data = fp_data.tolist()

    central, n_atoms, brute_formula = formula_signature(fragment)
    return {
        "sdf": fragment.to_string("sdf"),
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

### ---------- Fragmentation ---------- ###

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

def create_fragment(central_atom):
    frag = Molecule(identifier=f"{central_atom.label}_frag")
    atom_map = {central_atom: frag.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = frag.add_atom(neighbor)
    for bond in central_atom.bonds:
        a1, a2 = bond.atoms
        if a1 in atom_map and a2 in atom_map:
            try:
                frag.add_bond(bond.bond_type, atom_map[a1], atom_map[a2])
            except:
                pass
    return frag

def get_fragments(mol):
    return [create_fragment(atom) for atom in mol.atoms]

### ---------- Target Processing ---------- ###

def process_target(entry_id, threshold=0.999, with_charges=False):
    reader = io.EntryReader("CSD")
    mol = component_of_interest(reader.entry(entry_id).molecule)
    fragments = get_fragments(mol)
    grouped = defaultdict(list)
    for frag in fragments:
        data = generate_fp_data(frag, with_charges=with_charges)
        grouped[data["formula"]].append(data)
    unique_by_formula = defaultdict(list)
    for formula, frags in grouped.items():
        for frag in frags:
            if all(sim.compute_similarity_score(frag["fp"], other["fp"]) < threshold
                   for other in unique_by_formula[formula]):
                unique_by_formula[formula].append(frag)
    return unique_by_formula


### ---------- DB Query and Parsing ---------- ###

def _parse_fp_chunk(rows, with_sdf: bool):
    """Turn `fetchall()` rows into dicts (optionally with the SDF).

    Parameters
    ----------
    rows      : list[tuple]
        Result rows straight from DuckDB.
    with_sdf  : bool
        True ⟹ each row is (frag_id, fp_p[_pq], sdf)
        False ⟹ each row is (frag_id, fp_p[_pq])
    """
    parsed = []
    if with_sdf:                               # 3‑column row
        for frag_id, fp_str, sdf in rows:
            parsed.append(
                {"frag_id": frag_id,
                 "fp"     : _safe_fp_list(fp_str),
                 "sdf"    : sdf}
            )
    else:                                      # 2‑column row
        for frag_id, fp_str in rows:
            parsed.append(
                {"frag_id": frag_id,
                 "fp"     : _safe_fp_list(fp_str)}
            )
    return parsed


def _safe_fp_list(fp_str):
    """Turn the string representation of a fingerprint into a list of floats."""
    fp_str = str(fp_str).strip()
    if fp_str.startswith("[") and fp_str.endswith("]"):
        try:
            return [float(x) for x in fp_str.strip("[] ").split(",")]
        except ValueError:
            pass
    return []


def fetch_sdf_for_frag_ids(db_path, frag_ids):
    if not frag_ids:
        return {}
    con = duckdb.connect(db_path)
    con.execute("CREATE TEMP TABLE tmp_ids(frag_id INT)")
    con.executemany("INSERT INTO tmp_ids VALUES (?)", [(i,) for i in frag_ids])
    query = """
        SELECT fragments.rowid AS frag_id, fragments.sdf
        FROM fragments
        JOIN tmp_ids ON fragments.rowid = tmp_ids.frag_id
    """
    df = con.execute(query).fetchdf()
    con.close()
    return {row["frag_id"]: row["sdf"] for _, row in df.iterrows()}

### ---------- Chunked Streaming Comparison ---------- ###
def compare_formula_streaming(
        target_frags, db_path, pop_ids, formula,
        threshold=0.999, n_processes=4, chunk_limit=600_000,
        with_charges=False):

    # --- decide once whether we really need the SDF column -------------
    (central_atom, n_atoms, formula_str) = formula
    con = duckdb.connect(db_path)
    con.execute("CREATE TEMP TABLE tmp_pop_ids(entry_id TEXT)")
    con.executemany("INSERT INTO tmp_pop_ids VALUES (?)", [(eid,) for eid in pop_ids])
    need_sdf = target_frags[0]["n_atoms"] == 2          # bi‑atomic only
    fp_column = "fp_pq" if with_charges else "fp_p"
    select_cols = f"f.rowid, f.{fp_column}" + (", f.sdf" if need_sdf else "")

    base_query = f"""
        SELECT {select_cols}
        FROM fragments f
        JOIN tmp_pop_ids p ON f.entry_id = p.entry_id
        WHERE f.central_atom = '{central_atom}'
          AND f.n_atoms     = {n_atoms}
          AND f.formula_str = '{formula_str}'
    """

    offset = 0
    target_keys = {fingerprint_key(t["fp"]) for t in target_frags}
    unmatched_keys = set(target_keys)
    best_matches = defaultdict(list)
    fallback_scores = {}

    while True:
        chunk_query = base_query + f" LIMIT {chunk_limit} OFFSET {offset}"
        t0 = time.time()
        rows = con.execute(chunk_query).fetchall()

        # 🛑 Handle the case where nothing is found at all
        if offset == 0 and not rows:
            print(f"⚠️ No fragments found in DB for formula {formula}")
            for tdata in target_frags:
                tkey = fingerprint_key(tdata["fp"])
                fallback_scores[tkey] = (-1, "NO_MATCH_IN_DB")
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

        args = [(parsed_chunk, target_frags, threshold)]
        with mp.Pool(n_processes) as pool:
            batch_updates = pool.map(compare_one_formula, args)[0]
        for tkey, matches in batch_updates:
            best_matches[tkey].extend(matches)
            if matches and tkey in unmatched_keys:
                unmatched_keys.remove(tkey)

        for tdata in target_frags:
            tkey = fingerprint_key(tdata["fp"])
            current_best = fallback_scores.get(tkey, (-1, None))
            best_score, best_id = current_best
            for pop_row in parsed_chunk:
                pop_fp = pop_row["fp"]
                if not pop_fp:
                    continue
                score = sim.compute_similarity_score(tdata["fp"], pop_fp)
                if score > best_score:
                    fallback_scores[tkey] = (score, pop_row["frag_id"])

        print(f"🔎 Matched {len(target_keys) - len(unmatched_keys)}/{len(target_keys)}")
        if not unmatched_keys:
            print("✅ All fragments matched (above threshold).")
            break

        offset += chunk_limit

    con.close()

    # Final fallback handling
    for tdata in target_frags:
        tkey = fingerprint_key(tdata["fp"])
        matches = best_matches.get(tkey, [])
        if matches:
            matches.sort(key=lambda x: -x[0])
            best_matches[tkey] = matches[:3]
        else:
            fallback_score, fallback_id = fallback_scores.get(tkey, (-1, None))
            if fallback_id is not None:
                best_matches[tkey] = [(fallback_score, fallback_id)]

    return best_matches

### ---------- Comparison Helper ---------- ###

def compare_one_formula(args):
    population_chunk, target_data_list, threshold = args

    is_biatomic = target_data_list[0]["n_atoms"] == 2     
    best = defaultdict(lambda: (-1.0, None))            

    for pop_row in population_chunk:
        pop_fp  = pop_row["fp"]
        pop_id  = pop_row["frag_id"]
        pop_sdf = pop_row.get("sdf")   

        if not pop_fp:
            continue

        for tdata in target_data_list:
            tkey = fingerprint_key(tdata["fp"])

            # --- choose the scoring rule ---
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

            # keep the best score so far
            if score > best[tkey][0]:
                best[tkey] = (score, pop_id)

    return [(tkey, [best[tkey]]) for tkey in best]




### ---------- Main Analysis ---------- ###

def run_analysis(entry_id, population_file, idx,
                 db_path="all_fragments_3.duckdb",
                 target_threshold=0.999,
                 comparison_threshold=0.99):

    with_charges_targets = {'ABAYAF', 'RULJAM'}

    start = time.time()
    print(f"\n🔍 Target: {entry_id}")
    output_dir = "frag_comparison_results_duckdb_6"
    os.makedirs(output_dir, exist_ok=True)

    # --- 1) Generate target fragments ---
    t0 = time.time()
    target_frags_dict = process_target(
        entry_id,
        threshold=target_threshold,
        with_charges=(entry_id in with_charges_targets)
    )
    t1 = time.time()
    total_fragments = sum(len(v) for v in target_frags_dict.values())
    print(f"✅ Unique target fragments: {total_fragments} (retrieved in {t1 - t0:.2f}s)")
    print(f"✅ {len(target_frags_dict)} formulas: {list(target_frags_dict.keys())}")

    # --- 2) Load population IDs ---
    t2 = time.time()
    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f if line.strip()]
    t3 = time.time()
    print(f"⚙️ Initial population: {len(pop_ids)} molecules (loaded in {t3 - t2:.2f}s)")

    all_formulas_best = {}

    # --- 3) Per-formula streaming comparison ---
    for formula, tfrags in target_frags_dict.items():
        if not tfrags:
            continue
        print(f"⚙️ Comparing formula: {formula} ({len(tfrags)} targets)")
        t4 = time.time()
        best_matches = compare_formula_streaming(
            tfrags, db_path, pop_ids, formula,
            threshold=comparison_threshold, n_processes=8, chunk_limit=600_000, 
            with_charges=(entry_id in with_charges_targets)
        )
        all_formulas_best[formula] = best_matches
        t5 = time.time()
        print(f"✅ Finished comparisons for {formula} (in {t5 - t4:.2f}s)")

    # --- 4) Summary Counters ---
    t6 = time.time()
    matched_fragments = 0
    fallback_matches = 0
    no_match_in_db = 0
    total_fragments = 0

    for formula, tfrags in target_frags_dict.items():
        bests = all_formulas_best.get(formula, {})
        for tdata in tfrags:
            total_fragments += 1
            tkey = fingerprint_key(tdata["fp"])
            matches = bests.get(tkey, [])
            if not matches:
                fallback_matches += 1  # very rare
                continue

            top_score, frag_id = matches[0]
            is_biatomic = (tdata["n_atoms"] == 2)

            if frag_id == "NO_MATCH_IN_DB":
                print(f"🟥 No population fragment with formula {formula} found for fragment {tkey[:6]}")
                no_match_in_db += 1
            elif is_biatomic and (1.0 - top_score) <= 0.1:
                matched_fragments += 1
            elif not is_biatomic and top_score >= comparison_threshold:
                matched_fragments += 1
            else:
                fallback_matches += 1
    t7 = time.time()

    # --- 5) Write matched SDFs ---
    t9 = time.time()
    needed_frag_ids = set()
    for formula, tfrags in target_frags_dict.items():
        for tdata in tfrags:
            tkey = fingerprint_key(tdata["fp"])
            top_matches = all_formulas_best[formula].get(tkey, [])
            for (score, fid) in top_matches:
                if fid != "NO_MATCH_IN_DB":
                    needed_frag_ids.add(fid)

    sdf_map = fetch_sdf_for_frag_ids(db_path, list(needed_frag_ids))
    out_index = 0
    missing_sdf_ids = []

    for formula, tfrags in target_frags_dict.items():
        for tdata in tfrags:
            out_index += 1
            tkey = fingerprint_key(tdata["fp"])
            matches = all_formulas_best[formula].get(tkey, [])
            output_path = os.path.join(output_dir, f"{idx}_{entry_id}_frag{out_index}_matches.sdf")

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
                                match_entry.attributes["Fallback"] = "True"
                        except:
                            match_entry.attributes["DistanceDifference"] = "ERROR"
                    else:
                        match_entry.attributes["Similarity"] = f"{score:.4f}"
                        if score < comparison_threshold:
                            match_entry.attributes["Fallback"] = "True"

                    writer.write_entry(match_entry)

    if missing_sdf_ids:
        print(f"⚠️ Warning: {len(missing_sdf_ids)} matched fragment IDs had no SDF and were skipped.")
    print(f"📁 Wrote match SDFs in {time.time() - t9:.2f}s")

    # --- 6) Final Log Summary ---
    print(f"📊 Matched {matched_fragments}/{total_fragments} fragments (above threshold).")
    print(f"🟡 {fallback_matches} matched only with fallback (below threshold).")
    print(f"🟥 {no_match_in_db} had no population fragment with the same formula.")
    print(f"✅ Done with target {entry_id} in {time.time() - start:.2f}s")

### ---------- Main ---------- ###

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
            population_file=f"targets/init_populations_protons_0_5/{i}_{target}_init_pop.txt",
            idx=i,
            db_path="all_fragments_3.duckdb",
            target_threshold=0.999,
            comparison_threshold=0.98
        )

    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")