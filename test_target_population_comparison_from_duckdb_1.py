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

def get_array_from_ccdcmol(ccdcmol):
    coords = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return coords - coords.mean(axis=0)

def fingerprint_key(fp_obj):
    return tuple(fp_obj.tolist() if isinstance(fp_obj, np.ndarray) else fp_obj)

def formula_signature(fragment):
    """
    Return a tuple: (central_atom, n_atoms, formula_str).
    """
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    n_atoms = len(atoms)
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, n_atoms, formula)

def generate_fp_data(fragment):
    """
    Return a dict with 'sdf', 'fp', 'formula', 'n_atoms', 'central_atom' for the given CCDC fragment.
    """
    return {
        "sdf": fragment.to_string("sdf"),
        "fp": fp.generate_fingerprint_from_data(get_array_from_ccdcmol(fragment)),
        "formula": formula_signature(fragment),
        "n_atoms": len(fragment.atoms),
        "central_atom": fragment.atoms[0].atomic_symbol,
    }

def interatomic_distance(sdf_string):
    mol = Molecule.from_string(sdf_string, format="sdf")
    coords = [atom.coordinates for atom in mol.atoms]
    return np.linalg.norm(np.array(coords[0]) - np.array(coords[1]))

### ---------- Fragmentation Functions ---------- ###

def component_of_interest(molecule):
    """
    Return the 'main' component for the given molecule, skipping salts/solvents.
    """
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
    """
    Create a single-atom-centered fragment for each atom in 'mol'.
    """
    return [create_fragment(atom) for atom in mol.atoms]

### ---------- Target Processing ---------- ###

def process_target(entry_id, threshold=0.999):
    """
    From a CSD entry_id, extract its main component, generate all single-atom-centered fragments,
    group them by formula, and keep only unique ones above a similarity threshold (0.999).
    """
    reader = io.EntryReader("CSD")
    mol = component_of_interest(reader.entry(entry_id).molecule)
    fragments = get_fragments(mol)

    grouped = defaultdict(list)
    for frag in fragments:
        data = generate_fp_data(frag)
        grouped[data["formula"]].append(data)

    # remove near-duplicates within the target
    unique_by_formula = defaultdict(list)
    for formula, frags in grouped.items():
        for frag in frags:
            if all(sim.compute_similarity_score(frag["fp"], other["fp"]) < threshold
                   for other in unique_by_formula[formula]):
                unique_by_formula[formula].append(frag)

    return unique_by_formula

### ---------- Database Query (Formula by Formula) ---------- ###

def _parse_fp_chunk(rows):
    parsed = []
    for frag_id, fp_str in rows:
        fp_str = str(fp_str).strip()
        if fp_str.startswith("[") and fp_str.endswith("]"):
            try:
                fplist = [float(x) for x in fp_str.strip("[] ").split(",")]
            except:
                fplist = []
        else:
            fplist = []
        parsed.append({"frag_id": frag_id, "fp": fplist})
    return parsed

def load_population_fps_for_formula(db_path, pop_ids, formula, target_fps=None, n_processes=4, chunk_limit=6_000_000):
    """
    Load fragment fingerprints from DuckDB for a specific formula, in chunks (LIMIT/OFFSET).
    - target_fps: optional list of target fragments (with 'fp' keys) to support early stopping
    """

    (central_atom, n_atoms, formula_str) = formula

    con = duckdb.connect(db_path)

    # Create TEMP table for pop_ids
    con.execute("CREATE TEMP TABLE tmp_pop_ids(entry_id TEXT)")
    con.executemany("INSERT INTO tmp_pop_ids VALUES (?)", [(eid,) for eid in pop_ids])

    base_query = f"""
        SELECT
            f.rowid AS frag_id,  
            f.fp
        FROM fragments f
        JOIN tmp_pop_ids p ON f.entry_id = p.entry_id
        WHERE
            f.central_atom = '{central_atom}'
            AND f.n_atoms = {n_atoms}
            AND f.formula_str = '{formula_str}'
    """

    offset = 0
    all_results = []
    matched_keys = set()
    target_keys = set()

    if target_fps:
        target_keys = {fingerprint_key(t["fp"]) for t in target_fps}

    print(f"🔎 Executing DuckDB chunked query for {formula_str}...")
    while True:
        t_chunk_start = time.time()
        chunk_query = base_query + f" LIMIT {chunk_limit} OFFSET {offset}"
        rows = con.execute(chunk_query).fetchall()

        if not rows:
            break

        print(f"📦 Retrieved {len(rows):,} rows at offset {offset} (in {time.time() - t_chunk_start:.2f}s)")

        # Multiprocessing parse
        t_parse = time.time()
        chunk_size = max(1, len(rows) // n_processes)
        chunks = [rows[i:i+chunk_size] for i in range(0, len(rows), chunk_size)]
        with mp.Pool(n_processes) as pool:
            results = pool.map(_parse_fp_chunk, chunks)
        parsed = [item for sublist in results for item in sublist]
        print(f"✅ Parsed {len(parsed):,} fingerprints in {time.time() - t_parse:.2f}s")

        all_results.extend(parsed)

        # Early exit check
        if target_fps:
            for item in parsed:
                k = fingerprint_key(item["fp"])
                if k in target_keys:
                    matched_keys.add(k)
            unmatched = target_keys - matched_keys
            print(f"🔎 Matched {len(matched_keys)}/{len(target_keys)} target fragments so far.")
            if not unmatched:
                print("✅ All target fragments matched — stopping early.")
                break

        offset += chunk_limit

    con.close()
    print(f"🧹 Total parsed fragments: {len(all_results):,}")
    return all_results


# def load_population_fps_for_formula(db_path, pop_ids, formula):
#     """
#     ### NEW or CHANGED ###
#     Query only rowid (unique id) and fingerprint from the DB for a single formula.
#     We skip the SDF here to save memory.

#     Return a list of dicts:
#         [ { 'frag_id': rowid?, 'fp': [floats], ... }, ... ]
#     so that we can do fingerprint comparisons.
#     """
#     # pop_ids = set(pid.strip().upper() for pid in pop_ids)
#     # if not pop_ids:
#     #     return []

#     (central_atom, n_atoms, formula_str) = formula

#     con = duckdb.connect(db_path)

#     # Create TEMP table for pop_ids
#     con.execute("CREATE TEMP TABLE tmp_pop_ids(entry_id TEXT)")
#     pop_id_rows = [(eid,) for eid in pop_ids]
#     con.executemany("INSERT INTO tmp_pop_ids VALUES (?)", pop_id_rows)

#     query = f"""
#         SELECT
#         f.rowid AS frag_id,  
#         f.fp
#     FROM fragments f
#     JOIN tmp_pop_ids p ON f.entry_id = p.entry_id
#     WHERE
#         f.central_atom = '{central_atom}'
#         AND f.n_atoms = {n_atoms}
#         AND f.formula_str = '{formula_str}'
#     """
#     df = con.execute(query).fetchdf()
#     con.close()

#     results = []
#     for _, row in df.iterrows():
#         frag_id = row["frag_id"]
#         fp_str = str(row["fp"]).strip()
#         if fp_str.startswith("[") and fp_str.endswith("]"):
#             try:
#                 fplist = [float(x) for x in fp_str.strip("[] ").split(",")]
#             except:
#                 fplist = []
#         else:
#             fplist = []
#         results.append({"frag_id": frag_id, "fp": fplist})
#     return results

def fetch_sdf_for_frag_ids(db_path, frag_ids):
    """
    Given a list of frag_id (rowid), fetch corresponding SDF from the database.
    Return { frag_id: sdf_string }.
    """
    if not frag_ids:
        return {}

    con = duckdb.connect(db_path)

    # Create a temp table of IDs
    con.execute("CREATE TEMP TABLE tmp_ids(frag_id INT)")
    con.executemany("INSERT INTO tmp_ids VALUES (?)", [(i,) for i in frag_ids])

    query = """
        SELECT fragments.rowid AS frag_id, fragments.sdf
        FROM fragments
        JOIN tmp_ids ON fragments.rowid = tmp_ids.frag_id
    """
    df = con.execute(query).fetchdf()
    con.close()

    sdf_map = {row["frag_id"]: row["sdf"] for _, row in df.iterrows()}
    return sdf_map


### ---------- Similarity Comparison (per formula) ---------- ###

def compare_one_formula(args):
    """
    This function is used by multiprocessing to compare a chunk of population
    fingerprints against *all* target fragments for one formula, but only until
    those target fragments are matched.

    We re-introduce a biatomic distance comparison:
      If both fragments have 2 atoms, we compute the interatomic distance difference
      (<= 0.01 => match). Otherwise, we do the normal fingerprint similarity check.

    Returns: a list of (target_fp_key, [ (score, frag_id), ...topN... ]) updates
             so we can merge them back (i.e. store top-3 matches).
    """
    population_chunk, target_data_list, threshold = args

    # Dictionary: target_fp_key -> list of (score, frag_id)
    match_records = defaultdict(list)

    # Keep track of which target fragments (by fp_key) are still unmatched
    unmatched_fp_keys = {fingerprint_key(t["fp"]) for t in target_data_list}

    # Loop over each population fragment in this chunk
    for pop_row in population_chunk:
        pop_fp = pop_row["fp"]
        pop_id = pop_row["frag_id"]
        pop_n_atoms = pop_row.get("n_atoms", None)  # If you store n_atoms in pop_row

        # If all target fragments are matched, we can stop
        if not unmatched_fp_keys:
            break

        # Compare to each target fragment
        for tdata in target_data_list:
            tkey = fingerprint_key(tdata["fp"])
            if tkey not in unmatched_fp_keys:
                # Already matched in a previous iteration
                continue

            t_n_atoms = tdata["n_atoms"]

            # --- Reintroduce the distance-based comparison for 2-atom fragments --- #
            is_both_biatomic = (t_n_atoms == 2 and pop_n_atoms == 2)
            if is_both_biatomic:
                # We need the SDF for both fragments, so presumably pop_row also has "sdf"
                # or we fetch it from DB. Below assumes pop_row["sdf"] is available.
                try:
                    t_dist = interatomic_distance(tdata["sdf"])
                    p_dist = interatomic_distance(pop_row["sdf"])  # must exist in pop_row
                    diff = abs(t_dist - p_dist)
                    if diff <= 0.01:
                        # We'll store a "score" of (1 - diff) so that closer => higher score
                        score = 1.0 - diff
                        match_records[tkey].append((score, pop_id))
                        unmatched_fp_keys.remove(tkey)
                except:
                    pass
            else:
                # Normal similarity check for non-biatomic or if n_atoms is missing
                score = sim.compute_similarity_score(tdata["fp"], pop_fp)
                if score >= threshold:
                    match_records[tkey].append((score, pop_id))
                    unmatched_fp_keys.remove(tkey)
            # ----------------------------------------------------------------------- #
            #
            # We do NOT break after matching, because other target fragments
            # in this chunk might also need to compare with the same pop fragment.
            #

    # Now keep only top 3 matches per target
    final_updates = []
    for tdata in target_data_list:
        tkey = fingerprint_key(tdata["fp"])
        matches = match_records[tkey]
        # Sort descending by the "score" (similarity or 1 - distance)
        matches.sort(key=lambda x: -x[0])
        matches = matches[:3]
        final_updates.append((tkey, matches))

    return final_updates

def compare_formula_in_parallel(target_frags, population, threshold=0.999, n_processes=4):
    """
    ### NEW or CHANGED ###
    We take all target fragments that share the same formula (the 'target_frags' list),
    and compare them to the entire population list for that formula, in parallel chunks.

    We short-circuit as soon as all target fragments have found a match, but only
    after we've finished a chunk. (For finer-grained short-circuiting, you'd need
    shared memory or manager logic. This approach is simpler.)
    """
    if not target_frags or not population:
        return {}

    # We'll define how big each chunk is. Larger chunk => fewer calls but slower short-circuit.
    chunk_size = 5000
    # Create chunks
    chunks = [population[i:i+chunk_size] for i in range(0, len(population), chunk_size)]

    # Prepare a dict to store final best matches per target fragment
    # Key = fingerprint_key(t["fp"])
    # Value = list of (score, frag_id)
    best_matches = defaultdict(list)

    # We'll track how many are still unmatched
    unmatched_tkeys = set(fingerprint_key(t["fp"]) for t in target_frags)

    # We'll do multiprocess with e.g. n_processes
    with mp.Pool(n_processes) as pool:
        for chunk_idx, chunk_data in enumerate(chunks):
            if not unmatched_tkeys:
                # All matched => short-circuit
                break

            # Prepare arguments for compare_one_formula
            # We pass the entire target_frags, but it might be more memory
            # If that is large, we can pass just the unmatched subset
            args_list = [(chunk_data, target_frags, threshold)]  # just 1 task for 1 chunk

            results = pool.map(compare_one_formula, args_list)

            # results is a list with 1 item => a list of (tkey, [ (score, frag_id), ... ])
            # Merge them
            for batch_updates in results:
                for (tkey, new_matches) in batch_updates:
                    if new_matches:
                        # we combine with existing
                        best_matches[tkey].extend(new_matches)
                        # If we found matches, it means that tkey is matched
                        # But we may want top3 merges. Let's do that after
                        # or we just note that it's matched. Then remove from unmatched
                        if new_matches:
                            # remove from unmatched
                            if tkey in unmatched_tkeys:
                                unmatched_tkeys.remove(tkey)

            # After chunk, let's reduce top-3
            for tkey in best_matches:
                best_matches[tkey].sort(key=lambda x: -x[0])
                best_matches[tkey] = best_matches[tkey][:3]

    return best_matches

### ---------- Main Analysis Loop ---------- ###

def run_analysis(entry_id, population_file, idx,
                db_path="all_fragments.duckdb",
                threshold=0.999):
    start = time.time()
    print(f"\n🔍 Target: {entry_id}")
    output_dir = "frag_comparison_results_duckdb"
    os.makedirs(output_dir, exist_ok=True)

    # 1) Generate target fragments
    t0 = time.time()
    target_frags_dict = process_target(entry_id, threshold)    
    t1 = time.time()
    total_fragments = sum(len(v) for v in target_frags_dict.values())
    print(f"✅ Unique target fragments: {total_fragments} (retrieved in {t1 - t0:.2f}s)")
    print(f"✅ {len(target_frags_dict)} Unique formulas: {list(target_frags_dict.keys())}")

    # 2) Load population IDs
    t2 = time.time()
    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f if line.strip()]
    t3 = time.time()
    print(f"⚙️ Initial population: {len(pop_ids)} molecules (loaded in {t3 - t2:.2f}s)")

    # We’ll collect final top matches for *all* formulas in a structure:
    # { (formula): { tkey: [ (score, frag_id), ... ] } }
    all_formulas_best = {}

    # Now we handle each formula separately
    # For each formula in the target
    for formula, tfrags in target_frags_dict.items():
        
        if not tfrags:
            continue

        # 3) Load the population (only rowid + fp) that match this formula
        t4 = time.time()
        print(f'⚙️ Loading matching fragments fingerprints for: {formula}')
        population_fps = load_population_fps_for_formula(db_path, pop_ids, formula, n_processes=8)
        t5 = time.time()
        print(f'✅ Loaded {len(population_fps)} fingerprints ({t5 - t4:.2f}s)')
        if not population_fps:
            print(f"   No population fragments found for formula: {formula}")
            all_formulas_best[formula] = {}
            continue

        print(f"   >>> Processing formula={formula}, #fps={len(population_fps)}, #target={len(tfrags)}")
        t6 = time.time()
        # 4) Compare in parallel, short-circuit as soon as all matched
        best_matches_for_formula = compare_formula_in_parallel(
            tfrags, population_fps, threshold=0.99, n_processes=8
        )
        all_formulas_best[formula] = best_matches_for_formula
        t7 = time.time()
        print(f"   ✅ Finished comparisons for {formula} (in {t7 - t6:.2f}s)")
        

    # Summaries: how many got matched
    t8 = time.time()
    matched_fragments = 0
    for formula, tfrags in target_frags_dict.items():
        for tdata in tfrags:
            tkey = fingerprint_key(tdata["fp"])
            if tkey in all_formulas_best[formula]:
                if all_formulas_best[formula][tkey]:
                    matched_fragments += 1
    t9 = time.time()
    print(f"📊 Matched {matched_fragments}/{total_fragments} fragments from target {entry_id}. ({t9 - t8:.2f}s)")

    # 5) Write out SDF files
    # 5a) Write all target fragments
    t10 = time.time()
    target_sdf_path = os.path.join(output_dir, f"{idx}_{entry_id}_target_unique_fragments.sdf")
    with MoleculeWriter(target_sdf_path) as w:
        for group in target_frags_dict.values():
            for frag in group:
                w.write(Molecule.from_string(frag["sdf"], format="sdf"))
    print(f"🧪 Wrote target fragments SDF to {target_sdf_path} in {time.time() - t10:.2f}s")

    # 5b) Write matched SDFs
    # We need to fetch the SDF for each matched frag_id
    # Then we can build a result SDF per target fragment
    t11 = time.time()

    # Build a big set of all needed frag_ids
    needed_frag_ids = set()
    for formula, tfrags in target_frags_dict.items():
        for tdata in tfrags:
            tkey = fingerprint_key(tdata["fp"])
            top_matches = all_formulas_best[formula].get(tkey, [])
            for (score, fid) in top_matches:
                needed_frag_ids.add(fid)

    # Now fetch SDF for all needed frag_ids
    # (This is a single query if your DB can handle it.)
    sdf_map = fetch_sdf_for_frag_ids(db_path, list(needed_frag_ids))

    # Next, write out a separate SDF for each target fragment
    out_index = 0
    for formula, tfrags in target_frags_dict.items():
        for tdata in tfrags:
            out_index += 1
            tkey = fingerprint_key(tdata["fp"])
            matches = all_formulas_best[formula].get(tkey, [])
            # We'll create a file "fragX_matches.sdf"
            output_path = os.path.join(
                output_dir,
                f"{idx}_{entry_id}_frag{out_index}_matches.sdf"
            )
            with MoleculeWriter(output_path) as writer:
                # Write the target
                target_mol = Molecule.from_string(tdata["sdf"], format="sdf")
                writer.write(target_mol)

                # For each matched population, write
                for (score, fid) in matches:
                    pop_sdf = sdf_map.get(fid, None)
                    if pop_sdf is None:
                        continue
                    match_mol = Molecule.from_string(pop_sdf, format="sdf")
                    match_entry = Entry.from_molecule(match_mol)

                    # If both are diatomic, do distance difference check
                    if len(match_mol.atoms) == 2 and len(target_mol.atoms) == 2:
                        try:
                            d1 = interatomic_distance(target_mol.to_string("sdf"))
                            d2 = interatomic_distance(match_mol.to_string("sdf"))
                            match_entry.attributes["DistanceDifference"] = f"{abs(d1 - d2):.4f}"
                        except:
                            match_entry.attributes["DistanceDifference"] = "ERROR"
                    else:
                        match_entry.attributes["Similarity"] = f"{score:.4f}"
                    writer.write_entry(match_entry)

    print(f"📁 Wrote match SDFs in {time.time() - t11:.2f}s")
    print(f"✅ Done with target {entry_id} in {time.time() - start:.2f}s")

### ---------- Example Main ---------- ###
if __name__ == "__main__":
    overall_start = time.time()

    targets = [
        'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
        'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
        'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
        'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
        'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV'
    ]

    for i, target in enumerate(targets, start=1):
        run_analysis(
            entry_id=target,
            population_file=f"targets/init_populations_protons_2/{i}_{target}_init_pop.txt",
            idx=i,
            db_path="all_fragments.duckdb",
            threshold=0.999
        )

    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")
