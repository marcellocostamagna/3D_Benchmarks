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
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    n_atoms = len(atoms)
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, n_atoms, formula)

def generate_fp_data(fragment):
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

def process_target(entry_id, threshold=0.999):
    reader = io.EntryReader("CSD")
    mol = component_of_interest(reader.entry(entry_id).molecule)
    fragments = get_fragments(mol)
    grouped = defaultdict(list)
    for frag in fragments:
        data = generate_fp_data(frag)
        grouped[data["formula"]].append(data)
    unique_by_formula = defaultdict(list)
    for formula, frags in grouped.items():
        for frag in frags:
            if all(sim.compute_similarity_score(frag["fp"], other["fp"]) < threshold
                   for other in unique_by_formula[formula]):
                unique_by_formula[formula].append(frag)
    return unique_by_formula

### ---------- DB Query and Parsing ---------- ###

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

def compare_formula_streaming(target_frags, db_path, pop_ids, formula, threshold=0.999, n_processes=4, chunk_limit=6_000_000):
    (central_atom, n_atoms, formula_str) = formula
    con = duckdb.connect(db_path)
    con.execute("CREATE TEMP TABLE tmp_pop_ids(entry_id TEXT)")
    con.executemany("INSERT INTO tmp_pop_ids VALUES (?)", [(eid,) for eid in pop_ids])
    base_query = f"""
        SELECT f.rowid AS frag_id, f.fp
        FROM fragments f
        JOIN tmp_pop_ids p ON f.entry_id = p.entry_id
        WHERE f.central_atom = '{central_atom}'
          AND f.n_atoms = {n_atoms}
          AND f.formula_str = '{formula_str}'
    """

    offset = 0
    target_keys = {fingerprint_key(t["fp"]) for t in target_frags}
    unmatched_keys = set(target_keys)
    best_matches = defaultdict(list)

    while unmatched_keys:
        chunk_query = base_query + f" LIMIT {chunk_limit} OFFSET {offset}"
        t0 = time.time()
        rows = con.execute(chunk_query).fetchall()
        if not rows:
            break
        print(f"📦 Retrieved {len(rows):,} rows at offset {offset} (in {time.time()-t0:.2f}s)")

        # Multiprocessing parse
        chunk_size = max(1, len(rows) // n_processes)
        chunks = [rows[i:i+chunk_size] for i in range(0, len(rows), chunk_size)]
        with mp.Pool(n_processes) as pool:
            results = pool.map(_parse_fp_chunk, chunks)
        parsed_chunk = [item for sublist in results for item in sublist]

        # Compare
        args = [(parsed_chunk, target_frags, threshold)]
        with mp.Pool(n_processes) as pool:
            batch_updates = pool.map(compare_one_formula, args)[0]
        for tkey, matches in batch_updates:
            best_matches[tkey].extend(matches)
            if matches and tkey in unmatched_keys:
                unmatched_keys.remove(tkey)

        # Trim top 3, or fallback to best available match if no match above threshold
        for tdata in target_frags:
            tkey = fingerprint_key(tdata["fp"])
            matches = best_matches.get(tkey, [])
            if matches:
                matches.sort(key=lambda x: -x[0])
                best_matches[tkey] = matches[:3]
            elif parsed_chunk:  # fallback logic
                # Compute best match manually
                best_match = None
                best_score = -1
                for pop_row in parsed_chunk:
                    score = sim.compute_similarity_score(tdata["fp"], pop_row["fp"])
                    if score > best_score:
                        best_score = score
                        best_match = pop_row["frag_id"]
                if best_match is not None:
                    best_matches[tkey] = [(best_score, best_match)]


        print(f"🔎 Matched {len(target_keys) - len(unmatched_keys)}/{len(target_keys)}")
        if not unmatched_keys:
            print("✅ All fragments matched.")
            break

        offset += chunk_limit

    con.close()
    return best_matches

### ---------- Comparison Helper ---------- ###

def compare_one_formula(args):
    population_chunk, target_data_list, threshold = args
    match_records = defaultdict(list)
    unmatched_fp_keys = {fingerprint_key(t["fp"]) for t in target_data_list}
    for pop_row in population_chunk:
        pop_fp = pop_row["fp"]
        pop_id = pop_row["frag_id"]
        pop_n_atoms = pop_row.get("n_atoms", None)
        if not unmatched_fp_keys:
            break
        for tdata in target_data_list:
            tkey = fingerprint_key(tdata["fp"])
            if tkey not in unmatched_fp_keys:
                continue
            t_n_atoms = tdata["n_atoms"]
            is_both_biatomic = (t_n_atoms == 2 and pop_n_atoms == 2)
            if is_both_biatomic:
                try:
                    t_dist = interatomic_distance(tdata["sdf"])
                    p_dist = interatomic_distance(pop_row["sdf"])
                    if abs(t_dist - p_dist) <= 0.01:
                        score = 1.0 - abs(t_dist - p_dist)
                        match_records[tkey].append((score, pop_id))
                        unmatched_fp_keys.remove(tkey)
                except:
                    pass
            else:
                score = sim.compute_similarity_score(tdata["fp"], pop_fp)
                if score >= threshold:
                    match_records[tkey].append((score, pop_id))
                    unmatched_fp_keys.remove(tkey)
    final_updates = []
    for tdata in target_data_list:
        tkey = fingerprint_key(tdata["fp"])
        matches = match_records[tkey]
        matches.sort(key=lambda x: -x[0])
        final_updates.append((tkey, matches[:3]))
    return final_updates

### ---------- Main Analysis ---------- ###

def run_analysis(entry_id, population_file, idx,
                 db_path="all_fragments.duckdb",
                 target_threshold=0.999,
                 comparison_threshold=0.99):
    start = time.time()
    print(f"\n🔍 Target: {entry_id}")
    output_dir = "frag_comparison_results_duckdb_3"
    os.makedirs(output_dir, exist_ok=True)

    # --- 1) Generate target fragments ---
    t0 = time.time()
    target_frags_dict = process_target(entry_id, target_threshold)    
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
            threshold=comparison_threshold, n_processes=8, chunk_limit=6_000_000
        )
        all_formulas_best[formula] = best_matches
        t5 = time.time()
        print(f"✅ Finished comparisons for {formula} (in {t5 - t4:.2f}s)")

  # --- 4) Summary ---
    t6 = time.time()
    matched_fragments = 0
    total_fragments = 0
    for formula, tfrags in target_frags_dict.items():
        bests = all_formulas_best.get(formula, {})
        for tdata in tfrags:
            total_fragments += 1
            tkey = fingerprint_key(tdata["fp"])
            if tkey in bests and bests[tkey]:  # check if the fragment has any match
                matched_fragments += 1
    t7 = time.time()
    print(f"📊 Matched {matched_fragments}/{total_fragments} fragments. (checked in {t7 - t6:.2f}s)")

    # --- 5a) Write target SDF ---
    t8 = time.time()
    target_sdf_path = os.path.join(output_dir, f"{idx}_{entry_id}_target_unique_fragments.sdf")
    with MoleculeWriter(target_sdf_path) as w:
        for group in target_frags_dict.values():
            for frag in group:
                w.write(Molecule.from_string(frag["sdf"], format="sdf"))
    print(f"🧪 Wrote target SDF to {target_sdf_path} in {time.time() - t8:.2f}s")

    # --- 5b) Write matched SDFs ---
    t9 = time.time()
    needed_frag_ids = set()
    for formula, tfrags in target_frags_dict.items():
        for tdata in tfrags:
            tkey = fingerprint_key(tdata["fp"])
            top_matches = all_formulas_best[formula].get(tkey, [])
            for (score, fid) in top_matches:
                needed_frag_ids.add(fid)

    sdf_map = fetch_sdf_for_frag_ids(db_path, list(needed_frag_ids))

    out_index = 0
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
                    pop_sdf = sdf_map.get(fid, None)
                    if pop_sdf:
                        match_mol = Molecule.from_string(pop_sdf, format="sdf")
                        match_entry = Entry.from_molecule(match_mol)
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
    print(f"📁 Wrote match SDFs in {time.time() - t9:.2f}s")
    print(f"✅ Done with target {entry_id} in {time.time() - start:.2f}s")

### ---------- Main ---------- ###

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
            target_threshold=0.999,
            comparison_threshold=0.98
        )

    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")
