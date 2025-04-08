import os
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import defaultdict, Counter
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc import io
from ccdc.entry import Entry
from hsr import fingerprint as fp
from hsr import similarity as sim
import duckdb  # We will query our duckdb database.

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
    Return a tuple: (central_atom, n_atoms, formula_str)
    where formula_str is the counted formula for all atoms in 'fragment'.
    """
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    n_atoms = len(atoms)
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, n_atoms, formula)

def generate_fp_data(fragment):
    """
    Create a dict with 'sdf', 'fp', 'formula', 'n_atoms', 'central_atom' for the given CCDC fragment.
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

### ---------- DuckDB Loader ---------- ###

def load_population_subset_from_duckdb(pop_ids, formulas, db_path="all_fragments.duckdb"):
    """
    Connect to the DuckDB database `db_path`, and retrieve all rows from the `fragments` table
    that match:
      entry_id IN pop_ids
    AND
      (central_atom, n_atoms, formula_str) IN formulas

    Returns a dict: {formula: [list_of_frag_dicts], ...}
       where each frag_dict has { 'fp': [...], 'sdf': "...", 'n_atoms': int, 'formula': (central, n, form_str) }
    """

    # Convert pop_ids to a list or tuple for the SQL IN clause
    pop_ids = set(pid.strip().upper() for pid in pop_ids)
    # Convert formulas from e.g. [('C',4,'C2O2'), ...] to a set of tuples
    # We'll do a parameterized approach in SQL.

    # We'll store results in pop_group
    pop_group = defaultdict(list)

    if not pop_ids or not formulas:
        print("⚠️ No population IDs or formulas provided.")
        return pop_group

    # We'll build an IN clause or do a join approach. 
    # For many items, a join on a temp table can be more robust.
    # But for simplicity, let's do:
    #   central_atom, n_atoms, formula_str in (('C',4,'C2O2'), ...)
    # We'll need to build a big OR or a single "VALUES" join.
    # DuckDB can handle "IN" with a list of tuples, but let's do a small trick.

    # Connect to DuckDB
    con = duckdb.connect(db_path)

    # 1) We create a temporary table for pop_ids
    con.execute("CREATE TEMP TABLE tmp_pop_ids (entry_id TEXT)")
    pop_id_rows = [(eid,) for eid in pop_ids]
    # Use executemany to insert each row properly:
    con.executemany("INSERT INTO tmp_pop_ids VALUES (?)", pop_id_rows)

    # 2) We create a temporary table for formulas
    con.execute("""
    CREATE TEMP TABLE tmp_formulas(
        central_atom TEXT,
        n_atoms INT,
        formula_str TEXT
    )
    """)
    formula_rows = [(str(f[0]), int(f[1]), str(f[2])) for f in formulas]
    # Replace con.execute(..., formula_rows) with:
    con.executemany("INSERT INTO tmp_formulas VALUES (?, ?, ?)", formula_rows)


    # Now we can do a join:
    query = """
        SELECT f.central_atom, f.n_atoms, f.formula_str,
               f.entry_id, f.sdf, f.fp
        FROM fragments f
        JOIN tmp_pop_ids p ON f.entry_id = p.entry_id
        JOIN tmp_formulas tf
          ON f.central_atom = tf.central_atom
         AND f.n_atoms = tf.n_atoms
         AND f.formula_str = tf.formula_str
    """

    result_df = con.execute(query).fetchdf()
    con.close()

    # Now we parse each row into the format we want
    for _, row in result_df.iterrows():
        # The formula we want to group by
        form_tuple = (row["central_atom"], int(row["n_atoms"]), row["formula_str"])
        # Build the fragment dict
        # 'fp' is stored as text; parse it if it's a bracketed list
        # If it's something like "[0.1, 0.2, ...]", we parse:
        fp_str = str(row["fp"]).strip()
        try:
            # if it's bracketed, parse
            if fp_str.startswith("[") and fp_str.endswith("]"):
                fplist = [float(x) for x in fp_str.strip("[] ").split(",")]
            else:
                # fallback
                fplist = []
        except:
            fplist = []

        pop_group[form_tuple].append({
            "fp": fplist,
            "sdf": row["sdf"],
            "n_atoms": int(row["n_atoms"]),
            "formula": form_tuple
        })

    return pop_group

### ---------- Similarity Comparison ---------- ###

def compare_group(args):
    """
    Compare a target_group of fragments vs a pop_group (list of fragments).
    Keep top 3 matches by similarity or by distance difference if biatomic.
    """
    target_group, pop_group, threshold = args
    target_status = {
        fingerprint_key(t["fp"]): {
            "data": t,
            "top_matches": [],
            "matched": False
        } for t in target_group
    }
    unmatched_keys = set(target_status.keys())

    for pop in pop_group:
        pop_fp = pop["fp"]
        pop_sdf = pop["sdf"]
        pop_formula = pop["formula"]
        pop_atoms = pop["n_atoms"]

        to_remove = []

        for key in unmatched_keys:
            t = target_status[key]["data"]
            is_both_biatomic = (t["n_atoms"] == 2 and pop_atoms == 2)

            if is_both_biatomic and t["formula"] == pop_formula:
                # For 2-atom fragments, we compare distance difference
                try:
                    t_dist = interatomic_distance(t["sdf"])
                    p_dist = interatomic_distance(pop_sdf)
                    diff = abs(t_dist - p_dist)
                    if diff <= 0.01:
                        score = 1.0 - diff
                        target_status[key]["top_matches"].append((score, pop_sdf))
                        target_status[key]["matched"] = True
                        to_remove.append(key)
                except:
                    continue

            elif not is_both_biatomic:
                score = sim.compute_similarity_score(t["fp"], pop_fp)
                if score >= threshold:
                    target_status[key]["top_matches"].append((score, pop_sdf))
                    target_status[key]["matched"] = True
                    to_remove.append(key)

            # Keep top 3 matches
            target_status[key]["top_matches"].sort(key=lambda x: -x[0])
            target_status[key]["top_matches"] = target_status[key]["top_matches"][:3]

        unmatched_keys -= set(to_remove)
        if not unmatched_keys:
            break

    return [
        (key, {
            "target_sdf": info["data"]["sdf"],
            "top_matches": info["top_matches"],
            "matched": info["matched"]
        }) for key, info in target_status.items()
    ]

def compare_fragments_parallel(target_frags, pop_group, threshold=0.999, n_processes=8):
    """
    Compare each formula in target_frags vs the corresponding population list in pop_group.
    We do this in parallel if desired.
    """
    from time import perf_counter

    # Flatten tasks:  for each formula, we have a group of target fragments & population fragments.
    tasks = []
    for formula, tfrags in target_frags.items():
        # pop_group is a dict { formula: [list_of_pop_frags], ... }
        if formula in pop_group:
            tasks.append((tfrags, pop_group[formula], threshold))

    if not tasks:
        print("⚠️ No matching formulas found in the population.")
        return {}

    # We'll do parallel 'compare_group' calls
    print(f"🧩 Comparing {len(tasks)} formula groups in parallel with {n_processes} processes...")

    all_results = defaultdict(lambda: {
        "target_sdf": None,
        "top_matches": [],
        "matched": False
    })

    mp_start = perf_counter()
    with mp.Pool(n_processes) as pool:
        for batch in pool.imap_unordered(compare_group, tasks):
            for key, result in batch:
                # Merge results
                all_results[key]["target_sdf"] = result["target_sdf"]
                all_results[key]["top_matches"].extend(result["top_matches"])
                all_results[key]["matched"] |= result["matched"]
                # keep top 3
                all_results[key]["top_matches"].sort(key=lambda x: -x[0])
                all_results[key]["top_matches"] = all_results[key]["top_matches"][:3]

    mp_end = perf_counter()
    print(f"✅ Finished parallel comparison in {mp_end - mp_start:.2f}s")
    return dict(all_results)

### ---------- Main ---------- ###

def run_analysis(entry_id, population_file, idx,
                db_path="all_fragments.duckdb",
                threshold=0.999):
    """
    1) Build target fragments from 'entry_id'
    2) Load population IDs from population_file
    3) For each target fragment formula, retrieve matching fragments from DuckDB that have:
          - entry_id in pop_ids
          - same (central_atom, n_atoms, formula_str)
    4) Compare in parallel
    5) Write matched fragments to SDF
    """
    start = time.time()
    print(f"\n🔍 Target: {entry_id}")

    output_dir = "frag_comparison_results_duckdb"
    os.makedirs(output_dir, exist_ok=True)

    # 1) Generate target fragments
    t0 = time.time()
    target_frags = process_target(entry_id, threshold)
    t1 = time.time()
    total = sum(len(v) for v in target_frags.values())
    print(f"✅ Unique target fragments: {total} (retrieved in {t1 - t0:.2f}s)")
    print(f"✅ Unique formulas: {list(target_frags.keys())}")

    # 2) Load population IDs
    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f if line.strip()]
    print(f"⚙️ Population size: {len(pop_ids)}")

    # 3) Retrieve from DuckDB all fragments matching those formulas & pop IDs
    #    We'll pass all formula keys from target_frags to DuckDB
    formulas = list(target_frags.keys())
    pop_group = load_population_subset_from_duckdb(pop_ids, formulas, db_path=db_path)

    # 4) Compare in parallel
    comparisons = compare_fragments_parallel(target_frags, pop_group,
                                             threshold=threshold, n_processes=8)

    # Summaries
    matched = sum(1 for v in comparisons.values() if v["matched"])
    print(f"📊 Matched {matched}/{total} fragments from target {entry_id}.")

    # 5) Write out SDF files
    # 5a) Target SDF
    t5 = time.time()
    target_sdf_path = os.path.join(output_dir, f"{idx}_{entry_id}_target_unique_fragments.sdf")
    with MoleculeWriter(target_sdf_path) as w:
        for group in target_frags.values():
            for frag in group:
                w.write(Molecule.from_string(frag["sdf"], format="sdf"))
    print(f"🧪 Wrote target fragments SDF to {target_sdf_path} in {time.time() - t5:.2f}s")

    # 5b) Write matched SDFs
    t6 = time.time()
    for i, (fp_key, comp) in enumerate(comparisons.items(), start=1):
        output_path = os.path.join(output_dir, f"{idx}_{entry_id}_frag{i}_matches.sdf")
        with MoleculeWriter(output_path) as writer:
            target_mol = Molecule.from_string(comp["target_sdf"], format="sdf")
            writer.write(target_mol)

            for sim_score, sdf in comp["top_matches"]:
                match_mol = Molecule.from_string(sdf, format="sdf")
                match_entry = Entry.from_molecule(match_mol)
                if len(match_mol.atoms) == 2 and len(target_mol.atoms) == 2:
                    try:
                        d1 = interatomic_distance(target_mol.to_string("sdf"))
                        d2 = interatomic_distance(match_mol.to_string("sdf"))
                        match_entry.attributes["DistanceDifference"] = f"{abs(d1 - d2):.4f}"
                    except:
                        match_entry.attributes["DistanceDifference"] = "ERROR"
                else:
                    match_entry.attributes["Similarity"] = f"{sim_score:.4f}"
                writer.write_entry(match_entry)
    print(f"📁 Wrote match SDFs in {time.time() - t6:.2f}s")

    print(f"✅ Done with target {entry_id} in {time.time() - start:.2f}s")


if __name__ == "__main__":
    overall_start = time.time()

    # Example list of targets to analyze
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
            threshold=0.99
        )

    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")
