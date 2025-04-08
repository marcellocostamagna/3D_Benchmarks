import os
import glob
import hsr
import ccdc
from ccdc.io import MoleculeReader
import numpy as np
import re

def extract_frag_index(filename):
    """Extract numeric frag index from filename."""
    match = re.search(r'_frag(\d+)_matches\.sdf$', filename)
    return int(match.group(1)) if match else -1

def extract_custom_field(sdf_text, field_name):
    """Extract value for a custom SDF tag like <Similarity> or <DistanceDifference>."""
    pattern = rf"<{field_name}>\s*([\s\S]*?)\s*\$\$\$\$"
    match = re.search(pattern, sdf_text)
    if match:
        try:
            return float(match.group(1).strip())
        except ValueError:
            return None
    return None

def get_array_from_ccdcmol(ccdcmol):
    coords = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return coords - coords.mean(axis=0)

# Define your targets
targets = ['ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
           'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
           'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
           'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
           'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV']

folder_path = "frag_comparison_results_duckdb"

# Stats summary
total_all = 0
matched_all = 0

# Iterate over indexed targets
for idx, target in enumerate(targets, start=1):
    base_name = f"{idx}_{target}"
    target_file = os.path.join(folder_path, f"{base_name}_target_unique_fragments.sdf")

    if not os.path.exists(target_file):
        print(f"Target file missing: {target_file}")
        continue

    print(f"\nProcessing: {base_name}")

    # Find fragment files
    frag_pattern = os.path.join(folder_path, f"{base_name}_frag*_matches.sdf")
    frag_files = sorted(glob.glob(frag_pattern))

    if not frag_files:
        print(f"  No fragment files found for {base_name}")
        continue

    frag_files_sorted = sorted(frag_files, key=extract_frag_index)
    total_comparable = 0
    match_count = 0
    skipped_count = 0

    for frag_idx, frag_file in enumerate(frag_files_sorted):
        with open(frag_file, "r") as f:
            sdf_text = f.read()

        # Determine which field is present
        if "<DistanceDifference>" in sdf_text:
            distance_field = extract_custom_field(sdf_text, "DistanceDifference")
            use_distance = True
            use_similarity = False
            similarity_field = None
        elif "<Similarity>" in sdf_text:
            similarity_field = extract_custom_field(sdf_text, "Similarity")
            use_distance = False
            use_similarity = True
            distance_field = None
        else:
            similarity_field = None
            distance_field = None
            use_distance = use_similarity = False

        # Load molecules
        frags = MoleculeReader(frag_file)
        list_mols = list(frags)

        if len(list_mols) < 2:
            skipped_count += 1
            continue

        total_comparable += 1
        frag = frags[1]
        trgt = frags[0]

        # Fingerprints
        frg_array = get_array_from_ccdcmol(frag)
        trgt_array = get_array_from_ccdcmol(trgt)
        frg_fp = hsr.fingerprint.generate_fingerprint_from_data(frg_array)
        trgt_fp = hsr.fingerprint.generate_fingerprint_from_data(trgt_array)

        sim = hsr.similarity.compute_similarity_score(frg_fp, trgt_fp)

        # Matching logic
        if use_distance and distance_field is not None:
            match_count += 1
        elif use_similarity and similarity_field is not None:
            # Log mismatches between HSR and SDF similarity
            if abs(sim - similarity_field) > 0.001:
                print(f"  → Fragment {frag_idx+1}: SDF sim = {similarity_field:.4f}, HSR sim = {sim:.4f}")
            if sim >= 0.9799:
                match_count += 1
        else:
            if sim >= 0.9799:
                match_count += 1

    # Target summary
    print(f"  → Compared {total_comparable} fragments.")
    print(f"  → Matches (sim > 0.98 or DistanceDifference): {match_count}")
    print(f"  → Skipped due to missing fragments: {skipped_count}")
    match_rate = 100 * match_count / total_comparable if total_comparable > 0 else 0
    print(f"  → Match rate: {match_rate:.2f}%")

    # Global stats
    total_all += total_comparable
    matched_all += match_count

# Final summary
print("\n======= FINAL SUMMARY =======")
print(f"Total fragments compared: {total_all}")
print(f"Fragments matched (via sim > 0.98 or DistanceDifference): {matched_all}")
print(f"Overall match rate: {(100 * matched_all / total_all) if total_all else 0:.2f}%")
