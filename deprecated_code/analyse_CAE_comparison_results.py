import os
import glob
import re
from ccdc.io import MoleculeReader

def extract_frag_index(filename):
    match = re.search(r'_frag(\d+)_matches\.sdf$', filename)
    return int(match.group(1)) if match else -1

def extract_field(sdf_text, field):
    pattern = rf"<{field}>\s*(.*?)\s*(?:\n|\r|\$\$\$\$)"
    match = re.search(pattern, sdf_text)
    return match.group(1).strip() if match else None

# --- Target list ---
targets = ['ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
           'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
           'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
           'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
           'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV']

folder_path = "frag_comparison_results_duckdb_4"

total_all = 0
matched_all = 0

print()

for idx, target in enumerate(targets, start=1):
    base_name = f"{idx}_{target}"
    target_file = os.path.join(folder_path, f"{base_name}_target_unique_fragments.sdf")

    if not os.path.exists(target_file):
        print(f"❌ Target file missing: {target_file}")
        continue

    print(f"🔍 Processing: {base_name}")

    frag_pattern = os.path.join(folder_path, f"{base_name}_frag*_matches.sdf")
    frag_files = sorted(glob.glob(frag_pattern), key=extract_frag_index)

    total_files = len(frag_files)
    matched_fragments = 0
    skipped_fragments = []
    fallback_fragments = []

    for frag_file in frag_files:
        with open(frag_file, 'r') as f:
            sdf_text = f.read()

        mols = list(MoleculeReader(frag_file))
        if len(mols) < 2:
            skipped_fragments.append(os.path.basename(frag_file))
            continue

        fallback_value = extract_field(sdf_text, "Fallback")
        if fallback_value and fallback_value.strip().lower() == "true":
            sim = extract_field(sdf_text, "Similarity")
            dist = extract_field(sdf_text, "DistanceDifference")
            fallback_fragments.append((os.path.basename(frag_file), sim or dist))
        else:
            matched_fragments += 1

    num_skipped = len(skipped_fragments)
    num_fallbacks = len(fallback_fragments)
    num_matches = matched_fragments

    # Fragment-level reporting
    if skipped_fragments:
        print(f"  🚫 Missing formula (1 mol only):")
        for f in skipped_fragments:
            print(f"     - {f}")

    if fallback_fragments:
        print(f"  ⚠️ Fallback used (not counted as match):")
        for fname, val in fallback_fragments:
            print(f"     - {fname} (value: {val})")

    print(f"  → Number of fragments: {total_files}")
    print(f"  → Number of matches: {num_matches}")
    print(f"  → Number of fallbacks: {num_fallbacks}")
    print(f"  → Number of missing formulas: {num_skipped}")
    match_rate = 100 * num_matches / total_files if total_files else 0
    print(f"  → Match rate: {match_rate:.2f}%\n")

    # Global counter
    total_all += total_files
    matched_all += num_matches

# --- Final summary ---
print("======= FINAL SUMMARY =======")
print(f"Total fragments: {total_all}")
print(f"Matched fragments: {matched_all}")
overall_match_rate = (100 * matched_all / total_all) if total_all else 0
print(f"Overall match rate: {overall_match_rate:.2f}%")
