#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

def main():
    p = argparse.ArgumentParser(
        description="Scan .sdf files in a directory, report fallback stats, and file counts")
    p.add_argument(
        "--sdf-dir",
        type=Path,
        required=True,
        help="Directory containing .sdf files to scan")
    args = p.parse_args()

    # regex to detect a Fallback tag line
    re_fallback = re.compile(r"^> *<Fallback>\s*$", re.MULTILINE)
    # regex to extract the number following a <Similarity> tag
    re_sim = re.compile(r"^> *<Similarity>\s*\n([0-9]+\.[0-9]+)", re.MULTILINE)

    # gather all .sdf files
    sdf_files = sorted(args.sdf_dir.glob("*.sdf"))
    total_files = len(sdf_files)

    sims_by_file = {}
    single_mol_count = 0

    for sdf_path in sdf_files:
        text = sdf_path.read_text()
        # count how many molecule records by counting '$$$$'
        count_seps = text.count("$$$$")
        if count_seps == 1:
            single_mol_count += 1

        # only interested in ones with a <Fallback> tag
        if not re_fallback.search(text):
            continue

        # collect all similarity scores in that file
        sims = [float(m) for m in re_sim.findall(text)]
        if sims:
            sims_by_file[sdf_path] = sims

    fallback_file_count = len(sims_by_file)

    # flatten into list of (score, file) and sort ascending
    all_scores = [
        (score, fp)
        for fp, scores in sims_by_file.items()
        for score in scores
    ]
    all_scores.sort(key=lambda x: x[0])

    # Five lowest
    lowest5 = all_scores[:5]

    # Five highest & mean of those
    highest5 = sorted(all_scores, key=lambda x: x[0], reverse=True)[:5]
    mean_top5 = sum(score for score, _ in highest5) / len(highest5) if highest5 else 0.0

    # --- Output ---
    print(f"Total .sdf files scanned:           {total_files}")
    print(f"Files with one molecule NO MATCH:    {single_mol_count}")
    print(f"Files with FALLBACKS:         {fallback_file_count}\n")

    if fallback_file_count == 0:
        print("No .sdf files with <Fallback> found.")
        return

    print("Five lowest similarity scores (among those files):")
    for score, fp in lowest5:
        print(f"  {score:.4f}  →  {fp.name}")

    print("\nFive highest similarity scores (among those files):")
    for score, fp in highest5:
        print(f"  {score:.4f}  →  {fp.name}")

    print(f"\nMean of the top 5 highest scores:  {mean_top5:.4f}")

if __name__ == "__main__":
    main()
