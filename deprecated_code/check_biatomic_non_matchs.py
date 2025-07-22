#!/usr/bin/env python3
import argparse
from pathlib import Path
import re

def count_atoms(sdf_text):
    """Count number of atoms in an SDF molecule block."""
    # Find the counts line: it is the 4th line of the molecule block
    lines = sdf_text.splitlines()
    if len(lines) < 4:
        return 0  # not a valid molecule
    try:
        # Standard SDF: atom count is characters 0–3 (may have leading spaces)
        count_line = lines[3]
        atom_count = int(count_line[0:3])
        return atom_count
    except Exception:
        return 0

def main():
    p = argparse.ArgumentParser(
        description="Find single-molecule biatomic .sdf files in a directory")
    p.add_argument(
        "--sdf-dir",
        type=Path,
        required=True,
        help="Directory containing .sdf files to scan")
    args = p.parse_args()

    sdf_files = sorted(args.sdf_dir.glob("*.sdf"))
    total_files = len(sdf_files)
    biatomic_nonmatches = []

    for sdf_path in sdf_files:
        text = sdf_path.read_text()
        # Only single-molecule files (one $$$$)
        if text.count("$$$$") == 1:
            # Extract just the first molecule (entire file up to $$$$)
            mol_block = text.split("$$$$")[0]
            atom_count = count_atoms(mol_block)
            if atom_count == 2:
                biatomic_nonmatches.append(sdf_path)

    print(f"Total .sdf files scanned:        {total_files}")
    print(f"Biatomic non-match files found:  {len(biatomic_nonmatches)}\n")
    if biatomic_nonmatches:
        print("Files with a single, biatomic molecule (non-matches):")
        for fp in biatomic_nonmatches:
            print(f"  {fp.name}")

if __name__ == "__main__":
    main()
