#!/usr/bin/env python3

"""
Scan a directory of SDF files for single-molecule biatomic CAE files (unmatched CAEs)
and files containing at least one biatomic CAE marked as a distorted match.

- Reports files with a single biatomic molecule (non-matches).
- Reports files where any biatomic CAE is labeled as a distorted match (DistortedMatch tag).
"""

import argparse
from pathlib import Path
import re

def count_atoms(mol_block):
    """Return number of atoms from the 4th line of the molecule block."""
    lines = mol_block.splitlines()
    if len(lines) < 4:
        return 0
    try:
        count_line = lines[3]
        atom_count = int(count_line[0:3])
        return atom_count
    except Exception:
        return 0

def main():
    p = argparse.ArgumentParser(
        description="Report single-molecule biatomic SDF files (CAEs with no matches) and files with biatomic distorted matches"
    )
    p.add_argument(
        "--sdf-dir",
        type=Path,
        required=True,
        help="Directory containing .sdf files to scan"
    )
    args = p.parse_args()

    sdf_files = sorted(args.sdf_dir.glob("*.sdf"))
    total_files = len(sdf_files)
    biatomic_nonmatches = []
    biatomic_distorted = []

    # regex to split molecule blocks
    re_mol_split = re.compile(r"\$\$\$\$")
    # regex to detect a DistortedMatch tag (case-insensitive)
    re_distorted_tag = re.compile(r"^> *<DistortedMatch>\s*$", re.MULTILINE)

    for sdf_path in sdf_files:
        text = sdf_path.read_text()
        mol_blocks = re_mol_split.split(text)
        mol_blocks = [mb for mb in mol_blocks if mb.strip()]
        # --- Find single-molecule, biatomic non-matches ---
        if len(mol_blocks) == 1:
            atom_count = count_atoms(mol_blocks[0])
            if atom_count == 2:
                biatomic_nonmatches.append(sdf_path.name)
        # --- Find biatomic distorted match blocks ---
        for mol_block in mol_blocks:
            if re_distorted_tag.search(mol_block):
                atom_count = count_atoms(mol_block)
                if atom_count == 2:
                    biatomic_distorted.append(sdf_path.name)
                    break  # Only report each file once for distorted

    print(f"Total .sdf files scanned:        {total_files}")
    print(f"Biatomic non-match files found:  {len(biatomic_nonmatches)}")
    print(f"SDF files with a biatomic DistortedMatch: {len(biatomic_distorted)}\n")

    if biatomic_nonmatches:
        print("Files with a single, biatomic molecule (non-matches):")
        for fp in biatomic_nonmatches:
            print(f"  {fp}")
    else:
        print("No biatomic non-match files found.")

    print()  # extra line

    if biatomic_distorted:
        print("Files containing at least one distorted match for a biatomic CAE:")
        for fp in biatomic_distorted:
            print(f"  {fp}")
    else:
        print("No files with a biatomic DistortedMatch found.")

if __name__ == "__main__":
    main()
