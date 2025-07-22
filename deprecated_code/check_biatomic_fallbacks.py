#!/usr/bin/env python3
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
        description="Report single-molecule biatomic SDF files (non-matches) and files with biatomic fallbacks"
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
    biatomic_fallbacks = []

    # regex to split molecule blocks
    re_mol_split = re.compile(r"\$\$\$\$")
    # regex to detect a Fallback tag
    re_fallback_tag = re.compile(r"^> *<Fallback>\s*$", re.MULTILINE)

    for sdf_path in sdf_files:
        text = sdf_path.read_text()
        mol_blocks = re_mol_split.split(text)
        mol_blocks = [mb for mb in mol_blocks if mb.strip()]
        # --- Find single-molecule, biatomic non-matches ---
        if len(mol_blocks) == 1:
            atom_count = count_atoms(mol_blocks[0])
            if atom_count == 2:
                biatomic_nonmatches.append(sdf_path.name)
        # --- Find biatomic fallback blocks ---
        for mol_block in mol_blocks:
            if re_fallback_tag.search(mol_block):
                atom_count = count_atoms(mol_block)
                if atom_count == 2:
                    biatomic_fallbacks.append(sdf_path.name)
                    break  # Only report each file once for fallbacks

    print(f"Total .sdf files scanned:        {total_files}")
    print(f"Biatomic non-match files found:  {len(biatomic_nonmatches)}")
    print(f"SDF files with a biatomic Fallback: {len(biatomic_fallbacks)}\n")

    if biatomic_nonmatches:
        print("Files with a single, biatomic molecule (non-matches):")
        for fp in biatomic_nonmatches:
            print(f"  {fp}")
    else:
        print("No biatomic non-match files found.")

    print()  # extra line

    if biatomic_fallbacks:
        print("Files containing at least one distorted match (Fallback) for a biatomic CAE:")
        for fp in biatomic_fallbacks:
            print(f"  {fp}")
    else:
        print("No files with a biatomic Fallback found.")

if __name__ == "__main__":
    main()
