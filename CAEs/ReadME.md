# Connected Atom Environments (CAEs) 

This folder contains scripts to enumerate and compare Connected Atom Environments (CAEs) —substructures defined by a central atom and its directly bonded neighbors, with all atomic features and 3D geometry- between target molecules and their candidate starting populations.

The goal is to quantify whether a starting population contains all the molecular fragments (CAEs) needed to reconstruct a given target.

## Scripts

- `get_CAEs_from_viable_structures.py`:  
Builds a DuckDB database(`all_caes.duckdb`) of all CAEs from `viable_structures.csv`.  
*Note: Processing nearly a million structures may take several hours.*
- `CAE_comparison.py`:  
Compares the CAEs of each target to those present in its starting population. Outputs a log file and SDF files for all CAEs compared for downstream analysis.

**Note:** Scripts for summarizing and visualizing CAE matching results can be found in the `CAEs_analysis` folder.

## Usage
Run scripts directly as:  
```bash
python script.py
```
Replace `script.py` with the script you wish to run.
Ensure required input files (`viable_structures.csv`, `all_caes.duckdb`) are present.

**Tip:**  
To compare different starting populations, edit the `population_folder` variable at the top of `CAE_comparison.py` to point to the desired starting population folder, e.g.:
```bash
population_folder = "../Starting_populations/Starting_populations_0_5"  # <--- CHANGE THIS IF NEEDED
```





