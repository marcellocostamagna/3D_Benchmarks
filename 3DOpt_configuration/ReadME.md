# 3DOpt_configuration

This folder contains scripts to generate the necessary configuration files for building the 3DOpt benchmark inside MolScore. It also includes a script to generate 2D diagrams of the final target molecules chosen for the benchmark.

## Scripts

- `generate_json_tasks_files.py`:  
Generates `.json` configuration files defining each taks associated with a target, which together compose the 3DOpt benchmark.
- `generate_target_files.py`:  
Generates the `.txt`files defining the reference (taget) molecules for each task. Each file contains the CSD entry code, the SMILES, and the appropriate HSR fingerprint of the chosen target.
- `generate_targets_diagrams.py`:  
For each target, extracts the relevant component and creates a 2D diagram (SVG) using CCDC tools for visualization and reporting.