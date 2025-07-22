# Metal targets statistics
This folder contains summary data and visualization scripts for the selected 3DOpt benchmark targets that feature transition metals. The goal is to illustrate the diversity of both metal centers and their coordination geometries within the chosen set.

## Files

- `Metal_targets_stats.txt`:  
Tabulates the occurrence of each metal and its corresponding coordination geometry among the targets.

## Script

- `visualize_metal_targets_stats.py`:  
Reads `Metal_targets_stats.txt` and generates two horizontal bar charts:

    - One for the frequency of each transition metal present

    - One for the distribution of coordination geometries
    Plots are saved as `stats_metals.svg` and `stats_geometries.svg`.

## Usage
Simply run the command:
```bash
python visualize_metal_targets_stats.py
```
This will generate SVG visualizations summarizing the statistics.
