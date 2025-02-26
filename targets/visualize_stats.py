import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def read_stats_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    metals = defaultdict(int)
    geometries = defaultdict(int)
    populations = defaultdict(int)
    
    current_section = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        if line == "METALS":
            current_section = "metals"
            continue
        elif line == "GEOMETRIES":
            current_section = "geometries"
            continue
        elif line == "INITIAL POPULATIONS":
            current_section = "populations"
            continue
        
        if current_section == "metals" or current_section == "geometries":
            parts = line.split(": ")
            if len(parts) == 2:
                if current_section == "metals":
                    metals[parts[0]] = int(parts[1])
                else:
                    geometries[parts[0]] = int(parts[1])
        elif current_section == "populations":
            parts = line.split(": ")
            if len(parts) == 2:
                key = parts[0].split("-")[1]  # Extract name after number
                populations[key] = int(parts[1])
    
    return metals, geometries, populations

def plot_histogram(data, title):
    sorted_items = sorted(data.items(), key=lambda x: x[1], reverse=True)
    if title == "Geometries":
        labels, values = zip(*[(key.replace("_", "\u00A0"), val) for key, val in sorted_items])
    else:
        labels, values = zip(*sorted_items)
    
    max_value = max(values)
    min_value = min(values)
    avg_value = np.mean(values)
    
    normalized_values = [v / max_value * 5 for v in values]  # Scale down bar lengths
    
    plt.figure(figsize=(6, 6))
    if title == "Metals":
        bars = plt.barh(labels, normalized_values, color='dodgerblue')
    elif title == "Geometries":
        bars = plt.barh(labels, normalized_values, color='firebrick')
    else:
        bars = plt.barh(labels, normalized_values, color='forestgreen')
        
    
    for bar, value in zip(bars, values):
        plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/3, f"{value}", 
                 va='center', fontweight='bold')
    
    if title == "Initial Populations":
        plt.axvline(x=min_value / max_value * 5, color='green', linestyle='--', label=f'Min: {min_value}')
        plt.axvline(x=max_value / max_value * 5, color='red', linestyle='--', label=f'Max: {max_value}')
        plt.axvline(x=avg_value / max_value * 5, color='blue', linestyle='--', label=f'Avg: {avg_value:.0f}')
        
    plt.yticks(ticks=range(len(labels)), labels=[f"$\\it{{{label}}}$" for label in labels])
    plt.title(f"{title} Distribution", fontweight='bold')
    plt.gca().invert_yaxis()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.xticks([])
    if title == "Initial Populations":
        plt.legend(loc='lower right', frameon=True, facecolor='white', edgecolor='black', shadow=True, framealpha=1)
    plt.tight_layout()
    plt.show()

def main():
    filename = "stats.txt"
    metals, geometries, populations = read_stats_file(filename)
    
    plot_histogram(metals, "Metals")
    plot_histogram(geometries, "Geometries")
    plot_histogram(populations, "Initial Populations")

if __name__ == "__main__":
    main()
