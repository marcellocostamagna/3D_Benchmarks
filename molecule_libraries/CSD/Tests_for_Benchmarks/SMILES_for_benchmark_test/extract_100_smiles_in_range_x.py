import pandas as pd
import numpy as np

# text smiles with the entry and the similarity score
sim_text = '/Users/marcellocostamagna/3D_Benchmarks/CSD_similarity/similarity_results/similarity_results_grubbs_4d.txt'

# extract the names before ':' 
identifiers = []
scores = []

with open(sim_text, 'r') as f:
    for line in f:
        # Extract the identifier (name before ':')
        identifier = line.split(':')[0].strip()
        identifiers.append(identifier)
        
        # Get the list after the ':' in every line and parse it as a list of floats
        tmp_scores = line.split(':')[1].strip().strip('[]').split(',')
        if tmp_scores[0] == '':
            scores.append(0.0)
        else:
            tmp_scores = [float(score.strip()) for score in tmp_scores]
            # Extract the highest score from the list of scores
            highest_score = max(tmp_scores)
            scores.append(highest_score)

# Combine identifiers and scores into a list of tuples and sort by score
data = list(zip(identifiers, scores))
data.sort(key=lambda x: x[1])

# Select identifiers around the given score
selected_score = 0.3  # Example score
num_identifiers = 100  # Number of unique identifiers to gather
half_window = num_identifiers // 2

# Find the closest score index
distance = [abs(score - selected_score) for _, score in data]
closest_index = distance.index(min(distance))

# Gather identifiers around the selected score
start_index = max(0, closest_index - half_window)
end_index = min(len(data), closest_index + half_window)

# Adjust the start and end indices to ensure we gather the correct number of identifiers
while (end_index - start_index) < num_identifiers:
    if start_index > 0:
        start_index -= 1
    elif end_index < len(data):
        end_index += 1
    else:
        break

selected_identifiers = [identifier for identifier, _ in data[start_index:end_index]]

# csv file with the smiles
csv_file = '/Users/marcellocostamagna/3D_Benchmarks/molecule_libraries/CSD/final_structures_with_smiles.csv'

# Extract the smiles from the CSV file for the selected identifiers
df = pd.read_csv(csv_file)
smiles = []
identifier_score_dict = {}

for identifier in selected_identifiers:
    smile = df[df['Identifier'] == identifier]['SMILES'].values
    if len(smile) > 0 and isinstance(smile[0], str) and smile[0].strip() and not pd.isna(smile[0]):
        smiles.append(smile[0].strip())
        identifier_score_dict[smile[0].strip()] = dict(data)[identifier]

# Ensure uniqueness and gather until we have the desired number of unique smiles
unique_smiles = list(set(smiles))
extra_index = end_index

while len(unique_smiles) < num_identifiers:
    if extra_index < len(data):
        identifier, _ = data[extra_index]
        smile = df[df['Identifier'] == identifier]['SMILES'].values
        if len(smile) > 0 and isinstance(smile[0], str) and smile[0].strip() and not pd.isna(smile[0]) and smile[0] not in unique_smiles:
            unique_smiles.append(smile[0].strip())
            identifier_score_dict[smile[0].strip()] = dict(data)[identifier]
        extra_index += 1
    elif start_index > 0:
        start_index -= 1
        identifier, _ = data[start_index]
        smile = df[df['Identifier'] == identifier]['SMILES'].values
        if len(smile) > 0 and isinstance(smile[0], str) and smile[0].strip() and not pd.isna(smile[0]) and smile[0] not in unique_smiles:
            unique_smiles.append(smile[0].strip())
            identifier_score_dict[smile[0].strip()] = dict(data)[identifier]
    else:
        break

# Write the unique smiles to a file with the score in the filename
output_file = f'/Users/marcellocostamagna/3D_Benchmarks/molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_for_benchmark_test/smiles_for_grubbs_{selected_score:.2f}.txt'
with open(output_file, 'w') as f:
    for smile in unique_smiles:
        f.write(f'{smile}\n')

# Create a dictionary with smiles and their corresponding scores, sorted from highest to lowest
sorted_identifier_score_dict = dict(sorted(identifier_score_dict.items(), key=lambda item: item[1], reverse=True))

score_output_file = f'/Users/marcellocostamagna/3D_Benchmarks/molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_for_benchmark_test/smiles_with_scores_for_grubbs_{selected_score:.2f}.txt'
with open(score_output_file, 'w') as f:
    for smile, score in sorted_identifier_score_dict.items():
        f.write(f'{smile}: {score}\n')
