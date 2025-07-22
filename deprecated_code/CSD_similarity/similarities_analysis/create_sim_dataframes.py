import pandas as pd
import os
import re

# Function to parse the text file and process the data
def parse_and_process_file(file_path):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            match = re.match(r"(\w+): \[(.*)\]", line.strip())
            if match:
                identifier = match.group(1)
                values_str = match.group(2).split(',')
                # Filter out empty strings before converting to float
                values = [float(v) for v in values_str if v.strip()]
                if values and any(v != 0 for v in values):  # Exclude empty or zero-only lists
                    max_value = max(values)
                    data[identifier] = max_value
    return data

if __name__ == '__main__':
    results_dir = 'CSD_similarity/similarity_results'
    targets = ['grubbs', 'polymerization', 'PNP', 'salen']
    similarities = ['3d', '4d', 'csd', 'pybel']

    target_data = {}
    for target in targets:
        combined_data = {}
        identifier_files = {similarity: f'{results_dir}/similarity_results_{target}_{similarity}.txt' for similarity in similarities}
        
        for similarity, file in identifier_files.items():
            data = parse_and_process_file(file)
            for identifier, value in data.items():
                if identifier not in combined_data:
                    combined_data[identifier] = {}
                combined_data[identifier][similarity] = value

        target_data[target] = combined_data

    # Collect results and identifiers present in all four files
    shared_identifiers = set(target_data[targets[0]].keys())
    for target in targets:
        shared_identifiers.intersection_update(target_data[target].keys())

    results = {}
    non_shared_identifiers = {}

    for target in targets:
        target_results = []
        target_non_shared = {}
        for identifier, values in target_data[target].items():
            if identifier in shared_identifiers:
                target_results.append([identifier] + [values.get(similarity, None) for similarity in similarities])
            else:
                target_non_shared[identifier] = [similarity for similarity in similarities if similarity in values]
        results[target] = pd.DataFrame(target_results, columns=['Identifier'] + similarities)
        non_shared_identifiers[target] = target_non_shared

    # Write non-shared identifiers to text files for each target
    for target, identifiers in non_shared_identifiers.items():
        with open(f'non_shared_identifiers_{target}.txt', 'w') as f:
            for identifier, sims in identifiers.items():
                f.write(f"{identifier}: {', '.join(sims)}\n")

    # Save each target's DataFrame to a separate CSV file
    for target, df in results.items():
        df.to_csv(f'combined_similarity_results_{target}.csv', index=False)

