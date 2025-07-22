import os

def extract_unique_entries_from_populations(folder_path):
    unique_entries = set()

    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            print(filename)
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r', encoding='utf-8') as file:
                for line in file:
                    if line.strip():  # skip empty lines
                        entry = line.strip().split()[0]  # get first field
                        unique_entries.add(entry)

    print(f"✅ Total unique entries across all population files: {len(unique_entries)}")
    return unique_entries

folder = "targets/init_populations_protons/"
unique_mols = extract_unique_entries_from_populations(folder)

