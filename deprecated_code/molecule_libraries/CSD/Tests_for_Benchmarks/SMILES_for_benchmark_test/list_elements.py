import re


def extract_elements(smiles_list):
    # Define a regex pattern to match element symbols
    # Elements are represented by either a single uppercase letter or an uppercase letter followed by a lowercase letter
    pattern = r'\[[A-Z][a-z]?\]|\b[A-Z][a-z]?\b'

    elements = set()

    for smiles in smiles_list:
        # Find all elements using regex
        matches = re.findall(pattern, smiles)
        
        # Process each match, stripping brackets if necessary
        for match in matches:
            # If the match has brackets like [Fe], strip them
            element = match.strip('[]')
            elements.add(element)

    return elements

def main():
    smiles_file = './molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_for_benchmark_test/smiles.txt'

    # Read SMILES from file and store them in a list
    with open(smiles_file, 'r') as f:
        smiles_list = f.readlines()
        smiles_list = [smiles.strip() for smiles in smiles_list]

    # Extract elements from the SMILES list
    unique_elements = extract_elements(smiles_list)

    # Print out the unique elements
    # print(f'Unique elements found in the SMILES list:')
    # for element in sorted(unique_elements):
    #     print(element)
    
    # Save the unique elements to a file
    with open('./molecule_libraries/CSD/Tests_for_Benchmarks/SMILES_for_benchmark_test/elements.txt', 'w') as f:
        for element in sorted(unique_elements):
            f.write(f'{element}\n')

if __name__ == '__main__':
    main()