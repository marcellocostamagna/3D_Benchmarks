from ccdc import io
import os

input_txt_file = './CSD_similarity/similarity_results/similarity_results_grubbs_4d.txt'

with open(input_txt_file, 'r') as file:
    entries = {}
    lines_limit = 20
    for line in file:
        # limit the number of lines read
        if lines_limit == 0:
            break
        lines_limit -= 1
        # divide the line with the : character
        line = line.strip().split(':')
        
        key = line[0].strip()
        scores = line[1].strip()
        scores = scores.strip('[]')

        # extract the similarity scores and add the highest one to the dictionary
        sim_scores = [float(x) for x in scores.split(',')]
        entries[key] = max(sim_scores)
        
# cretate a directory for storing molecules files (sdf format)
output_dir = './CSD_similarity/molecules_for_benchmark_test'
os.makedirs(output_dir, exist_ok=True)

for entry in entries:
    # get the molecule with the highest similarity score
    csd_reader = io.EntryReader('CSD')
    molecule = csd_reader.entry(entry).molecule.heaviest_component
    # save the molecule to a file
    with io.MoleculeWriter(f'{output_dir}/{entry}.sdf') as w:
        w.write(molecule)
