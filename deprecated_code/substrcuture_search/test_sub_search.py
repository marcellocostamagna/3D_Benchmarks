
# Substructure searching with 3D measurements
from ccdc.search import SubstructureSearch, MoleculeSubstructure, SMARTSSubstructure
import pandas as pd

# Strting dataframe 
df = pd.read_csv('CSD_similarity/final_structures_with_fp.csv')
df_identifiers = df['Identifier'].tolist()


## POLYMERIZATION ##
# 1. Zr with two Cp ligands and two Cl ligands
smarts_query_1 = SMARTSSubstructure('[ZrD12](~[#9,#17,#35,#53])~[#9,#17,#35,#53]')

smarts_query_2 = SMARTSSubstructure('[Zr]1234(~[#6]~5~[#6]~1~[#6]~2~[#6]~3~[#6]~45)5678(~[#6]~9~[#6]~5~[#6]~6~[#6]~7~[#6]~89)(-[#17])-[#17]')

# 2. Ligannd with Silicon
# smarts_query = SMARTSSubstructure('C1C=CC=C1[Si][Si]C1C=CC=C1')


# smarts_query = SMARTSSubstructure("[Ru]")
substructure_search_1 = SubstructureSearch()
substructure_search_2 = SubstructureSearch()
sub_id_1 = substructure_search_1.add_substructure(smarts_query_1)
sub_id_2 = substructure_search_2.add_substructure(smarts_query_2)
csd_hits_1 = substructure_search_1.search()
csd_hits_2 = substructure_search_2.search()
print(f'Number of hits in first search: {len(csd_hits_1)}')
print(f'Number of hits in second search: {len(csd_hits_2)}')

# Get the hits that are in the dataframe
hits_1 = [hit.identifier for hit in csd_hits_1 if hit.identifier in df_identifiers]
hits_2 = [hit.identifier for hit in csd_hits_2 if hit.identifier in df_identifiers]

# Find the non common hits and specify in which search they were found
non_common_hits = [hit for hit in hits_1 if hit not in hits_2]
print(non_common_hits)
print(f'Number of hits in first search: {len(hits_1)}')
print(f'Number of hits in second search: {len(hits_2)}')



#  There is a minor extension to Daylight SMARTS to allow the representation of
#  quadruple, delocalised and pi bonds, using the characters '_', '"' and '|' respectively.