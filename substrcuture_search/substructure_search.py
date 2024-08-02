# Substructure searching with 3D measurements
from ccdc.search import SubstructureSearch, MoleculeSubstructure, SMARTSSubstructure
import pandas as pd

# Strting dataframe 
df = pd.read_csv('CSD_similarity/final_structures_with_fp.csv')
df_identifiers = df['Identifier'].tolist()

### SMARTS ###

## GRUBBS ##

# 1. Ru with 5 ligands, 2 Carbons, 1 Nitrogen, 2 Chlorines 
# smarts_query = SMARTSSubstructure("[RuD5](~[C])(~[C])(~[N])(~[Cl])(~[Cl])")

# 2. Carbene
# Persistent carbene
# smarts_query = SMARTSSubstructure("[#6]-[#6]1:[#6]:[#6](-[#6]):[#6](:[#6](-[#6]):[#6]:1)-[#7]1-[#6]-[#6]-[#7](-[#6]-1~[*;!#1])-[#6]1:[#6](-[#6]):[#6]:[#6](-[#6]):[#6]:[#6]:1-[#6]")
# General NHC
# smarts_query = SMARTSSubstructure("[#7;D3]1(~[*;!#1])-[#6]-[#6]-[#7;D3](~[*;!#1])-[#6]-1-[#6]")


## PNP ##
# 1. Unbonded PNP
# smarts_query = SMARTSSubstructure("[#6]1:[#6]:[#6]:[#6]2:[#6]:[#6]3:[#6]:[#6]:[#6]:[#6](-[#6]-[#15](~[*])(~[*])~[*]):[#6]:3:[#7]:[#6]:2:[#6]:1-[#6]-[#15](~[*])(~[*])~[*]")

# 2. P bounded PNP
# smarts_query = SMARTSSubstructure('[#6]1:[#6]:[#6]:[#6]2:[#6]:[#6]3:[#6]:[#6]:[#6]:[#6]4:[#6]:3:[#7]:[#6]:2:[#6]:1-[#6]-[#15]([*;!#1])([*;!#1])-[*;!#1]-[#15]([*;!#1])([*;!#1])-[#6]-4')

# 3. PNP skeleton   
# smarts_query = SMARTSSubstructure('[#6]1~[#6]~[#7]~[#6]~[#6]-[#6]-[#15]-[*;!#1]-[#15]-[#6]-1')


## POLYMERIZATION ##
# 1. Zr with two Cl ligands and 12 bonds
# smarts_query = SMARTSSubstructure('[ZrD12](~[#9,#17,#35,#53])~[#9,#17,#35,#53]')

# 2. Zr with two Cp ligands and two Cl ligands
# smarts_query = SMARTSSubstructure('[Zr]1234(~[#6]~5~[#6]~1~[#6]~2~[#6]~3~[#6]~45)5678(~[#6]~9~[#6]~5~[#6]~6~[#6]~7~[#6]~89)(~[#9,#17,#35,#53])~[#9,#17,#35,#53]')

# 3. Ligannd with Silicon
# smarts_query = SMARTSSubstructure('[#6]1~[#6]~[#6]~[#6]~[#6]~1~[Si]~[Si]~[#6]1~[#6]~[#6]~[#6]~[#6]~1') 


## SALEN TYPE ##
# 1. Cu with two O and two N ligands
# smarts_query = SMARTSSubstructure('[CuD4](~[#8])(~[#8])(~[#7])~[#7]')

# 2. Salen specific ligand
# smarts_query = SMARTSSubstructure('[#6]1~[#6][#6]~[#6]~[#6]2~[#6]~1~[#6]1~[#6](~[#8]~[*;!#1]34~[#7](~[#6]~[#6]~[#6]~[#7]~3~[#6](~[#6])~[#6]3~[#6](~[#6]~[#6]~[#6]~[#6]~3)~[#8]~4)~[#6]~1)~[#6]~[#6]~2')

# 3. General Salen ligand
smarts_query = SMARTSSubstructure('[#8]1~[Cu]23~[#7](~[#6]~[#6]~[#6]~1)~[#6]~[#6]~[#7]~2~[#6]~[#6]~[#6]~[#8]~3')

# smarts_query = SMARTSSubstructure("[Ru]")
substructure_search = SubstructureSearch()
sub_id = substructure_search.add_substructure(smarts_query)
csd_hits = substructure_search.search()
print(f'Number of hits: {len(csd_hits)}')

# Get the hits that are in the dataframe
hits = [hit.identifier for hit in csd_hits if hit.identifier in df_identifiers]
print(hits)
print(f'Number of hits: {len(hits)}')

#  There is a minor extension to Daylight SMARTS to allow the representation of
#  quadruple, delocalised and pi bonds, using the characters '_', '"' and '|' respectively.