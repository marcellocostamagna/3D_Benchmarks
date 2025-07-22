
import pandas as pd

import nltk
import numpy as np

#########################
# smiles grammar
gram = """smiles -> chain
atom -> bracket_atom
atom -> aliphatic_organic
atom -> aromatic_organic
aliphatic_organic -> 'B'
aliphatic_organic -> 'C'
aliphatic_organic -> 'F'
aliphatic_organic -> 'H'
aliphatic_organic -> 'I'
aliphatic_organic -> 'N'
aliphatic_organic -> 'O'
aliphatic_organic -> 'P'
aliphatic_organic -> 'S'
aliphatic_organic -> 'Cl'
aliphatic_organic -> 'Br'
aliphatic_organic -> 'Si'
aliphatic_organic -> 'Se'

aliphatic_organic -> 'Fe'  
aliphatic_organic -> 'Cu' 
aliphatic_organic -> 'Ni'
aliphatic_organic -> 'Ru'

aromatic_organic -> 'b'
aromatic_organic -> 'c'
aromatic_organic -> 'n'
aromatic_organic -> 'o'
aromatic_organic -> 'p'
aromatic_organic -> 's'
aromatic_organic -> 'se'
bracket_atom -> '[' BAI ']'
BAI -> isotope symbol BAC
BAI -> symbol BAC
BAI -> isotope symbol
BAI -> symbol
BAC -> chiral BAH
BAC -> BAH
BAC -> chiral
BAH -> hcount BACH
BAH -> BACH
BAH -> hcount
BACH -> charge
symbol -> aliphatic_organic
symbol -> aromatic_organic
isotope -> DIGIT
isotope -> DIGIT DIGIT
isotope -> DIGIT DIGIT DIGIT
DIGIT -> '1'
DIGIT -> '2'
DIGIT -> '3'
DIGIT -> '4'
DIGIT -> '5'
DIGIT -> '6'
DIGIT -> '7'
DIGIT -> '8'
DIGIT -> '9'
DIGIT -> '0'
chiral -> '@'
chiral -> '@@'
hcount -> 'H'
hcount -> 'H' DIGIT
charge -> '-'
charge -> '-' DIGIT
charge -> '-' DIGIT DIGIT
charge -> '+'
charge -> '+' DIGIT
charge -> '+' DIGIT DIGIT
bond -> '-'
bond -> '='
bond -> '#'
bond -> '/'
bond -> '\\'
ringbond -> DIGIT
ringbond -> bond DIGIT
branched_atom -> atom
branched_atom -> atom RB
branched_atom -> atom BB
branched_atom -> atom RB BB
RB -> RB ringbond
RB -> ringbond
BB -> BB branch
BB -> branch
branch -> '(' chain ')'
branch -> '(' bond chain ')'
chain -> branched_atom
chain -> chain branched_atom
chain -> chain bond branched_atom
Nothing -> None"""

# form the CFG and get the start symbol
GCFG = nltk.CFG.fromstring(gram)

#############

def get_smiles_tokenizer(cfg):
    long_tokens = [a for a in cfg._lexical_index.keys() if len(a) > 1]
    # there are currently 6 double letter entities in the grammar (7 with a new metal)
    # these are their replacement, with no particular meaning
    # they need to be ascii and not part of the SMILES symbol vocabulary
    replacements = ['!', '?', '.', ',', ';', '$', '_', '>', '<', '§'] #(one symbol added)
    # print(f'length tokens: {len(long_tokens)}')
    # print(f'length replacements: {len(replacements)}')
    assert len(long_tokens) == len(replacements)
    for token in replacements:
        assert token not in cfg._lexical_index

    def tokenize(smiles):
        for i, token in enumerate(long_tokens):
            smiles = smiles.replace(token, replacements[i])
        tokens = []
        for token in smiles:
            try:
                ix = replacements.index(token)
                tokens.append(long_tokens[ix])
            except Exception:
                tokens.append(token)
        return tokens
    return tokenize


def encode(smiles):
    # GCFG = smiles_grammar.GCFG
    tokenize = get_smiles_tokenizer(GCFG)
    tokens = tokenize(smiles)
    parser = nltk.ChartParser(GCFG)
    try:
        parse_tree = parser.parse(tokens).__next__()
    except StopIteration:
        # print(f'Failed to parse {smiles}')
        return None
    productions_seq = parse_tree.productions()
    productions = GCFG.productions()
    prod_map = {}
    for ix, prod in enumerate(productions):
        prod_map[prod] = ix
    indices = np.array([prod_map[prod] for prod in productions_seq], dtype=int)
    return indices


def prods_to_eq(prods):
    seq = [prods[0].lhs()]
    for prod in prods:
        if str(prod.lhs()) == 'Nothing':
            break
        for ix, s in enumerate(seq):
            if s == prod.lhs():
                seq = seq[:ix] + list(prod.rhs()) + seq[ix + 1:]
                break
    try:
        return ''.join(seq)
    except Exception:
        return ''


def decode(rule):
    productions = GCFG.productions()
    prod_seq = [productions[i] for i in rule]
    return prods_to_eq(prod_seq)

#######################
successful_encodes = []
successful_smiles = []

df = pd.read_csv('./molecule_libraries/CSD/final_structures_with_smiles.csv')

smiles = df['SMILES'].tolist()
smiles = [smile for smile in smiles if type(smile) == str]
# smiles = [smile for smile in smiles if 'Fe' in smile or 'Cu' in smile or 'Ni' in smile or 'Ru' in smile]
smiles = [smile for smile in smiles if 'Ru' in smile]

count = 0
for smile in smiles:
    try:
        encoded = encode(smile)
    except Exception as e:
        # print(f"Failed to encode {smile}: {e}")
        continue
    
    if encoded is not None:
        successful_encodes.append((smile, encoded))
        print(f"Successful encode for {smile}: {encoded}")
        successful_smiles.append(smile)
        count += 1
        print(f"Count: {count}")
        
    if len(successful_encodes) >= 100:
        break
        
with open('./molecule_libraries/CSD/Tests_for_Benchmarks/Find_encodable_smiles/successful_smiles.txt', 'w') as file:
    for smile in successful_smiles:
        file.write(f"{smile}\n")