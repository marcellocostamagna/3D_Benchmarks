        
import ccdc

# TODO: Check there are all the metals and non-metals
METALS = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                 'Sn', 'Sb', 'Te', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
                 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
                 'Lr']
NON_METALS = ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Se', 'Br', 'I', 'At']

# METAL IN RING QUERY
def metal_in_ring(metals):
    # Create a new query substructure
    metal_in_ring_query = ccdc.search.QuerySubstructure()
    # Add a metal query atom with a cyclic constraint
    metal = ccdc.search.QueryAtom(metals)
    metal.cyclic = True  # Add cyclic constraint
    # Add the QueryAtom to the QuerySubstructure
    metal_in_ring_query.add_atom(metal)
    return metal_in_ring_query


# # TWO METALS IN A RING, bridged complexes (Iterative search)
def two_metals_in_a_ring(metals, ligands, num_intermediates):
    # Create a new query substructure
    query = ccdc.search.QuerySubstructure()
    # Add the first metal atom
    metal1 = query.add_atom(ccdc.search.QueryAtom(metals))
    metal1.cyclic = True 
    # Add intermediate ligand atoms
    previous_atom = metal1
    for _ in range(num_intermediates):
        ligand = query.add_atom(ccdc.search.QueryAtom(ligands))
        ligand.cyclic = True  
        query.add_bond('Any', previous_atom, ligand) 
        previous_atom = ligand  
    # Add the second metal atom
    metal2 = query.add_atom(ccdc.search.QueryAtom(metals))
    metal2.cyclic = True 
    query.add_bond('Any', previous_atom, metal2)  
    return query

# DATIVE BONDS (both in rings and not in rings)
def dative_bond(metals, non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
    # Add pi bond between the metal and the ligand
    bond = ccdc.search.QueryBond('Pi')
    query.add_bond(bond, metal, ligand)
    return query

# DATIVE BONDS only in rings (Iterative search)
def dative_in_ring(metals, non_metals, num_intermediates):
    query = ccdc.search.QuerySubstructure()
    non_pi_bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    previous_atom = metal
    # Add the ligand atom(s) #start with at least one ligand
    for _ in range(num_intermediates):
        ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
        bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
        query.add_bond(bond, previous_atom, ligand)
        previous_atom = ligand 
    # insert the element in the ring
    non_metal = query.add_atom(ccdc.search.QueryAtom(non_metals)) 
    query.add_bond(non_pi_bond, previous_atom, non_metal)
    # Add dative bond between the metal and the non-metal
    dative_bond = ccdc.search.QueryBond('Pi')
    query.add_bond(dative_bond , metal, non_metal)
    return query

# SIGMA BONDS (AGOSTIC INTERACTIONS) (both in rings and not in rings)
def sigma_bond(metals, non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
    # Add sigma bond between the metal and the ligand
    bond = ccdc.search.QueryBond('Sigma')
    query.add_bond(bond, metal, ligand)
    return query

# SIGMA BONDS (AGOSTIC INTERACTIONS) only in rings
def sigma_bonds_in_ring(metals, non_metals, num_intermediates):
    query = ccdc.search.QuerySubstructure()
    non_pi_bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
    
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    previous_atom = metal
    # Add the ligand atom(s) #start with at least one ligand
    for _ in range(num_intermediates):
        ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
        bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
        query.add_bond(bond, previous_atom, ligand)
        previous_atom = ligand 
    # insert the element in the ring
    hydrogen = query.add_atom(ccdc.search.QueryAtom('H')) 
    query.add_bond(non_pi_bond, previous_atom, hydrogen)
    # Add dative bond between the metal and the non-metal
    dative_bond = ccdc.search.QueryBond('Pi')
    query.add_bond(dative_bond , metal, hydrogen)
    return query
