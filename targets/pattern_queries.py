# Script defining all the queries for the pattern search in the CSD
      
import ccdc

# TODO: Check there are all the metals and non-metals
METALS = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                 'Sn', 'Sb', 'Te', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
                 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
                 'Lr']
NON_METALS = ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Se', 'Br', 'I', 'At']

### FOR METAL CONTAINING COMPOUNDS ###

### METAL IN RINGS QUERIES ###
# ONE METAL IN A RING, chelating ligands
def chelating_ligands(metals):
    query = ccdc.search.QuerySubstructure()
    # Add a metal query atom with a cyclic constraint
    metal = ccdc.search.QueryAtom(metals)
    metal.cyclic = True  
    # Add the QueryAtom to the QuerySubstructure
    query.add_atom(metal)
    return query

# MONO BRIDGED COMPLEXES
def mono_bridged_complexes(metals, non_metals, num_intermediates):
    # Iterate from 1 to num_intermediates
    queries = []
    for i in range(1, num_intermediates+1):
        query = ccdc.search.QuerySubstructure()
        # Add the metal atom
        metal1 = query.add_atom(ccdc.search.QueryAtom(metals))
        metal2 = query.add_atom(ccdc.search.QueryAtom(metals))
        # Add intermediate ligand atoms
        previous_atom = metal1
        for _ in range(i):
            ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
            ligand.cyclic = False
            query.add_bond('Any', previous_atom, ligand) 
            previous_atom = ligand  
        # Add the second metal atom
        query.add_bond('Any', previous_atom, metal2)  
        # Append the query to the list of queries
        queries.append(query)
    return queries

# TWO METALS IN A RING, multi-bridged complexes (Iterative search)
def multi_bridged_complexes(metals, ligands, num_intermediates):
    # iterate from 1 to num_intermediates
    queries = []
    for i in range(1, num_intermediates+1):
        # Create a new query substructure
        query = ccdc.search.QuerySubstructure()
        # Add the first metal atom
        metal1 = query.add_atom(ccdc.search.QueryAtom(metals))
        metal1.cyclic = True 
        # Add intermediate ligand atoms
        previous_atom = metal1
        for _ in range(i):
            ligand = query.add_atom(ccdc.search.QueryAtom(ligands))
            ligand.cyclic = True  
            query.add_bond('Any', previous_atom, ligand) 
            previous_atom = ligand  
        # Add the second metal atom
        metal2 = query.add_atom(ccdc.search.QueryAtom(metals))
        metal2.cyclic = True 
        query.add_bond('Any', previous_atom, metal2)  
        # Append the query to the list of queries
        queries.append(query)
    return queries

# LINEAR LIGANDS
def linear_ligands(metals, no_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    atom1 = query.add_atom(ccdc.search.QueryAtom(no_metals))
    atom1.cyclic = False
    # Add sigma bond between the metal and the ligand
    query.add_bond('Any', metal, atom1)
    return query

# METAL ALLENES
def metal_allenes(metals, non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    atom1 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2 = query.add_atom(ccdc.search.QueryAtom('C'))
    bond1 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    bond2 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    query.add_bond(bond1, metal, atom1)
    query.add_bond(bond2, atom1, atom2)
    return query

### DATIVE BONDS QUERIES ###
# DATIVE BONDS (General presence of a Pi bond between a metal and a non-metal)
def dative_bonds(metals, non_metals):
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
def dative_in_rings(metals, non_metals, num_intermediates):
    # Iterate from 1 to num_intermediates
    queries = []
    for i in range(1, num_intermediates+1):
        query = ccdc.search.QuerySubstructure()
        non_pi_bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
        # Add the metal atom
        metal = query.add_atom(ccdc.search.QueryAtom(metals))
        ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
        # The first bond between the metal and the ligand can be any bond (also a pi bond)
        query.add_bond('Any', metal, ligand)
        previous_atom = ligand
        # Add the ligand atom(s) #start with at least one ligand
        for _ in range(i):
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
        # Append the query to the list of queries
        queries.append(query)
    return queries

# DATIVE CONJUGATED CHAINS
def dative_in_conjugated_chains(metals, non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    atom1 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2.cyclic = False
    atom3 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom3.cyclic = False
    atom4 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom4.cyclic = False
    bond1 = ccdc.search.QueryBond('Pi')
    bond2 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    bond3 = ccdc.search.QueryBond('Single')
    bond4 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    query.add_bond(bond1, metal, atom1)
    query.add_bond(bond2, atom1, atom2)
    query.add_bond(bond3, atom2, atom3)
    query.add_bond(bond4, atom3, atom4)
    return query

# DATIVE CONJUGATED RINGS
def dative_in_conjugated_rings(metals, non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    atom1 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom1.cyclic = True
    atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2.cyclic = True
    atom3 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom3.cyclic = True
    atom4 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom4.cyclic = True
    bond1 = ccdc.search.QueryBond('Pi')
    bond2 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    bond3 = ccdc.search.QueryBond('Single')
    bond4 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    query.add_bond(bond1, metal, atom1)
    query.add_bond(bond2, atom1, atom2)
    query.add_bond(bond3, atom2, atom3)
    query.add_bond(bond4, atom3, atom4)
    return query

# CYCLIC MULTIHAPTO LIGANDS
def cyclic_multihapto_ligands(metals, non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom (connected with a pi bond)
    ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
    ligand.cyclic = True
    ligand2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    ligand2.cyclic = True
    # Add pi bond between the metal and the ligand
    bond = ccdc.search.QueryBond('Pi')
    bond2 = ccdc.search.QueryBond('Pi')
    query.add_bond(bond, metal, ligand)
    query.add_bond(bond2, metal, ligand2)
    # Ensure the cyclic is with a aromatic bond
    query.add_bond('Aromatic', ligand, ligand2)
    return query

# SIGMA BONDS (AGOSTIC INTERACTIONS) (both in rings and not in rings)
def sigma_bonds(metals):
    query = ccdc.search.QuerySubstructure()
    # Add the metal atom
    metal = query.add_atom(ccdc.search.QueryAtom(metals))
    # Add the ligand atom
    ligand = query.add_atom(ccdc.search.QueryAtom('H'))
    # Add sigma bond between the metal and the ligand
    bond = ccdc.search.QueryBond('Pi')
    query.add_bond(bond, metal, ligand)
    return query

# SIGMA BONDS (AGOSTIC INTERACTIONS) only in rings
def sigma_bonds_in_rings(metals, non_metals, num_intermediates):
    # Iterate from 1 to num_intermediates
    queries = []
    for i in range(1, num_intermediates+1):
        query = ccdc.search.QuerySubstructure()
        non_pi_bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
        # Add the metal atom
        metal = query.add_atom(ccdc.search.QueryAtom(metals))
        ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
        # The first bond between the metal and the ligand can be any bond (also a pi bond)
        query.add_bond('Any', metal, ligand)
        previous_atom = ligand
        # Add the ligand atom(s) #start with at least one ligand
        for _ in range(i):
            ligand = query.add_atom(ccdc.search.QueryAtom(non_metals))
            bond = ccdc.search.QueryBond(['Single', 'Double', 'Triple', 'Delocalized', 'Aromatic'])
            query.add_bond(bond, previous_atom, ligand)
            previous_atom = ligand 
        # insert the Hydrogen in the ring (sigma complex)
        hydrogen = query.add_atom(ccdc.search.QueryAtom('H')) 
        query.add_bond(non_pi_bond, previous_atom, hydrogen)
        # Add dative bond between the metal and the non-metal
        dative_bond = ccdc.search.QueryBond('Pi')
        query.add_bond(dative_bond , metal, hydrogen)
        queries.append(query)
    return queries


### FOR ANY COMPOUND ###

### RINGS QUERIES ###
# SATURATED RINGS
def saturated_rings(non_metals, num_intermediates):
    # Iterate from 1 to num_intermediates
    queries = []
    for i in range(2, num_intermediates+1):
        query = ccdc.search.QuerySubstructure()
        # Add the first atom
        atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom1.cyclic = True
        previous_atom = atom1
        for _ in range(i):
            atom = query.add_atom(ccdc.search.QueryAtom(non_metals))
            atom.cyclic = True
            query.add_bond('Single', previous_atom, atom)
            previous_atom = atom
        # Close the ring 
        query.add_bond('Single', previous_atom, atom1)
        queries.append(query)
    return queries

# UNSATURATED RINGS
def unsaturated_rings(non_metals, num_intermediates):
    # Iterate from 1 to num_intermediates
    queries = []
    for i in range(1, num_intermediates+1):
        query = ccdc.search.QuerySubstructure()
        unsat_bond = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
        # Add the first atom
        atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom1.cyclic = True
        atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom2.cyclic = True
        query.add_bond(unsat_bond, atom1, atom2)
        previous_atom = atom2
        for _ in range(i):
            atom = query.add_atom(ccdc.search.QueryAtom(non_metals))
            atom.cyclic = True
            query.add_bond('Any', previous_atom, atom)
            previous_atom = atom
        # Close the ring
        query.add_bond('Any', previous_atom, atom1)
        queries.append(query)
    return queries

# UNSATURATED CONJUGATED RINGS
def unsaturated_conjugated_rings(non_metals, num_intermediates):
    # Iterate from 1 to num_intermediates
    queries = []
    for i in range(1, num_intermediates+1):
        query = ccdc.search.QuerySubstructure()
        unsat_bond1 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
        unsat_bond2 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
        # Add the first atom
        atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom1.cyclic = True
        atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom2.cyclic = True
        atom3 = query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom3.cyclic = True
        atom4 = query.add_atom(ccdc.search.QueryAtom(non_metals))
        atom4.cyclic = True
        query.add_bond(unsat_bond1, atom1, atom2)
        query.add_bond('Single', atom2, atom3)
        query.add_bond(unsat_bond2, atom3, atom4)
        previous_atom = atom4
        for _ in range(i):
            atom = query.add_atom(ccdc.search.QueryAtom(non_metals))
            atom.cyclic = True
            query.add_bond('Any', previous_atom, atom)
            previous_atom = atom
        # Close the ring
        query.add_bond('Any', previous_atom, atom1)
        queries.append(query)
        return query

# AROMATIC RINGS
def aromatic_rings(non_metals):
    query = ccdc.search.QuerySubstructure()
    aromatic_bond = ccdc.search.QueryBond('Aromatic')
    # Add the first atom
    atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom1.cyclic = True
    atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2.cyclic = True
    query.add_bond(aromatic_bond, atom1, atom2)
    return query

### CHAINS QUERIES ###

# SATURATED CHAINS
def saturated_chains(non_metals):
    query = ccdc.search.QuerySubstructure()
    atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2.cyclic = False
    bond1 = ccdc.search.QueryBond(['Single'])
    query.add_bond(bond1, atom1, atom2)
    return query

# UNSATURATED CHAINS
def unsaturated_chains(non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the first atom
    atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2.cyclic = False
    atom3 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom3.cyclic = False
    bond = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    query.add_bond(bond, atom1, atom2)
    query.add_bond('Any', atom2, atom3)
    return query

# UNSATURATED CONJUGATED CHAINS
def unsaturated_conjugated_chains(non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the first atom
    atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2.cyclic = False
    atom3 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom3.cyclic = False
    atom4 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    bond1 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    bond3 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    query.add_bond(bond1, atom1, atom2)
    query.add_bond('Single', atom2, atom3)
    query.add_bond(bond3, atom3, atom4)
    return query

# ALLENES
def allenes(non_metals):
    query = ccdc.search.QuerySubstructure()
    # Add the first atom
    atom1= query.add_atom(ccdc.search.QueryAtom(non_metals))
    atom2 = query.add_atom(ccdc.search.QueryAtom('C'))
    atom3 = query.add_atom(ccdc.search.QueryAtom(non_metals))
    bond1 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    bond2 = ccdc.search.QueryBond(['Double', 'Triple', 'Quadruple'])
    query.add_bond(bond1, atom1, atom2)
    query.add_bond(bond2, atom2, atom3)
    return query


        
    