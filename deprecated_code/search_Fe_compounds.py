import pandas as pd
from ccdc import io
from ccdc.molecule import Molecule

# Load CSV
df = pd.read_csv("targets/viable_structures.csv")

# Set up CSD Reader
csd_reader = io.EntryReader("CSD")

# Lists to store identifiers
fe2_plus = []
fe3_plus = []

# Function to check for Fe2+ or Fe3+
def check_fe_oxidation_states(identifier):
    try:
        entry = csd_reader.entry(identifier)
        mol = entry.molecule

        for atom in mol.atoms:
            if atom.atomic_symbol == "Fe":
                ox = atom.formal_charge
                if ox == 2:
                    fe2_plus.append(identifier)
                elif ox == 3:
                    fe3_plus.append(identifier)
    except Exception as e:
        print(f"Failed to process {identifier}: {e}")

# Apply to each identifier
for identifier in df["Identifier"]:
    check_fe_oxidation_states(identifier)

# Print or save results
print("Fe2+ compounds:", fe2_plus)
print("Fe3+ compounds:", fe3_plus)

# Optionally save to file
with open("fe2_plus.txt", "w") as f2:
    f2.write("\n".join(fe2_plus))

with open("fe3_plus.txt", "w") as f3:
    f3.write("\n".join(fe3_plus))
