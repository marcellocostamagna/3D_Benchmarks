# Script to perform the cretaion of a initial population 
# form the choice of targets and a similarity cutoff and 
# check if the initial popultaion is viable; meaning all the
# representatitve molecular patterns present in the targets 
# are also present in the initial population

# PIPELINE:

# 1. Definition of patterns of interest

# 2. Targets, read from a file (only manual input)

# 3. Viable molecules, read from a file

# 4. Initial population (potential) : filtering of viable molecules
#    - Only molecules with a similarity to the targets below a cutoff are kept

# 5. Check if the initial population is viable
#    - All the representative molecular patterns are at least present in one molecule of the initial population
#    - (Optional) Check the frequency of the patterns in the initial population

# 6. Alternatives:
#   If the initial population is not viable:
#   - Highlight the patterns (and the associated targets) that are not present in the initial population
#   - Demand a new choice of targets
#   
#   If the initial population is viable:
#   - Return the initial population
#   - (Optional) Return the frequency of the patterns in the initial population

