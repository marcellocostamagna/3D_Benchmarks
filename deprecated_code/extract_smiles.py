# Define input file
file_path = "initial_population_test_ABAHIW.txt"

# Read the file and process each line
with open(file_path, "r") as file:
    lines = [line.split(":", 1)[1].strip() + "\n" for line in file if ":" in line]

output_file = "population_6.txt"

# Write the modified content back to the file
with open(output_file, "w") as file:
    file.writelines(lines)

print("File processed successfully!")
