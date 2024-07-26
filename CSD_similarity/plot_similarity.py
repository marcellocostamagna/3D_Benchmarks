import matplotlib.pyplot as plt
import re

def read_data(file_path):
    """Read the file and extract maximum values for each entry."""
    max_values = []
    with open(file_path, 'r') as file:
        for line in file:
            try:
                # Extract everything within the brackets
                numbers = re.search(r'\[(.*?)\]', line).group(1)
                if numbers:  # Check if the list is not empty
                    # Convert numbers to float and find the maximum
                    max_value = max(map(float, numbers.split(',')))
                    max_values.append(max_value)
            except AttributeError:
                # Skip lines where no numbers are found
                continue
    return max_values

def plot_max_values(max_values):
    """Plot the maximum values."""
    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(max_values, marker='o', linestyle='-', color='b')
    plt.title('Maximum Values per Entry')
    plt.xlabel('Entries')
    plt.ylabel('Max Value')
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    # file_path = 'similarity_results_grubbs.txt'  
    # file_path = 'similarity_results_PNP.txt'
    # file_path = 'similarity_results_salen.txt'
    file_path = 'similarity_results_polymerization.txt'
    max_values = read_data(file_path)
    plot_max_values(max_values)
