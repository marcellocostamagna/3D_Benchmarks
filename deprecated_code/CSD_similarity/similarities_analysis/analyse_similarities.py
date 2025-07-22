import pandas as pd
import os

# Directory containing the CSV files
results_dir = 'CSD_similarity/similarities_analysis/'

# Targets and similarity measures
targets = ['grubbs', 'polymerization', 'PNP', 'salen']
similarities = ['3d', '4d', 'csd', 'pybel']

# Function to load DataFrames from CSV files and filter valid rows
def load_and_filter_dataframes(results_dir, targets):
    dataframes = {}
    for target in targets:
        file_path = os.path.join(results_dir, f'combined_similarity_results_{target}.csv')
        df = pd.read_csv(file_path)
        # Filter out rows with any missing or non-numeric values
        df = df.dropna()
        df = df[df[similarities].applymap(lambda x: isinstance(x, (int, float))).all(axis=1)]
        dataframes[target] = df
    return dataframes

# Load the DataFrames
dataframes = load_and_filter_dataframes(results_dir, targets)

# Perform Descriptive Statistics
def descriptive_statistics(df):
    return df.describe()

# Perform Correlation Analysis
def correlation_analysis(df):
    return df.corr()

# Perform Distribution Analysis (Histogram)
def plot_distributions(df, target):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle(f'Distribution of Similarity Measures for {target}')
    
    for i, col in enumerate(similarities):
        ax = axes[i // 2, i % 2]
        ax.hist(df[col], bins=15)
        ax.set_xlabel('Similarity Measure Value')
        ax.set_ylabel('Frequency')
        ax.set_title(f'Distribution of {col}')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

# Perform Outlier Detection
def detect_outliers(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    outliers = df[((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))).any(axis=1)]
    return outliers

# Analysis for each target
analysis_results = {}
for target, df in dataframes.items():
    print(f"\nAnalysis for {target}")
    
    # Descriptive Statistics
    print("Descriptive Statistics:")
    descriptive_stats = descriptive_statistics(df)
    print(descriptive_stats)
    
    # Correlation Analysis
    print("\nCorrelation Analysis:")
    correlation_matrix = correlation_analysis(df)
    print(correlation_matrix)
    
    # Distribution Analysis
    print("\nDistribution Analysis:")
    plot_distributions(df, target)
    
    # Outlier Detection
    print("\nOutliers:")
    outliers = detect_outliers(df)
    print(outliers)
    
    # Save analysis results
    analysis_results[target] = {
        'descriptive_statistics': descriptive_stats,
        'correlation_analysis': correlation_matrix,
        'outliers': outliers
    }

# # Ensure 'openpyxl' is installed
# try:
#     import openpyxl
# except ImportError:
#     import subprocess
#     subprocess.check_call([sys.executable, "-m", "pip", "install", "openpyxl"])
#     import openpyxl

# Save analysis results to Excel
# with pd.ExcelWriter('analysis_results.xlsx', engine='openpyxl') as writer:
#     for target, results in analysis_results.items():
#         results['descriptive_statistics'].to_excel(writer, sheet_name=f'{target}_desc_stats')
#         results['correlation_analysis'].to_excel(writer, sheet_name=f'{target}_corr_analysis')
#         results['outliers'].to_excel(writer, sheet_name=f'{target}_outliers')
