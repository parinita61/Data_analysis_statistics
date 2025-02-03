import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the dataset
data = pd.read_csv("PupilBioTest_PMP_revA.csv")

# Filter data by tissue types
tissue1 = data[data['Tissue'] == 'cfDNA']  
tissue2 = data[data['Tissue'] == 'Islet']  

# Sum the CpG columns to get the coverage for each row
cpg_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
data['Coverage'] = data[cpg_columns].sum(axis=1)

# Calculate median and CV for single CpG coverage in each tissue
def calculate_coverage_stats(tissue_data):
    median_coverage = tissue_data['Coverage'].median()
    cv_coverage = tissue_data['Coverage'].std() / tissue_data['Coverage'].mean()
    return median_coverage, cv_coverage

tissue1_median, tissue1_cv = calculate_coverage_stats(tissue1)
tissue2_median, tissue2_cv = calculate_coverage_stats(tissue2)

# Save results
coverage_stats = pd.DataFrame({
    "Tissue": ["cfDNA", "Islet"],  # Use actual tissue names
    "Median Coverage": [tissue1_median, tissue2_median],
    "CV": [tissue1_cv, tissue2_cv]
})

# Plot median and CV
plt.figure(figsize=(10, 5))

# Median Coverage Plot
plt.subplot(1, 2, 1)
plt.bar(coverage_stats["Tissue"], coverage_stats["Median Coverage"], color=['blue', 'green'])
plt.title("Median Coverage by Tissue")
plt.ylabel("Median Coverage")

# Coefficient of Variation Plot
plt.subplot(1, 2, 2)
plt.bar(coverage_stats["Tissue"], coverage_stats["CV"], color=['blue', 'green'])
plt.title("Coefficient of Variation (CV)")
plt.ylabel("CV")

plt.tight_layout()
plt.savefig("coverage_stats_plot.png")
plt.show()
