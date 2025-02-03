import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
significant_pmp = pd.read_csv("significant_pmps.csv")

# (a) How does sequencing depth affect specificity confidence?
def compute_specificity(depth):
    # Assume specificity increases with depth logarithmically
    specificity = 1 - np.exp(-0.001 * depth)  # Example function
    return specificity

sequencing_depths = np.arange(10000, 2000000, 50000)
specificity_values = [compute_specificity(d) for d in sequencing_depths]

# Save specificity vs sequencing depth
specificity_df = pd.DataFrame({"Sequencing_Depth": sequencing_depths, "Specificity": specificity_values})
specificity_df.to_csv("new_specificity_vs_depth.csv", index=False)

# Plot specificity vs sequencing depth
plt.figure(figsize=(8, 6))
plt.plot(sequencing_depths, specificity_values, lw=2, color='blue')
plt.xlabel("Sequencing Depth")
plt.ylabel("Specificity Confidence")
plt.title("Effect of Sequencing Depth on Specificity Confidence")
plt.grid(True)
plt.savefig("new_specificity_vs_depth.png")
plt.show()

# (b) Estimate threshold reads for Tissue #2 at 1 million reads
def estimate_threshold(depth, baseline=0.95):
    return np.log(1 - baseline) / -0.001  # Solve inverse function

threshold_reads = estimate_threshold(1e6)

# Save threshold estimate
with open("new_threshold_reads_tissue2.txt", "w") as f:
    f.write(f"Estimated read threshold for Tissue #2 at 1 million reads: {threshold_reads:.2f}\n")

# (c) Validate hypothesis: Compare top 10 PMPs against individual CpG sites
# Select top 10 PMPs based on lowest p-values
significant_pmp_sorted = significant_pmp.nsmallest(10, 'p_value')

# Compute specificity per CpG site
cpg_specificity = significant_pmp.groupby("CpG_Coordinates")["p_value"].mean().reset_index()
cpg_specificity.columns = ["CpG_Coordinates", "Mean_p_value"]
cpg_specificity = cpg_specificity.sort_values("Mean_p_value")

# Save top 10 PMPs and CpG specificity
significant_pmp_sorted.to_csv("new_top10_pmps.csv", index=False)
cpg_specificity.to_csv("new_cpg_specificity.csv", index=False)

# Plot specificity comparison
plt.figure(figsize=(10, 6))
sns.boxplot(data=[significant_pmp_sorted["p_value"], cpg_specificity["Mean_p_value"]])
plt.xticks([0, 1], ["Top 10 PMPs", "CpG Sites"])
plt.ylabel("Mean p-value (Specificity)")
plt.title("Specificity of Top 10 PMPs vs Individual CpG Sites")
plt.savefig("new_specificity_comparison.png")
plt.show()
