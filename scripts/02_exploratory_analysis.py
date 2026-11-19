
# =====================================================
# Project: Diabetic Wound Bioinformatics
# Week 2: Exploratory Data Analysis
# Step 1: Load biological dataset
# =====================================================

import pandas as pd
import numpy as np

# Define dataset path
file_path = "data/raw/GSE199939/GSE199939_gene.tpm.matrix.annot_-_21.txt.gz"

# Load expression matrix
data = pd.read_csv(
    file_path,
    sep="\t",
    compression="gzip"
)

# Display dataset dimensions
print("\nDataset dimensions:")
print(data.shape)

# Display first five rows
print("\nFirst five rows:")
print(data.head())

# Display column names
print("\nColumn names:")
print(data.columns.tolist())

# Display data types
print("\nData types:")
print(data.dtypes)

# =====================================================
# END OF STEP 1
# =====================================================
# =====================================================
# STEP 2: Check missing values and zero expression
# =====================================================

# Select only expression columns
sample_columns = data.columns[4:]

expression_data = data[sample_columns]

# Check missing values
print("\nMissing values per sample:")
print(expression_data.isnull().sum())

print("\nTotal missing values:")
print(expression_data.isnull().sum().sum())

# Check zero expression values
print("\nZero expression values per sample:")
print((expression_data == 0).sum())

print("\nTotal zero expression values:")
print((expression_data == 0).sum().sum())

# Check genes with zero expression in every sample
all_zero_genes = (expression_data == 0).all(axis=1).sum()

print("\nGenes with zero expression in all samples:")
print(all_zero_genes)

# =====================================================
# END OF STEP 2
# =====================================================
# =====================================================
# STEP 3: Check gene identifiers
# =====================================================

print("\nUnique gene IDs:")
print(data["gene_id"].nunique())

print("\nDuplicate gene IDs:")
print(data["gene_id"].duplicated().sum())

print("\nFirst three gene IDs:")
print(data["gene_id"].head(3).tolist())

print("\nLast three gene IDs:")
print(data["gene_id"].tail(3).tolist())

print("\nMissing gene names:")
print(data["gene_name"].isna().sum())

# =====================================================
# STEP 4: Expression distribution analysis
# =====================================================

import matplotlib.pyplot as plt
import seaborn as sns

# Extract expression data
sample_columns = data.columns[4:]
expression_data = data[sample_columns]

# Calculate descriptive statistics
print("\nExpression summary statistics:")
print(expression_data.describe())

# Log2 transformation for visualization
log_expression = np.log2(expression_data + 1)

# Plot expression distributions
plt.figure(figsize=(14, 6))

sns.boxplot(data=log_expression)

plt.title("Gene Expression Distribution Across Samples")
plt.xlabel("Samples")
plt.ylabel("log2(TPM + 1)")
plt.xticks(rotation=45)
plt.tight_layout()

# Save plot
plt.savefig(
    "results/plots/expression_distribution.png",
    dpi=300
)

plt.show()

# =====================================================
# END OF STEP 4
# =====================================================


# =====================================================
# STEP 5: Principal Component Analysis (PCA)
# =====================================================

from sklearn.decomposition import PCA
import os

# Select expression columns
sample_columns = data.columns[4:]
expression_data = data[sample_columns]

# Filter genes:
# Keep genes with TPM > 1 in at least 3 samples
gene_filter = (expression_data > 1).sum(axis=1) >= 3
filtered_expression = expression_data.loc[gene_filter]

print("\nPCA analysis:")
print("Genes before filtering:", expression_data.shape[0])
print("Genes after filtering:", filtered_expression.shape[0])

# Log2 transformation
log_expression = np.log2(filtered_expression + 1)

# Transpose:
# Rows = samples, columns = genes
sample_matrix = log_expression.T

# Perform PCA (features are centered automatically)
pca = PCA(n_components=2)
pca_results = pca.fit_transform(sample_matrix)

# Create PCA results dataframe
pca_df = pd.DataFrame(
    pca_results,
    columns=["PC1", "PC2"]
)

pca_df["Sample"] = sample_columns
pca_df["Group"] = [
    "Diabetic Wound" if sample.startswith("DW") else "Control"
    for sample in sample_columns
]

# Print explained variance
print("\nExplained variance:")
print("PC1:", round(pca.explained_variance_ratio_[0] * 100, 2), "%")
print("PC2:", round(pca.explained_variance_ratio_[1] * 100, 2), "%")
print(
    "Total:",
    round(pca.explained_variance_ratio_[:2].sum() * 100, 2),
    "%"
)

# Save PCA results
os.makedirs("results/plots", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

pca_df.to_csv("results/tables/pca_scores.csv", index=False)

# Plot PCA
plt.figure(figsize=(10, 7))

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="Group",
    style="Group",
    s=120
)

# Add sample labels
for _, row in pca_df.iterrows():
    plt.text(
        row["PC1"],
        row["PC2"],
        row["Sample"],
        fontsize=8,
        ha="right"
    )

plt.title("PCA of Diabetic Wound and Control Samples")
plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)"
)
plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)"
)

plt.tight_layout()

plt.savefig(
    "results/plots/PCA_diabetic_wound.png",
    dpi=300
)

plt.show()

# =====================================================
# END OF STEP 5
# =====================================================
