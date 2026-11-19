
# =====================================================
# Project: Diabetic Wound Bioinformatics
# Step 1: Download biological dataset
# GEO Accession: GSE199939
# =====================================================

# Load library
library(GEOquery)

# Download supplementary files from GEO
getGEOSuppFiles(
  "GSE199939",
  makeDirectory = TRUE
)

# Find expression matrix
files <- list.files(
  "GSE199939",
  pattern = "gene.tpm.matrix",
  full.names = TRUE
)

# Display downloaded file
print(files)


# Read expression matrix
data <- read.delim(
  files[1],
  check.names = FALSE
)

# Inspect dataset
cat("\nDataset dimensions:\n")
print(dim(data))

cat("\nFirst six rows:\n")
print(head(data))

cat("\nColumn names:\n")
print(colnames(data))

cat("\nFirst six gene identifiers:\n")
print(head(data[[1]]))

# 7. Inspect first ten gene identifiers
cat("\nFirst ten gene identifiers:\n")
print(head(data[[1]], 10))

# 8. Inspect data types
cat("\nColumn data types:\n")
print(sapply(data, class))

# 9. Check missing values
cat("\nMissing values per column:\n")
print(colSums(is.na(data)))

# 10. Check duplicate gene identifiers
cat("\nDuplicate first-column identifiers:\n")
print(sum(duplicated(data[[1]])))

# Save a copy of the inspection report
sink("reports/dataset_inspection.txt")

cat("GSE199939 Dataset Inspection Report\n")
cat("===================================\n")
cat("Dimensions:", dim(data), "\n\n")
cat("Column names:\n")
print(colnames(data))
cat("\nFirst six rows:\n")
print(head(data))

sink()

cat("\nDataset download and inspection completed!\n")



# =====================================================
# FINAL WEEK 1 DATASET INSPECTION REPORT
# =====================================================

sink("reports/dataset_inspection.txt")

cat("WEEK 1: DATASET INSPECTION REPORT\n")
cat("====================================\n\n")

cat("1. DATASET INFORMATION\n")
cat("GEO Accession: GSE199939\n")
cat("Source: NCBI GEO\n")
cat("Organism: Homo sapiens\n")
cat("Expression measure: TPM\n\n")

cat("2. DOWNLOAD INFORMATION\n")
cat("File: GSE199939_gene.tpm.matrix.annot_-_21.txt.gz\n")
cat("Download status: Successful\n\n")

cat("3. DATASET STRUCTURE\n")
cat("Total rows:", nrow(data), "\n")
cat("Total columns:", ncol(data), "\n")
cat("Biological samples:", ncol(data) - 4, "\n")
cat("Gene identifier: Ensembl Gene ID\n")
cat("Gene symbol: Available\n")
cat("Gene description: Available\n\n")

cat("4. SAMPLE INFORMATION\n")
cat("Diabetic wound samples: 10 (DW1-DW10)\n")
cat("Non-diabetic controls: 11 (N1-N11)\n\n")

cat("5. DATA QUALITY INSPECTION\n")
cat("Duplicate gene identifiers:",
    sum(duplicated(data$gene_id)), "\n")
cat("Total missing values:", sum(is.na(data)), "\n")
cat("Affected gene: SLC1A5\n")
cat("Affected gene ID: ENSG00000105281\n")
cat("Missing values in diabetic samples: 0\n")
cat("Missing values in control samples: 11\n\n")

cat("6. INITIAL OBSERVATIONS\n")
cat("All 11 missing expression values occur in the\n")
cat("SLC1A5 gene row across non-diabetic controls.\n")
cat("No expression values were modified or replaced.\n\n")

cat("7. CONCLUSION\n")
cat("The GSE199939 expression matrix was successfully\n")
cat("downloaded and inspected. It contains 3069 gene\n")
cat("records and 21 biological samples.\n")

sink()

cat("\nWeek 1 inspection report generated successfully!\n")

