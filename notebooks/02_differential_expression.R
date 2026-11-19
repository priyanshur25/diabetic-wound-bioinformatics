# =====================================================
# Project: Diabetic Wound Bioinformatics
# Week 2: Differential Gene Expression Analysis
# Step 1: Prepare expression matrix
# =====================================================

# Install BiocManager if required
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install limma if required
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

# Load package
library(limma)

# =====================================================
# 1. Load original dataset
# =====================================================

file_path <- "data/raw/GSE199939/GSE199939_gene.tpm.matrix.annot_-_21.txt.gz"

data <- read.delim(
  gzfile(file_path),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

cat("Original dataset dimensions:\n")
print(dim(data))

# =====================================================
# 2. Extract gene information
# =====================================================

gene_info <- data[, 1:4]

# Extract expression values
expression_data <- data[, 5:25]

# Convert to matrix
expression_matrix <- as.matrix(expression_data)

# Store gene IDs as row names
rownames(expression_matrix) <- data$gene_id

# =====================================================
# 3. Define sample groups
# =====================================================

sample_names <- colnames(expression_matrix)

group <- factor(
  ifelse(grepl("^DW", sample_names),
         "DiabeticWound",
         "Control"),
  levels = c("Control", "DiabeticWound")
)

# Display group information
cat("\nSample groups:\n")
print(table(group))

# =====================================================
# 4. Filter low-expression genes
# =====================================================

# Retain genes with TPM > 1 in at least 3 samples
keep_genes <- rowSums(expression_matrix > 1) >= 3

filtered_expression <- expression_matrix[keep_genes, ]

cat("\nGenes before filtering:\n")
print(nrow(expression_matrix))

cat("\nGenes after filtering:\n")
print(nrow(filtered_expression))

# =====================================================
# 5. Log2 transformation
# =====================================================

log_expression <- log2(filtered_expression + 1)

cat("\nFinal expression matrix dimensions:\n")
print(dim(log_expression))

# =====================================================
# 6. Check sample group alignment
# =====================================================

cat("\nSample names:\n")
print(colnames(log_expression))

cat("\nGroup assignments:\n")
print(group)

# =====================================================
# END OF STEP 1
# =====================================================

# =====================================================
# STEP 2: Create statistical design matrix
# =====================================================

# Create design matrix without an intercept
design <- model.matrix(~ 0 + group)

# Rename columns
colnames(design) <- levels(group)

# Display design matrix
cat("\nStatistical design matrix:\n")
print(design)

# Check dimensions
cat("\nDesign matrix dimensions:\n")
print(dim(design))

# Check group sizes
cat("\nNumber of samples per group:\n")
print(colSums(design))

# =====================================================
# END OF STEP 2
# =====================================================

# =====================================================
# STEP 3: Define group comparison
# =====================================================

contrast_matrix <- makeContrasts(
  DiabeticWound - Control,
  levels = design
)

cat("\nContrast matrix:\n")
print(contrast_matrix)

# =====================================================
# END OF STEP 3
# =====================================================


# =====================================================
# STEP 4: Fit limma linear model
# =====================================================

# Fit linear model for each gene
fit <- lmFit(log_expression, design)

# Apply the diabetic wound vs control comparison
fit_contrasts <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes moderation
fit_ebayes <- eBayes(fit_contrasts)

# Display model summary
cat("\nLimma model fitted successfully.\n")

cat("\nNumber of genes analysed:\n")
print(nrow(fit_ebayes$coefficients))

cat("\nNumber of samples:\n")
print(ncol(log_expression))

cat("\nContrast used:\n")
print(colnames(fit_ebayes$coefficients))

# =====================================================
# END OF STEP 4
# =====================================================

# =====================================================
# STEP 5: Extract differential expression results
# =====================================================

deg_results <- topTable(
  fit_ebayes,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Add gene IDs
deg_results$gene_id <- rownames(deg_results)

# Add gene annotations
deg_results <- merge(
  deg_results,
  gene_info,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

# Display the first 10 results
cat("\nTop 10 differentially expressed genes:\n")
print(head(deg_results[, c(
  "gene_id",
  "gene_name",
  "logFC",
  "AveExpr",
  "P.Value",
  "adj.P.Val"
)], 10))

# Display total number of genes tested
cat("\nTotal genes tested:\n")
print(nrow(deg_results))

# =====================================================
# END OF STEP 5
# =====================================================

# =====================================================
# STEP 6: Identify significant DEGs
# =====================================================

# Apply significance thresholds
significant_degs <- deg_results[
  deg_results$adj.P.Val < 0.05 &
    abs(deg_results$logFC) > 1,
]

# Count significant genes
cat("\nTotal significant DEGs:\n")
print(nrow(significant_degs))

# Count upregulated genes
upregulated <- significant_degs[
  significant_degs$logFC > 1,
]

cat("\nUpregulated genes:\n")
print(nrow(upregulated))

# Count downregulated genes
downregulated <- significant_degs[
  significant_degs$logFC < -1,
]

cat("\nDownregulated genes:\n")
print(nrow(downregulated))

# Display top significant genes
cat("\nTop 10 significant DEGs:\n")
print(head(significant_degs[, c(
  "gene_id",
  "gene_name",
  "logFC",
  "P.Value",
  "adj.P.Val"
)], 10))

# =====================================================
# END OF STEP 6
# =====================================================

# =====================================================
# STEP 7: Export differential expression results
# =====================================================

# Create output directory if it doesn't exist
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save all tested genes
write.csv(
  deg_results,
  "results/tables/all_differential_expression_results.csv",
  row.names = FALSE
)

# Save all significant DEGs
write.csv(
  significant_degs,
  "results/tables/significant_DEGs.csv",
  row.names = FALSE
)

# Save upregulated genes
write.csv(
  upregulated,
  "results/tables/upregulated_genes.csv",
  row.names = FALSE
)

# Save downregulated genes
write.csv(
  downregulated,
  "results/tables/downregulated_genes.csv",
  row.names = FALSE
)

cat("\nAll DEG result files saved successfully.\n")

# =====================================================
# END OF STEP 7
# =====================================================
