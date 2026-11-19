
# =====================================================
# STEP 4.1: Volcano Plot
# =====================================================

library(ggplot2)

# Load differential expression results
deg_results <- read.csv(
  "results/tables/all_differential_expression_results.csv"
)

# Classify genes
deg_results$regulation <- "Not Significant"

deg_results$regulation[
  deg_results$adj.P.Val < 0.05 &
    deg_results$logFC > 1
] <- "Upregulated"

deg_results$regulation[
  deg_results$adj.P.Val < 0.05 &
    deg_results$logFC < -1
] <- "Downregulated"

# Avoid log10(0)
deg_results$adj.P.Val[
  deg_results$adj.P.Val == 0
] <- .Machine$double.xmin

# Create volcano plot
volcano_plot <- ggplot(
  deg_results,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = regulation
  )
) +
  geom_point(alpha = 0.65, size = 1.2) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey70"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Gene Expression in Diabetic Wounds",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Gene Regulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# Display plot
print(volcano_plot)

# Save plot
ggsave(
  "results/plots/volcano_plot.png",
  volcano_plot,
  width = 10,
  height = 7,
  dpi = 300
)

cat("\nVolcano plot saved successfully.\n")

# =====================================================
# END OF STEP 4.1
# =====================================================


# =====================================================
# STEP 4.2: Heatmap of Top 50 Differentially Expressed Genes
# =====================================================

# Install pheatmap if it is not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# Load significant DEG results
significant_degs <- read.csv(
  "results/tables/significant_DEGs.csv"
)

# Select top 50 genes by adjusted p-value
top_genes <- significant_degs[
  order(significant_degs$adj.P.Val),
  "gene_id"
][1:50]

# Load original expression dataset
file_path <- "data/raw/GSE199939/GSE199939_gene.tpm.matrix.annot_-_21.txt.gz"

data <- read.delim(
  gzfile(file_path),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Extract expression matrix
expression_matrix <- as.matrix(data[, 5:25])
rownames(expression_matrix) <- data$gene_id

# Apply the same gene filtering used in differential expression
keep_genes <- rowSums(expression_matrix > 1) >= 3

filtered_expression <- expression_matrix[keep_genes, ]

# Log2 transformation
log_expression <- log2(filtered_expression + 1)

# Extract expression values for top 50 genes
heatmap_matrix <- log_expression[top_genes, ]

# Create sample annotation
sample_annotation <- data.frame(
  Group = ifelse(
    grepl("^DW", colnames(heatmap_matrix)),
    "Diabetic Wound",
    "Control"
  )
)

rownames(sample_annotation) <- colnames(heatmap_matrix)

# Generate heatmap
pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = sample_annotation,
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 50 Differentially Expressed Genes",
  fontsize = 10,
  filename = "results/plots/heatmap_top50_DEGs.png",
  width = 11,
  height = 10
)

cat("\nHeatmap generated and saved successfully.\n")

# =====================================================
# END OF STEP 4.2
# =====================================================


# =====================================================
# STEP 4.3: Identify Candidate Genes
# =====================================================

# Load significant DEG results
significant_degs <- read.csv(
  "results/tables/significant_DEGs.csv"
)

# Select top 10 upregulated genes
top_upregulated <- significant_degs[
  significant_degs$logFC > 1,
]

top_upregulated <- top_upregulated[
  order(top_upregulated$adj.P.Val),
][1:10, ]

# Select top 10 downregulated genes
top_downregulated <- significant_degs[
  significant_degs$logFC < -1,
]

top_downregulated <- top_downregulated[
  order(top_downregulated$adj.P.Val),
][1:10, ]

# Combine both groups
candidate_genes <- rbind(
  top_upregulated,
  top_downregulated
)

# Select important columns
candidate_genes <- candidate_genes[
  ,
  c(
    "gene_id",
    "gene_name",
    "description",
    "logFC",
    "adj.P.Val",
    "AveExpr"
  )
]

# Display candidate genes
cat("\nTop 10 Upregulated Genes:\n")
print(top_upregulated[, c(
  "gene_id", "gene_name", "logFC", "adj.P.Val"
)])

cat("\nTop 10 Downregulated Genes:\n")
print(top_downregulated[, c(
  "gene_id", "gene_name", "logFC", "adj.P.Val"
)])

# Save candidate gene table
write.csv(
  candidate_genes,
  "results/tables/candidate_genes_top10_up_down.csv",
  row.names = FALSE
)

cat("\nCandidate gene table saved successfully.\n")

# =====================================================
# END OF STEP 4.3
# =====================================================

