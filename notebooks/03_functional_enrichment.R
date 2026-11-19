
# =====================================================
# Project: Diabetic Wound Bioinformatics
# Week 3: Functional Enrichment Analysis
# Step 1: Install and load packages
# =====================================================

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Required Bioconductor packages
packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

# Install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Display package loading status
cat("\nAll enrichment packages loaded successfully.\n")

cat("\nclusterProfiler version:\n")
print(packageVersion("clusterProfiler"))

cat("\norg.Hs.eg.db version:\n")
print(packageVersion("org.Hs.eg.db"))

cat("\nenrichplot version:\n")
print(packageVersion("enrichplot"))

# =====================================================
# END OF STEP 1
# =====================================================


# =====================================================
# STEP 2: Load DEG tables and convert gene IDs
# =====================================================

# Load exported differential expression results
all_genes <- read.csv(
  "results/tables/all_differential_expression_results.csv"
)

significant_degs <- read.csv(
  "results/tables/significant_DEGs.csv"
)

# Extract unique Ensembl gene IDs
all_gene_ids <- unique(all_genes$gene_id)

significant_gene_ids <- unique(significant_degs$gene_id)

cat("\nTotal genes in background:\n")
print(length(all_gene_ids))

cat("\nTotal significant genes:\n")
print(length(significant_gene_ids))

# Convert Ensembl IDs to Entrez IDs
# Use human gene annotation database

background_mapping <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

significant_mapping <- bitr(
  significant_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate Entrez IDs
background_mapping <- background_mapping[
  !duplicated(background_mapping$ENTREZID),
]

significant_mapping <- significant_mapping[
  !duplicated(significant_mapping$ENTREZID),
]

# Display mapping statistics
cat("\nBackground genes mapped to Entrez IDs:\n")
print(nrow(background_mapping))

cat("\nSignificant genes mapped to Entrez IDs:\n")
print(nrow(significant_mapping))

cat("\nBackground mapping percentage:\n")
print(round(
  nrow(background_mapping) / length(all_gene_ids) * 100,
  2
))

cat("\nSignificant gene mapping percentage:\n")
print(round(
  nrow(significant_mapping) / length(significant_gene_ids) * 100,
  2
))

# Save converted gene ID tables
write.csv(
  background_mapping,
  "results/tables/background_gene_mapping.csv",
  row.names = FALSE
)

write.csv(
  significant_mapping,
  "results/tables/significant_gene_mapping.csv",
  row.names = FALSE
)

cat("\nGene ID conversion completed.\n")

# =====================================================
# END OF STEP 2
# =====================================================

# =====================================================
# STEP 2.4: Check mapping quality
# =====================================================

# Recreate mappings before duplicate removal
background_mapping_raw <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

significant_mapping_raw <- bitr(
  significant_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Count unique Ensembl IDs successfully mapped
background_mapped_ensembl <- unique(
  background_mapping_raw$ENSEMBL
)

significant_mapped_ensembl <- unique(
  significant_mapping_raw$ENSEMBL
)

# Calculate actual mapping coverage
background_coverage <- length(background_mapped_ensembl) /
  length(all_gene_ids) * 100

significant_coverage <- length(significant_mapped_ensembl) /
  length(significant_gene_ids) * 100

cat("\nActual background mapping coverage:\n")
print(round(background_coverage, 2))

cat("\nActual significant DEG mapping coverage:\n")
print(round(significant_coverage, 2))

# Count unmapped Ensembl IDs
cat("\nUnmapped background genes:\n")
print(length(all_gene_ids) - length(background_mapped_ensembl))

cat("\nUnmapped significant genes:\n")
print(length(significant_gene_ids) - length(significant_mapped_ensembl))

# Check one-to-many mappings
background_duplicates <- sum(
  duplicated(background_mapping_raw$ENSEMBL)
)

significant_duplicates <- sum(
  duplicated(significant_mapping_raw$ENSEMBL)
)

cat("\nAdditional mapping rows due to repeated Ensembl IDs:\n")
print(background_duplicates)

cat("\nAdditional mapping rows among significant DEGs:\n")
print(significant_duplicates)

# =====================================================
# END OF STEP 2.4
# =====================================================

# =====================================================
# STEP 2.5: Prepare final Entrez gene lists
# =====================================================

# Remove duplicate Entrez IDs
background_mapping <- background_mapping_raw[
  !duplicated(background_mapping_raw$ENTREZID),
]

significant_mapping <- significant_mapping_raw[
  !duplicated(significant_mapping_raw$ENTREZID),
]

# Create unique Entrez ID lists
background_entrez <- unique(background_mapping$ENTREZID)

significant_entrez <- unique(significant_mapping$ENTREZID)

# Check final list sizes
cat("\nFinal background Entrez IDs:\n")
print(length(background_entrez))

cat("\nFinal significant DEG Entrez IDs:\n")
print(length(significant_entrez))

# Check for duplicate Entrez IDs
cat("\nDuplicate background Entrez IDs:\n")
print(sum(duplicated(background_entrez)))

cat("\nDuplicate significant Entrez IDs:\n")
print(sum(duplicated(significant_entrez)))

# Save final mapping tables
write.csv(
  background_mapping,
  "results/tables/background_gene_mapping.csv",
  row.names = FALSE
)

write.csv(
  significant_mapping,
  "results/tables/significant_gene_mapping.csv",
  row.names = FALSE
)

# =====================================================
# END OF STEP 2.5
# =====================================================

# =====================================================
# STEP 3.1: GO Biological Process Enrichment
# =====================================================

go_bp <- enrichGO(
  gene = significant_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Display enrichment summary
cat("\nGO Biological Process enrichment results:\n")

print(go_bp)

# Number of enriched biological processes
cat("\nNumber of enriched GO terms:\n")

print(nrow(as.data.frame(go_bp)))

# Display top 10 enriched terms
go_bp_results <- as.data.frame(go_bp)

head(
  go_bp_results[
    order(go_bp_results$p.adjust),
    c("ID", "Description", "GeneRatio",
      "BgRatio", "p.adjust", "Count")
  ],
  10
)

# =====================================================
# END OF STEP 3.1
# =====================================================
# =====================================================
# STEP 3.2: Export and visualize GO BP results
# =====================================================

# Export complete GO enrichment results
write.csv(
  go_bp_results,
  "results/tables/GO_BP_enrichment_results.csv",
  row.names = FALSE
)

# Create dot plot of top 20 enriched terms
go_bp_plot <- enrichplot::dotplot(
  go_bp,
  showCategory = 20,
  orderBy = "p.adjust",
  font.size = 10,
  title = "GO Biological Process Enrichment"
)

# Display plot
print(go_bp_plot)

# Save plot
ggplot2::ggsave(
  filename = "results/plots/GO_BP_dotplot.png",
  plot = go_bp_plot,
  width = 12,
  height = 9,
  dpi = 300
)

cat("\nGO results and plot saved successfully.\n")

# =====================================================
# END OF STEP 3.2
# =====================================================
# =====================================================
# STEP 3.3: KEGG Pathway Enrichment
# =====================================================

kegg_results <- enrichKEGG(
  gene = significant_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Display KEGG enrichment results
cat("\nKEGG enrichment results:\n")

print(kegg_results)

# Convert results into a data frame
kegg_df <- as.data.frame(kegg_results)

# Count enriched pathways
cat("\nNumber of enriched KEGG pathways:\n")

print(nrow(kegg_df))

# Display top 10 pathways
head(
  kegg_df[
    order(kegg_df$p.adjust),
    c("ID", "Description", "GeneRatio",
      "BgRatio", "p.adjust", "Count")
  ],
  10
)

# =====================================================
# END OF STEP 3.3
# =====================================================

# =====================================================
# STEP 3.4: Export and visualize KEGG results
# =====================================================

# Export complete KEGG results
write.csv(
  kegg_df,
  "results/tables/KEGG_enrichment_results.csv",
  row.names = FALSE
)

# Create dot plot of top 20 pathways
kegg_plot <- enrichplot::dotplot(
  kegg_results,
  showCategory = 20,
  orderBy = "p.adjust",
  font.size = 10,
  title = "KEGG Pathway Enrichment"
)

# Display plot
print(kegg_plot)

# Save plot
ggplot2::ggsave(
  filename = "results/plots/KEGG_dotplot.png",
  plot = kegg_plot,
  width = 12,
  height = 9,
  dpi = 300
)

cat("\nKEGG results and plot saved successfully.\n")

# =====================================================
# END OF STEP 3.4
# =====================================================

# =====================================================
# STEP 3.5: Prepare Upregulated and Downregulated
# Entrez Gene Lists
# =====================================================

# Read upregulated and downregulated DEG tables
upregulated_degs <- read.csv(
  "results/tables/upregulated_genes.csv"
)

downregulated_degs <- read.csv(
  "results/tables/downregulated_genes.csv"
)

# Extract unique Ensembl IDs
upregulated_ids <- unique(upregulated_degs$gene_id)

downregulated_ids <- unique(downregulated_degs$gene_id)

# Convert Ensembl IDs to Entrez IDs
upregulated_mapping <- bitr(
  upregulated_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

downregulated_mapping <- bitr(
  downregulated_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate Entrez IDs
upregulated_mapping <- upregulated_mapping[
  !duplicated(upregulated_mapping$ENTREZID),
]

downregulated_mapping <- downregulated_mapping[
  !duplicated(downregulated_mapping$ENTREZID),
]

# Create final gene lists
upregulated_entrez <- unique(upregulated_mapping$ENTREZID)

downregulated_entrez <- unique(downregulated_mapping$ENTREZID)

# Display counts
cat("\nOriginal upregulated genes:\n")
print(length(upregulated_ids))

cat("\nMapped upregulated Entrez IDs:\n")
print(length(upregulated_entrez))

cat("\nOriginal downregulated genes:\n")
print(length(downregulated_ids))

cat("\nMapped downregulated Entrez IDs:\n")
print(length(downregulated_entrez))

# Save mapping tables
write.csv(
  upregulated_mapping,
  "results/tables/upregulated_gene_mapping.csv",
  row.names = FALSE
)

write.csv(
  downregulated_mapping,
  "results/tables/downregulated_gene_mapping.csv",
  row.names = FALSE
)

# =====================================================
# END OF STEP 3.5
# =====================================================

# =====================================================
# STEP 3.6: Verify directional gene mappings
# =====================================================

# Identify successfully mapped Ensembl IDs
up_mapped_ensembl <- unique(
  upregulated_mapping$ENSEMBL
)

down_mapped_ensembl <- unique(
  downregulated_mapping$ENSEMBL
)

# Calculate mapping coverage
up_coverage <- length(up_mapped_ensembl) /
  length(upregulated_ids) * 100

down_coverage <- length(down_mapped_ensembl) /
  length(downregulated_ids) * 100

cat("\nUpregulated Ensembl mapping coverage:\n")
print(round(up_coverage, 2))

cat("\nDownregulated Ensembl mapping coverage:\n")
print(round(down_coverage, 2))

# Check overlap between Entrez lists
common_entrez <- intersect(
  upregulated_entrez,
  downregulated_entrez
)

cat("\nEntrez IDs shared between both lists:\n")
print(length(common_entrez))

# Check whether all significant Entrez IDs are represented
combined_directional_entrez <- union(
  upregulated_entrez,
  downregulated_entrez
)

cat("\nUnique Entrez IDs across both lists:\n")
print(length(combined_directional_entrez))

# =====================================================
# END OF STEP 3.6
# =====================================================

# =====================================================
# STEP 3.7: Separate GO enrichment analyses
# =====================================================

# GO enrichment for upregulated genes
go_up <- enrichGO(
  gene = upregulated_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# GO enrichment for downregulated genes
go_down <- enrichGO(
  gene = downregulated_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Convert results into data frames
go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)

# Display number of enriched terms
cat("\nUpregulated GO Biological Processes:\n")
print(nrow(go_up_df))

cat("\nDownregulated GO Biological Processes:\n")
print(nrow(go_down_df))

# Display top 10 upregulated GO terms
cat("\nTop 10 GO terms for upregulated genes:\n")

head(
  go_up_df[
    order(go_up_df$p.adjust),
    c("ID", "Description", "GeneRatio",
      "p.adjust", "Count")
  ],
  10
)

# Display top 10 downregulated GO terms
cat("\nTop 10 GO terms for downregulated genes:\n")

head(
  go_down_df[
    order(go_down_df$p.adjust),
    c("ID", "Description", "GeneRatio",
      "p.adjust", "Count")
  ],
  10
)

# =====================================================
# END OF STEP 3.7
# =====================================================

# =====================================================
# STEP 3.8: Export and visualize directional GO results
# =====================================================

# Save enrichment tables
write.csv(
  go_up_df,
  "results/tables/GO_BP_upregulated_results.csv",
  row.names = FALSE
)

write.csv(
  go_down_df,
  "results/tables/GO_BP_downregulated_results.csv",
  row.names = FALSE
)

# Create upregulated GO dot plot
go_up_plot <- enrichplot::dotplot(
  go_up,
  showCategory = 15,
  orderBy = "p.adjust",
  font.size = 10,
  title = "GO BP - Upregulated Genes"
)

# Create downregulated GO dot plot
go_down_plot <- enrichplot::dotplot(
  go_down,
  showCategory = 15,
  orderBy = "p.adjust",
  font.size = 10,
  title = "GO BP - Downregulated Genes"
)

# Display plots
print(go_up_plot)
print(go_down_plot)

# Save plots
ggplot2::ggsave(
  "results/plots/GO_BP_upregulated.png",
  go_up_plot,
  width = 12,
  height = 9,
  dpi = 300
)

ggplot2::ggsave(
  "results/plots/GO_BP_downregulated.png",
  go_down_plot,
  width = 12,
  height = 9,
  dpi = 300
)

cat("\nDirectional GO results and plots saved successfully.\n")

# =====================================================
# END OF STEP 3.8
# =====================================================

# =====================================================
# STEP 3.9: Separate KEGG enrichment analyses
# =====================================================

# KEGG enrichment for upregulated genes
kegg_up <- enrichKEGG(
  gene = upregulated_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# KEGG enrichment for downregulated genes
kegg_down <- enrichKEGG(
  gene = downregulated_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Convert results into data frames
kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

# Display number of enriched pathways
cat("\nUpregulated KEGG pathways:\n")
print(nrow(kegg_up_df))

cat("\nDownregulated KEGG pathways:\n")
print(nrow(kegg_down_df))

# Display top 10 upregulated pathways
cat("\nTop 10 KEGG pathways - Upregulated genes:\n")

head(
  kegg_up_df[
    order(kegg_up_df$p.adjust),
    c("ID", "Description", "GeneRatio",
      "p.adjust", "Count")
  ],
  10
)

# Display top 10 downregulated pathways
cat("\nTop 10 KEGG pathways - Downregulated genes:\n")

head(
  kegg_down_df[
    order(kegg_down_df$p.adjust),
    c("ID", "Description", "GeneRatio",
      "p.adjust", "Count")
  ],
  10
)

# =====================================================
# END OF STEP 3.9
# =====================================================

# =====================================================
# STEP 3.10: Export and visualize directional KEGG
# =====================================================

# Save enrichment tables
write.csv(
  kegg_up_df,
  "results/tables/KEGG_upregulated_results.csv",
  row.names = FALSE
)

write.csv(
  kegg_down_df,
  "results/tables/KEGG_downregulated_results.csv",
  row.names = FALSE
)

# Create upregulated KEGG dot plot
kegg_up_plot <- enrichplot::dotplot(
  kegg_up,
  showCategory = 15,
  orderBy = "p.adjust",
  font.size = 10,
  title = "KEGG - Upregulated Genes"
)

# Create downregulated KEGG dot plot
kegg_down_plot <- enrichplot::dotplot(
  kegg_down,
  showCategory = 15,
  orderBy = "p.adjust",
  font.size = 10,
  title = "KEGG - Downregulated Genes"
)

# Display plots
print(kegg_up_plot)
print(kegg_down_plot)

# Save plots
ggplot2::ggsave(
  "results/plots/KEGG_upregulated.png",
  kegg_up_plot,
  width = 12,
  height = 9,
  dpi = 300
)

ggplot2::ggsave(
  "results/plots/KEGG_downregulated.png",
  kegg_down_plot,
  width = 12,
  height = 9,
  dpi = 300
)

cat("\nDirectional KEGG results and plots saved successfully.\n")

# =====================================================
# END OF STEP 3.10
# =====================================================
