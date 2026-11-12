#!/usr/bin/env Rscript
# Server-side: Extract enzyme × cell-type data for Fig2 Panel G (Writer_Heatmap) + Panel L (Composition Drift)
# =============================================================================
# Outputs:
#   1. cell_proportions.csv          — real cell counts per type, Sham vs Day1
#   2. enzyme_celltype_expression.csv — 28 enzymes × 12 cell types, mean expr
#   3. enzyme_primary_celltype.csv   — primary cell type per enzyme (highest expr)
#   4. enzyme_celltype_zscore.csv    — row-Z-scored for heatmap
# =============================================================================

library(Seurat)
library(dplyr)
suppressPackageStartupMessages(library(data.table))

so <- readRDS('~/scRNA_analysis/results/seurat_merged_annotated.rds')
DefaultAssay(so) <- 'RNA'

# ---- 1. Cell Type Proportions (Sham vs Day1) ----
cat("Computing cell proportions...\n")
# Filter to Sham and Day1
so_sub <- subset(so, subset = orig.ident %in% c("Sham", "Day1"))
cell_counts <- table(so_sub$celltype, so_sub$orig.ident)
cell_props_dt <- data.table(
  CellType = rownames(cell_counts),
  Sham_Count = as.integer(cell_counts[, "Sham"]),
  Day1_Count = as.integer(cell_counts[, "Day1"])
)
cell_props_dt[, Sham_Pct := round(Sham_Count / sum(Sham_Count) * 100, 1)]
cell_props_dt[, Day1_Pct := round(Day1_Count / sum(Day1_Count) * 100, 1)]
cell_props_dt[, Delta_Pct := Day1_Pct - Sham_Pct]
fwrite(cell_props_dt, '~/scRNA_analysis/results/cell_proportions.csv')
cat(sprintf("Cell proportions saved. Sham n=%d, Day1 n=%d\n",
            sum(cell_props_dt$Sham_Count), sum(cell_props_dt$Day1_Count)))
print(cell_props_dt[, .(CellType, Sham_Count, Day1_Count, Delta_Pct)])

# ---- 2. Enzyme Expression per Cell Type ----
cat("\nComputing enzyme expression per cell type...\n")
Idents(so) <- so$celltype

# 28 modification enzymes
enzyme_genes <- c(
  # m6A Writers
  "Mettl3", "Mettl14", "Wtap", "Rbm15", "Rbm15b", "Zc3h13", "Virma", "Cbll1",
  # m6A Readers
  "Hnrnpa2b1", "Hnrnpc", "Ythdc1", "Ythdc2", "Ythdf1", "Ythdf2", "Ythdf3",
  # m6A Erasers
  "Fto", "Alkbh5",
  # m5C Writers
  "Nsun2", "Nsun3", "Nsun5", "Trdmt1",
  # Nm Writers
  "Fbl", "Nop56", "Nop58", "Dkc1", "Ftsj3",
  # PSU Synthases
  "Pus1", "Pus3", "Pus7", "Pus10",
  # DNA/RNA Demethylase
  "Tet2"
)

# Category mapping
enzyme_cat <- list(
  "m6A Writer"      = c("Mettl3","Mettl14","Wtap","Rbm15","Rbm15b","Zc3h13","Virma","Cbll1"),
  "m6A Reader"      = c("Hnrnpa2b1","Hnrnpc","Ythdc1","Ythdc2","Ythdf1","Ythdf2","Ythdf3"),
  "m6A Eraser"      = c("Fto","Alkbh5"),
  "m5C Writer"      = c("Nsun2","Nsun3","Nsun5","Trdmt1"),
  "Nm Writer"  = c("Fbl","Nop56","Nop58","Dkc1","Ftsj3"),
  "PSU Synthase"    = c("Pus1","Pus3","Pus7","Pus10"),
  "Demethylase"     = c("Tet2")
)
cat_map <- setNames(rep(names(enzyme_cat), lengths(enzyme_cat)), unlist(enzyme_cat))

# Get AverageExpression
avg_expr <- AverageExpression(so, assays = 'RNA', slot = 'data', group.by = 'celltype')$RNA
# Try to match gene symbols (mouse: capitalize first letter)
# Seurat uses gene symbols; mouse = Mettl3, human = METTL3. Try both.
found_genes <- intersect(enzyme_genes, rownames(avg_expr))
# If none found, try case-insensitive
if (length(found_genes) == 0) {
  found_genes <- intersect(tolower(enzyme_genes), tolower(rownames(avg_expr)))
}
cat(sprintf("Found %d/%d enzyme genes in scRNA data\n", length(found_genes), length(enzyme_genes)))
cat("Missing:", paste(setdiff(enzyme_genes, found_genes), collapse=", "), "\n")

# Subset to found genes
enzyme_expr <- avg_expr[found_genes, , drop = FALSE]

# Add category annotation
enzyme_dt <- as.data.table(enzyme_expr, keep.rownames = "Gene")
enzyme_dt[, Category := cat_map[Gene]]

# Save raw expression
fwrite(enzyme_dt, '~/scRNA_analysis/results/enzyme_celltype_expression.csv')

# Z-score normalize within each gene (row)
enzyme_mat <- as.matrix(enzyme_dt[, -c("Gene", "Category")])
rownames(enzyme_mat) <- enzyme_dt$Gene
enzyme_z <- t(scale(t(enzyme_mat)))

enzyme_z_dt <- as.data.table(enzyme_z, keep.rownames = "Gene")
enzyme_z_dt[, Category := cat_map[Gene]]
fwrite(enzyme_z_dt, '~/scRNA_analysis/results/enzyme_celltype_zscore.csv')

# ---- 3. Primary Cell Type per Enzyme ----
cat("\nDetermining primary cell type per enzyme...\n")
primary_ct <- data.table(
  enzyme = rownames(enzyme_mat),
  primary_celltype = colnames(enzyme_mat)[apply(enzyme_mat, 1, which.max)],
  max_expression = apply(enzyme_mat, 1, max)
)
primary_ct[, Category := cat_map[enzyme]]
fwrite(primary_ct, '~/scRNA_analysis/results/enzyme_primary_celltype.csv')

cat("\nPrimary cell types:\n")
print(primary_ct)

cat("\n=== All done ===\n")
cat("Output files:\n")
cat("  ~/scRNA_analysis/results/cell_proportions.csv\n")
cat("  ~/scRNA_analysis/results/enzyme_celltype_expression.csv\n")
cat("  ~/scRNA_analysis/results/enzyme_celltype_zscore.csv\n")
cat("  ~/scRNA_analysis/results/enzyme_primary_celltype.csv\n")
