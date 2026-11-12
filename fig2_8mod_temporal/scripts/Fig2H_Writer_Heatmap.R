# Fig2 Panel G: Writer/Reader/Enzyme x Cell Type Expression Heatmap
# =============================================================================
# AUDIT 2026-06-19: This panel previously had output but no generation script.
# This script reads server-computed enzyme-celltype Z-scores from scRNA-seq
# (GSE274819 Seurat object). Run run_enzyme_celltype_server.R on server first
# to generate the required input CSV.
# =============================================================================

library(data.table)
library(ggplot2)
suppressPackageStartupMessages(library(pheatmap))

source("/mnt/d/desktop/DRS/DRS_TechShowcase/shared/palette/color_palette.R")

out_dir <- "../output"
data_dir <- "../data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load enzyme-celltype Z-scores ----
zscore_file <- file.path(data_dir, "enzyme_celltype_zscore_server.csv")

if (file.exists(zscore_file)) {
  cat("Loading server-computed enzyme-celltype Z-scores...\n")
  enzyme_z <- fread(zscore_file)
} else {
  cat("========================================\n")
  cat("Server enzyme-celltype data NOT FOUND.\n")
  cat("Run on server: run_enzyme_celltype_server.R\n")
  cat("Then copy enzyme_celltype_zscore_server.csv\n")
  cat("to fig2_8mod_temporal/data/\n")
  cat("========================================\n\n")

  # Generate from DRS TPM + enrichment data as best-available proxy
  # Load enzyme Day1 FC data (real DRS TPM)
  enzyme_file <- file.path(data_dir, "Table_Enzyme_Day1_FC.csv")
  if (!file.exists(enzyme_file)) {
    stop("Table_Enzyme_Day1_FC.csv not found. Cannot generate heatmap.")
  }

  cat("Using DRS TPM expression data as proxy (per-group mean, not cell-type-resolved).\n")
  cat("WARNING: This is bulk RNA expression, not cell-type-specific.\n")
  cat("Replace with scRNA server output when available.\n\n")

  enzyme_dt <- fread(enzyme_file)

  # Also load the full TPM matrix to get per-sample values
  tpm_file <- file.path(data_dir, "gene.TPM.matrix.xls")
  if (!file.exists(tpm_file)) {
    stop("gene.TPM.matrix.xls not found.")
  }

  drs_tpm <- fread(tpm_file)
  setnames(drs_tpm, 1, "gene_id")

  # Map enzyme symbols to gene IDs (using add_info)
  info_file <- file.path(data_dir, "gene.TPM.matrix.add_info.xls")
  drs_info <- fread(info_file, fill = TRUE)
  gn_map <- setNames(drs_info$ID, drs_info$gene_name)

  # Get per-group TPM for each enzyme
  sample_cols <- setdiff(colnames(drs_tpm), "gene_id")
  groups <- gsub("s[0-9]$", "", sample_cols)
  groups[groups == "I"] <- "IS"
  unique_groups <- c("NC", "IS", "Day1", "Day3", "Day7")

  enzyme_expr <- rbindlist(lapply(enzyme_dt$enzyme, function(sym) {
    gid <- gn_map[sym]
    if (is.na(gid) || !gid %in% drs_tpm$gene_id) return(NULL)
    vals <- as.numeric(drs_tpm[gene_id == gid, ..sample_cols])
    data.table(
      enzyme = sym,
      NC = mean(vals[groups == "NC"]),
      IS = mean(vals[groups == "IS"]),
      Day1 = mean(vals[groups == "Day1"]),
      Day3 = mean(vals[groups == "Day3"]),
      Day7 = mean(vals[groups == "Day7"])
    )
  }))

  if (nrow(enzyme_expr) == 0) stop("No enzyme expression data found in TPM matrix.")

  # Log2(TPM+1) and Z-score
  expr_mat <- as.matrix(enzyme_expr[, -"enzyme"])
  rownames(expr_mat) <- enzyme_expr$enzyme
  expr_mat_log <- log2(expr_mat + 1)
  expr_z <- t(scale(t(expr_mat_log)))

  enzyme_z <- as.data.table(expr_z, keep.rownames = "Gene")
  enzyme_z[, Category := enzyme_dt$Category[match(Gene, enzyme_dt$enzyme)]]

  # Write the interim data
  fwrite(enzyme_z, file.path(data_dir, "enzyme_expression_zscore_DRS_TPM.csv"))
  cat("Saved interim data to enzyme_expression_zscore_DRS_TPM.csv\n")
}

# ---- Prepare heatmap matrix ----
if ("Gene" %in% colnames(enzyme_z)) {
  gene_col <- "Gene"
} else {
  gene_col <- "enzyme"
}

heat_rows <- enzyme_z[[gene_col]]
heat_mat <- as.matrix(enzyme_z[, -c(gene_col, "Category"), with = FALSE])
rownames(heat_mat) <- heat_rows

# Handle duplicate gene names (e.g., Dkc1 in both Nm and PSU)
dup_genes <- heat_rows[duplicated(heat_rows)]
if (length(dup_genes) > 0) {
  cat("Duplicate gene names found:", paste(dup_genes, collapse=", "), "\n")
  cat("Making unique by appending category suffix...\n")
  for (i in seq_along(heat_rows)) {
    is_dup <- sum(heat_rows == heat_rows[i]) > 1
    if (is_dup) {
      rownames(heat_mat)[i] <- paste0(heat_rows[i], " (", enzyme_z$Category[i], ")")
    }
  }
}

# Category annotation
cat_vec <- setNames(enzyme_z$Category, rownames(heat_mat))

# Remove duplicate row names (keep first occurrence)
if (any(duplicated(rownames(heat_mat)))) {
  dup_idx <- which(duplicated(rownames(heat_mat)))
  cat("Removing", length(dup_idx), "duplicate row(s):", paste(rownames(heat_mat)[dup_idx], collapse=", "), "\n")
  heat_mat <- heat_mat[-dup_idx, , drop = FALSE]
  cat_vec <- cat_vec[rownames(heat_mat)]
}

# Remove rows with all NA
keep_rows <- rowSums(is.na(heat_mat)) == 0
heat_mat <- heat_mat[keep_rows, , drop = FALSE]
cat_vec <- cat_vec[rownames(heat_mat)]

cat_colors <- list(Category = c(
  "m6A Writer"     = "#CC5A4C",  # terracotta (m6A)
  "m6A Reader"     = "#E8948A",  # light red
  "m6A Eraser"     = "#F4A460",  # sandy brown
  "m5C Writer"     = "#7B3F8C",  # purple (distinct from reds)
  "Nm Writer" = "#3B4F8A",  # deep blue (Gm)
  "PSU Synthase"   = "#00A087",  # green (PSU)
  "Demethylase"    = "#5C6B9A"   # slate blue (Um)
))

# Ensure cat_vec names match heat_mat rownames
cat_vec <- cat_vec[rownames(heat_mat)]
ann_row <- data.frame(Category = factor(cat_vec, levels = names(cat_colors$Category)),
                      row.names = rownames(heat_mat))

# ---- Plot heatmap ----
is_scRNA <- any(grepl("^(PT|TAL|DCT|Macrophage|B_cell)$", colnames(heat_mat)))
title_msg <- if (is_scRNA) {
  "Writer/Reader Enzyme Expression x Cell Type (scRNA-seq Z-score)"
} else {
  "Writer/Reader Enzyme Expression Across IRI Timecourse (DRS TPM Z-score)\n[WARNING: bulk TPM, not cell-type-resolved. Replace with scRNA data.]"
}

p <- pheatmap(heat_mat,
  annotation_row = ann_row,
  annotation_colors = cat_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  fontsize_row = 8,
  fontsize_col = 9,
  color = colorRampPalette(c("#3B4F8A", "white", "#CC5A4C"))(100),
  main = title_msg,
  angle_col = 45,
  border_color = "white",
  silent = TRUE
)

# canonical Fig2H output (letter per composite v5 A-M roster)
dir.create("../outputs", showWarnings = FALSE, recursive = TRUE)
for (od in c(out_dir, "../outputs")) {
  pdf(file.path(od, "Fig2H_Writer_Heatmap.pdf"), width = 12, height = 9)
  grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
  png(file.path(od, "Fig2H_Writer_Heatmap.png"), width = 12, height = 9, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
}
cat("\nFig2H_Writer_Heatmap saved (output/ + outputs/).\n")
