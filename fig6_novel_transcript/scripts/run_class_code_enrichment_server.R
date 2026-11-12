#!/usr/bin/env Rscript
# Server-side: Per-Class-Code Novel Transcript Cell-Type Enrichment
# Run on: sw2-dynamic.xiyoucloud.pro
# ============================================================
# Computes scRNA-seq expression enrichment (Z-score) for each
# novel transcript class code (i/j/k/m/n/u/x) across 12 cell types.
# ============================================================

library(Seurat)
library(dplyr)
suppressPackageStartupMessages(library(data.table))

# ---- Load Seurat ----
so <- readRDS('~/scRNA_analysis/results/seurat_merged_annotated.rds')
DefaultAssay(so) <- 'RNA'
Idents(so) <- so$celltype

# ---- Load class code gene sets ----
class_codes <- c("i", "j", "k", "m", "n", "u", "x")
gene_sets <- list()
for (cc in class_codes) {
  fname <- paste0("geneset_class_", cc, "_symbols.txt")
  genes <- scan(fname, what = "", quiet = TRUE)
  gene_sets[[paste0("Class_", cc)]] <- genes
  cat(sprintf("Class %s: %d genes\n", cc, length(genes)))
}

# ---- Compute average expression per cell type ----
avg_expr <- AverageExpression(so, assays = 'RNA', slot = 'data', group.by = 'celltype')$RNA

enrich_raw <- lapply(names(gene_sets), function(gs_name) {
  genes <- intersect(gene_sets[[gs_name]], rownames(avg_expr))
  if (length(genes) == 0) return(NULL)
  colMeans(avg_expr[genes, , drop = FALSE])
})
names(enrich_raw) <- names(gene_sets)
enrich_raw <- enrich_raw[!sapply(enrich_raw, is.null)]

# ---- Z-score normalize across cell types (row-wise) ----
enrich_mat <- do.call(rbind, enrich_raw)
enrich_z <- t(scale(t(enrich_mat)))  # row z-score

# ---- Also add: Novel_All (all novel host genes combined) ----
novel_all_file <- "geneset_Novel_HostGenes_symbols.txt"
if (file.exists(novel_all_file)) {
  novel_all_genes <- scan(novel_all_file, what = "", quiet = TRUE)
  novel_all_genes <- intersect(novel_all_genes, rownames(avg_expr))
  novel_all_mean <- colMeans(avg_expr[novel_all_genes, , drop = FALSE])
  novel_all_z <- (novel_all_mean - mean(novel_all_mean)) / sd(novel_all_mean)
  enrich_z <- rbind(enrich_z, Novel_All = novel_all_z)
}

# ---- Save ----
saveRDS(list(
  enrich_mat = enrich_mat,
  enrich_z = enrich_z,
  gene_sets = gene_sets
), '~/scRNA_analysis/results/enrichment_class_code_novel.rds')

fwrite(as.data.table(enrich_z, keep.rownames = "Layer"),
       '~/scRNA_analysis/results/enrichment_class_code_zscore.csv')

cat('\nClass code enrichment complete.\n')
print(round(enrich_z, 3))
