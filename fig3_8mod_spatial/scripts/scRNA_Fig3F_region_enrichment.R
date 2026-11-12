# ============================================================
# Server-side: Fig3F region-stratified m6A-UP sets x cell-type enrichment
# 2026-07-23 — mirrors scRNA_5mod_Enrichment_Server.R EXACTLY
# (score = mean(colSums)/k; Z = scale per set; dual 1,000 permutations).
# Sets: gene_sets_fig3f/ (m6A_UP_CDS / m6A_UP_3UTR / m6A_UP_5UTR)
#       + M6A_UP full set (gene_sets_5mod) as reference row.
# Run: nohup Rscript scRNA_Fig3F_region_enrichment.R > fig3f_enrich.log 2>&1 &
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(parallel)
})
setwd("~/scRNA_analysis")

cat("Loading Seurat object...\n")
so <- readRDS("results/seurat_merged_annotated.rds")
DefaultAssay(so) <- "RNA"; Idents(so) <- so$celltype
cat(sprintf("Cells: %d, Genes: %d, CellTypes: %d\n",
            ncol(so), nrow(so), length(unique(so$celltype))))

load_sets <- function(dir) {
  out <- list()
  for (f in list.files(dir, pattern = "\\.txt$", full.names = TRUE)) {
    nm <- sub("\\.txt$", "", basename(f))
    gs <- trimws(readLines(f))
    gs <- gs[nzchar(gs) & !grepl("^#", gs)]
    gs <- intersect(gs, rownames(so))
    if (length(gs) >= 5) {
      out[[nm]] <- gs
      cat(sprintf("  %-16s: %d genes in scRNA\n", nm, length(gs)))
    } else cat(sprintf("  %-16s: SKIPPED (%d in scRNA)\n", nm, length(gs)))
  }
  out
}
cat("Loading Fig3F region sets:\n")
gene_sets <- load_sets("gene_sets_fig3f")
# reference row: full validated M6A_UP
ref <- trimws(readLines("gene_sets_5mod/M6A_UP.txt"))
ref <- intersect(ref[nzchar(ref) & !grepl("^#", ref)], rownames(so))
gene_sets[["M6A_UP_full"]] <- ref
cat(sprintf("  %-16s: %d genes in scRNA (reference)\n", "M6A_UP_full", length(ref)))

ct_list <- levels(so$celltype); ct_list <- ct_list[table(so$celltype)[ct_list] > 0]
mat <- GetAssayData(so, assay = "RNA", slot = "data")
ct_idx <- lapply(ct_list, function(ct) WhichCells(so, idents = ct)); names(ct_idx) <- ct_list

score_set <- function(gs) {
  vapply(ct_idx, function(cells)
    mean(colSums(mat[gs, cells, drop = FALSE])) / length(gs), numeric(1))
}

cat("Scoring", length(gene_sets), "gene sets...\n")
score_mat <- sapply(gene_sets, score_set)
rownames(score_mat) <- ct_list
z_mat <- scale(score_mat)
write.csv(round(z_mat, 4), "results/enrichment_fig3f_zscore.csv")
write.csv(round(score_mat, 6), "results/enrichment_fig3f_meanexpr.csv")
cat("Z-score matrix written.\n")

set.seed(7)
BG <- rownames(so)
hc_names <- names(gene_sets)[vapply(gene_sets, length, numeric(1)) >= 30]
NPERM <- 1000
p_glob <- matrix(NA, length(ct_list), length(hc_names), dimnames = list(ct_list, hc_names))
p_spec <- p_glob
for (nm in hc_names) {
  k <- length(gene_sets[[nm]])
  obs <- score_mat[, nm]; obs_dev <- abs(obs - mean(obs))
  perm_scores <- mclapply(seq_len(NPERM), function(i) score_set(sample(BG, k)), mc.cores = 16)
  perm_mat <- do.call(cbind, perm_scores)
  cnt_g <- numeric(length(ct_list)); cnt_s <- numeric(length(ct_list))
  for (j in seq_len(NPERM)) {
    rs <- perm_mat[, j]
    cnt_g <- cnt_g + (rs >= obs)
    cnt_s <- cnt_s + (abs(rs - mean(rs)) >= obs_dev)
  }
  p_glob[, nm] <- (cnt_g + 1) / (NPERM + 1)
  p_spec[, nm] <- (cnt_s + 1) / (NPERM + 1)
  cat(sprintf("  perm done: %s (k=%d)\n", nm, k), file = stderr())
}
write.csv(round(p_glob, 4), "results/enrichment_fig3f_permP_global.csv")
write.csv(round(p_spec, 4), "results/enrichment_fig3f_permP_specific.csv")
write.csv(round(matrix(p.adjust(p_glob, "BH"), nrow(p_glob), dimnames = dimnames(p_glob)), 4),
          "results/enrichment_fig3f_permFDR_global.csv")
write.csv(round(matrix(p.adjust(p_spec, "BH"), nrow(p_spec), dimnames = dimnames(p_spec)), 4),
          "results/enrichment_fig3f_permFDR_specific.csv")
cat("== FIG3F ENRICHMENT DONE ==\n")
