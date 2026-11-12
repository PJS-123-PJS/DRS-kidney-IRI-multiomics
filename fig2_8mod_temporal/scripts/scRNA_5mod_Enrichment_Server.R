# ============================================================
# Server-side: HC gene-set x cell-type enrichment (GSE274819)
# 2026-07-22 — mirrors scRNA_8mod_enrichment.R EXACTLY for the
# score/Z-score; adds dual permutation tests (global + specificity).
# Sets: 4 new canonical-path HC sets (gene_sets_5mod/) + legacy sets
# (gene_sets/) recomputed under one uniform pipeline.
# Run: nohup Rscript scRNA_HC_enrichment.R > hc_enrich.log 2>&1 &
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(parallel)
})
setwd("~/scRNA_analysis")

cat("Loading Seurat object...\n")
so <- readRDS("results/seurat_merged_annotated.rds")
DefaultAssay(so) <- "RNA"
Idents(so) <- so$celltype
cat(sprintf("Cells: %d, Genes: %d, CellTypes: %d\n",
            ncol(so), nrow(so), length(unique(so$celltype))))

# ---- ID mapping (legacy sets may be Ensembl) --------------------------------
map_data <- read.table("data/ShamA/features.tsv.gz", header = FALSE, sep = "\t",
                       col.names = c("ensembl", "symbol", "type"))
ens2sym <- setNames(map_data$symbol, map_data$ensembl)

load_sets <- function(dir) {
  out <- list()
  for (f in list.files(dir, pattern = "\\.txt$", full.names = TRUE)) {
    nm <- sub("\\.txt$", "", basename(f))
    gs <- trimws(readLines(f))
    gs <- gs[nzchar(gs) & !grepl("^#", gs)]
    if (length(gs) && grepl("^ENSMUSG", gs[1])) gs <- na.omit(ens2sym[gs])
    gs <- intersect(gs, rownames(so))
    if (length(gs) >= 5) {
      out[[nm]] <- gs
      cat(sprintf("  %-22s: %d genes in scRNA\n", nm, length(gs)))
    } else {
      cat(sprintf("  %-22s: SKIPPED (%d in scRNA)\n", nm, length(gs)))
    }
  }
  out
}
cat("Loading HC gene sets (gene_sets_5mod/):\n")
gene_sets <- load_sets("gene_sets_5mod")
# 5-mod run: validated sets only, no legacy noise sets

# ---- score matrix (identical math to original: mean(colMeans(data[gs,ct]))) --
ct_list <- levels(so$celltype)
ct_list <- ct_list[table(so$celltype)[ct_list] > 0]
mat <- GetAssayData(so, assay = "RNA", slot = "data")
ct_idx <- lapply(ct_list, function(ct) WhichCells(so, idents = ct))
names(ct_idx) <- ct_list

score_set <- function(gs) {
  # mean(colMeans(x)) == mean(colSums(x))/k  — exact, faster on sparse
  vapply(ct_idx, function(cells)
    mean(colSums(mat[gs, cells, drop = FALSE])) / length(gs), numeric(1))
}

cat("Scoring", length(gene_sets), "gene sets...\n")
score_mat <- sapply(gene_sets, score_set)   # ct x set
rownames(score_mat) <- ct_list
z_mat <- scale(score_mat)                   # per-set Z across cell types
write.csv(round(z_mat, 4), "results/enrichment_5mod_zscore.csv")
write.csv(round(score_mat, 6), "results/enrichment_5mod_meanexpr.csv")
cat("Z-score matrix written.\n")

# ---- dual permutation test (1,000x) for the 4 HC sets -------------------------
set.seed(7)
BG <- rownames(so)
hc_names <- names(gene_sets)[vapply(gene_sets, length, numeric(1)) >= 30]
NPERM <- 1000
p_glob <- matrix(NA, length(ct_list), length(hc_names),
                 dimnames = list(ct_list, hc_names))
p_spec <- p_glob

for (nm in hc_names) {
  k <- length(gene_sets[[nm]])
  obs <- score_mat[, nm]
  obs_dev <- abs(obs - mean(obs))
  perm_scores <- mclapply(seq_len(NPERM), function(i) {
    score_set(sample(BG, k))
  }, mc.cores = 16)
  perm_mat <- do.call(cbind, perm_scores)          # ct x NPERM
  cnt_g <- numeric(length(ct_list)); cnt_s <- numeric(length(ct_list))
  for (j in seq_len(NPERM)) {
    rs <- perm_mat[, j]
    cnt_g <- cnt_g + (rs >= obs)                              # global high expr
    cnt_s <- cnt_s + (abs(rs - mean(rs)) >= obs_dev)          # ct specificity
  }
  p_glob[, nm] <- (cnt_g + 1) / (NPERM + 1)
  p_spec[, nm] <- (cnt_s + 1) / (NPERM + 1)
  cat(sprintf("  perm done: %s (k=%d)\n", nm, k), file = stderr())
}
write.csv(round(p_glob, 4), "results/enrichment_5mod_permP_global.csv")
write.csv(round(p_spec, 4), "results/enrichment_5mod_permP_specific.csv")
write.csv(round(matrix(p.adjust(p_glob, "BH"), nrow(p_glob), dimnames = dimnames(p_glob)), 4),
          "results/enrichment_5mod_permFDR_global.csv")
write.csv(round(matrix(p.adjust(p_spec, "BH"), nrow(p_spec), dimnames = dimnames(p_spec)), 4),
          "results/enrichment_5mod_permFDR_specific.csv")
cat("== HC ENRICHMENT DONE ==\n")
