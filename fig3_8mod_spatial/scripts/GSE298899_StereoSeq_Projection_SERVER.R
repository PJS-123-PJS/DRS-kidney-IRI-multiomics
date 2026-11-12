# =============================================================================
# GSE298899 Stereo-seq projection (sham n=3 vs bIRI n=3, 14d) — SERVER version
# 7 validated direction sets x author compartments (figure_idents incl. FR-PT).
# Module score = mean(SCT set) - mean(SCT all genes) per bin; Z pooled across
# all bins; stats at SAMPLE level (n=3/group, Welch t) — animal-level spirit.
# Output: ~/GSE298899/stereoseq_projection_{sample_means,stats}.csv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(Matrix) })

setwd("/home/pjs/GSE298899")
gs_dir <- "/home/pjs/scRNA_analysis/gene_sets_5mod"
mod_files <- c(m6A_UP = "M6A_UP.txt", m6A_DOWN_resid = "M6A_DOWN.txt",
               m5C_DOWN = "M5C_DOWN.txt", PSU_UP = "PSU_UP.txt",
               Inosine_UP = "INOSINE_UP.txt", Nm_UP = "NM_UP.txt",
               Nm_DOWN = "NM_DOWN.txt", Um_UP = "UM_UP.txt")

md  <- fread("GSE298899_metadata.csv.gz")
bc  <- trimws(fread("GSE298899_barcodes.tsv.gz", header = FALSE)[[1]])
ft  <- trimws(fread("GSE298899_sct_features.tsv.gz", header = FALSE)[[1]])
                        # sct_matrix row names: 22277 detected genes (NOT raw_features 27770)
cat(sprintf("meta=%d barcodes=%d features=%d\n", nrow(md), length(bc), length(ft)))

cat("reading mtx via fread (few minutes) ...\n")
tri <- fread("GSE298899_sct_matrix.mtx.gz", skip = 3,
             col.names = c("i", "j", "x"))
hdr <- readLines(gzfile("GSE298899_sct_matrix.mtx.gz"), n = 5)
dims <- as.integer(strsplit(hdr[!grepl("^%", hdr)][1], " ")[[1]])
cat("mtx dims:", dims, "\n")          # rows cols nnz  (rows=features, cols=barcodes per mtx convention)
M <- sparseMatrix(i = tri$i, j = tri$j, x = tri$x, dims = dims[1:2])
rm(tri); gc(verbose = FALSE)
if (nrow(M) == length(ft) && ncol(M) == length(bc)) {
  M <- t(M)                              # to bins x genes
} else stopifnot(nrow(M) == length(bc))
stopifnot(ncol(M) == length(ft))
cat("matrix oriented:", nrow(M), "bins x", ncol(M), "genes\n")

md <- as.data.frame(md); rownames(md) <- md$barcodes
md <- md[bc, ]; stopifnot(all(rownames(md) == bc))

gmean <- Matrix::rowMeans(M)
gidx  <- setNames(seq_along(ft), ft)

score_one <- function(genes) {
  g <- intersect(trimws(genes), names(gidx))
  if (length(g) < 5) return(NULL)
  list(sc = Matrix::rowMeans(M[, gidx[g], drop = FALSE]) - gmean, n = length(g))
}

S <- data.frame(row.names = bc)
for (m in names(mod_files)) {
  r <- score_one(readLines(file.path(gs_dir, mod_files[[m]])))
  if (is.null(r)) { cat(m, "skipped\n"); next }
  S[[m]] <- as.numeric(scale(r$sc))      # pooled Z across all bins
  cat(sprintf("  %s: %d genes\n", m, r$n))
}
res <- cbind(md[, c("barcodes", "sample", "figure_idents")], S)
res$group <- ifelse(grepl("^sham", res$sample), "Sham", "IRI_14d")

long <- melt(as.data.table(res), id.vars = c("barcodes", "sample", "group", "figure_idents"),
             variable.name = "set", value.name = "z")
summ <- long[, .(z = mean(z)), by = .(set, figure_idents, group, sample)]
fwrite(summ, "stereoseq_projection_sample_means.csv")

pv <- summ[, {
  a <- z[group == "Sham"]; b <- z[group == "IRI_14d"]
  if (length(a) >= 2 && length(b) >= 2)
    .(Sham_meanZ = mean(a), IRI14d_meanZ = mean(b), deltaZ = mean(b) - mean(a),
      welch_p = tryCatch(t.test(b, a)$p.value, error = function(e) NA_real_))
  }, by = .(set, figure_idents)]
setorder(pv, set, deltaZ)
fwrite(pv, "stereoseq_projection_stats.csv")

key <- c("PT-S1/S2","PT-S3","FR-PT","TAL","CNT","Injured tubule")
cat("\n=== key compartments (Sham vs IRI-14d mean Z) ===\n")
print(as.data.frame(pv[figure_idents %in% key]))
cat("\n=== DONE ===\n")
