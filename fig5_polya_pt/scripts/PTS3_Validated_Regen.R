# =============================================================================
# Regenerate PTS3 layer tables on the VALIDATED multi-layer basis
# 2026-07-22 — old PTS3_Layer_{Matrix,Coverage,Enrichment}.csv were built from
# the retired wave layers (m6A↓/m5C↓/PSU↑/Um↑). Rebuilt from
# Table_MultiLayer_Genes_validated.csv (9 validated layers).
# Products (fig7_multiomics/output/PTS3/):
#   PTS3_Layer_Matrix.csv       Gene x 9 layers (0/1), PT-S3 markers only
#   PTS3_Layer_Coverage.csv     per-layer: global N, PT-S3 hits, pct, genes
#   PTS3_Layer_Enrichment.csv   per-layer: global% vs PT-S3% + fold
#   PTS3_vs_Global_LayerDist.csv  N_Layers distribution PT-S3 vs global
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

ml <- fread("../../fig7_multiomics/data/Table_MultiLayer_Genes_validated.csv")
out_dir <- "../../fig7_multiomics/output/PTS3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PT_S3 <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1",
           "Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc13a1","Slc22a7","Slc13a3","Lrp2","Cubn","Miox")
LAYERS <- c("m6A↑","m5C↓","PSU↑","inosine↑","Nm↑","Nm↓","AS_diff","PolyA_diff","Novel_Tx")

pt <- ml[Gene %in% PT_S3]
cat(sprintf("PT-S3 markers found in validated table: %d / %d\n", nrow(pt), length(PT_S3)))

# ---- matrix -----------------------------------------------------------------
mat <- as.data.table(sapply(LAYERS, function(ln)
  as.integer(vapply(strsplit(pt$Layers, "; ", fixed = TRUE), function(x) ln %in% x, logical(1)))))
mat <- cbind(Gene = pt$Gene, mat)
fwrite(mat, file.path(out_dir, "PTS3_Layer_Matrix.csv"))

# ---- coverage + enrichment ---------------------------------------------------
total_bg <- nrow(ml)   # genes with >=1 validated layer
cov_list <- list()
for (ln in LAYERS) {
  in_layer <- ml[vapply(strsplit(Layers, "; ", fixed = TRUE), function(x) ln %in% x, logical(1))]
  pt_hits <- intersect(in_layer$Gene, PT_S3)
  cov_list[[ln]] <- data.table(
    Layer = ln,
    Global_N = nrow(in_layer),
    PTS3_N = length(pt_hits),
    PTS3_Pct = round(100 * length(pt_hits) / length(PT_S3), 2),
    PTS3_Genes = paste(pt_hits, collapse = ";"),
    Global_Pct = round(100 * nrow(in_layer) / total_bg, 2),
    Fold_Enrichment = round((length(pt_hits) / length(PT_S3)) / (nrow(in_layer) / total_bg), 2)
  )
}
cov_dt <- rbindlist(cov_list)
fwrite(cov_dt[, .(Layer, Global_N, PTS3_N, PTS3_Pct, PTS3_Genes)],
       file.path(out_dir, "PTS3_Layer_Coverage.csv"))
fwrite(cov_dt[, .(Layer, Global_N, PTS3_N, Global_Pct, PTS3_Pct, Fold_Enrichment)],
       file.path(out_dir, "PTS3_Layer_Enrichment.csv"))
print(cov_dt[, .(Layer, Global_N, PTS3_N, PTS3_Pct, Fold_Enrichment)])

# ---- PT-S3 vs global layer-count distribution --------------------------------
g_dist <- ml[, .N, by = N_Layers][, `:=`(pct = 100 * N / sum(N), group = "Global")][, N := NULL]
p_dist <- ml[Gene %in% PT_S3, .N, by = N_Layers][, `:=`(pct = 100 * N / sum(N), group = "PT-S3")][, N := NULL]
dist_dt <- rbind(g_dist, p_dist)[order(N_Layers, group)]
fwrite(dist_dt, file.path(out_dir, "PTS3_vs_Global_LayerDist.csv"))
cat("PTS3 VALIDATED REGEN DONE\n")
