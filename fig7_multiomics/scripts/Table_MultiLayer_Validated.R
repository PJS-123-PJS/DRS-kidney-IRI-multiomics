# =============================================================================
# Multi-layer table, VALIDATED layers (replaces archived wave-based generator)
# 2026-07-22 — layers = validated direction sets + AS/PolyA/Novel:
#   m6A↑ (DRACH+ HC50 2,829) | m5C↓ (gate+ HC50 1,140) | PSU↑ (HC50 106)
#   inosine↑ (HC50 246) | Nm↑ (union 331) | Nm↓ (union 289)
#   AS_diff | PolyA_diff | Novel_Tx  (unchanged existing sets)
# Products: Table_MultiLayer_Genes_validated.csv / Table_Layer_Sizes_validated.csv
#           Table_Layer_Jaccard_validated.csv  (fig7_multiomics/data)
# Star genes = genes in >=3 layers.
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

vdir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/upload_5mod_validated"
f2d <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig2_8mod_temporal/data"
f4d <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig4_as_landscape/data"
out_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig7_multiomics/data"

rd <- function(f) { g <- readLines(f); unique(g[!g %in% c("-", "", "NA")]) }

all_layers <- list(
  "m6A↑"      = rd(file.path(vdir, "M6A_UP.txt")),
  "m5C↓"      = rd(file.path(vdir, "M5C_DOWN.txt")),
  "PSU↑"      = rd(file.path(vdir, "PSU_UP.txt")),
  "inosine↑"  = rd(file.path(vdir, "INOSINE_UP.txt")),
  "Nm↑"       = rd(file.path(vdir, "NM_UP.txt")),
  "Nm↓"       = rd(file.path(vdir, "NM_DOWN.txt")),
  "AS_diff"   = rd(file.path(f4d, "geneset_AS_diff.txt")),
  "PolyA_diff"= rd(file.path(f4d, "geneset_PolyA_diff.txt")),
  "Novel_Tx"  = rd(file.path(out_dir, "geneset_Novel_Tx_converted.txt"))
)
layer_sizes <- data.table(Layer = names(all_layers),
                          N_Genes = vapply(all_layers, length, integer(1)))
print(layer_sizes)

gene_summary <- rbindlist(lapply(names(all_layers), function(ln)
  data.table(Gene = all_layers[[ln]], Layer = ln)))
gene_hits <- gene_summary[, .(N_Layers = uniqueN(Layer),
                              Layers = paste(sort(unique(Layer)), collapse = "; ")), by = Gene]
setorder(gene_hits, -N_Layers, Gene)

ln <- names(all_layers)
jaccard_mat <- matrix(NA_real_, length(ln), length(ln), dimnames = list(ln, ln))
for (i in seq_along(ln)) for (j in seq_along(ln))
  jaccard_mat[i, j] <- length(intersect(all_layers[[i]], all_layers[[j]])) /
                       length(union(all_layers[[i]], all_layers[[j]]))

fwrite(gene_hits, file.path(out_dir, "Table_MultiLayer_Genes_validated.csv"))
fwrite(layer_sizes, file.path(out_dir, "Table_Layer_Sizes_validated.csv"))
fwrite(as.data.table(jaccard_mat, keep.rownames = "Layer"),
       file.path(out_dir, "Table_Layer_Jaccard_validated.csv"))

star <- gene_hits[N_Layers >= 3]
cat(sprintf("\ntotal genes in >=1 layer: %s\n", format(nrow(gene_hits), big.mark = ",")))
cat(sprintf("star genes (>=3 layers): %d (%.1f%%)\n", nrow(star), 100 * nrow(star) / nrow(gene_hits)))
cat("top 25 star genes:\n")
print(head(star, 25))
# legacy comparison
old <- fread(file.path(out_dir, "Table_MultiLayer_Genes.csv"))
old_star <- old[N_Layers >= 3]$Gene
cat(sprintf("\nlegacy star genes: %d | validated star genes: %d | overlap: %d\n",
            length(old_star), nrow(star), length(intersect(old_star, star$Gene))))
cat("MULTILAYER VALIDATED DONE\n")
