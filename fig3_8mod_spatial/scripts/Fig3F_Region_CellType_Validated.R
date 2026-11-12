# Fig3 Panel F v2 — region-stratified m6A-UP x cell-type enrichment (validated)
# ============================================================================
# 2026-07-23 — replaces wave-era Fig3F_WaveRegion_CellType (archived
# _archive_pre_v7.0). Sets built by build_Fig3F_region_sets.R (validated
# M6A_UP stratified by canonical Gene_Struc of HC UP sites), scored on
# GSE274819 by scRNA_Fig3F_region_enrichment.R (mean(colSums)/k, Z=scale,
# 1,000 dual permutations — same math as Fig2J pipeline).
# Reads server CSVs from data/fig3f_region_sets/. HARD FAIL if missing.
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

DD <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/fig3f_region_sets"
OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
for (f in c("enrichment_fig3f_zscore.csv", "enrichment_fig3f_permFDR_specific.csv"))
  if (!file.exists(file.path(DD, f))) stop("HARD FAIL: missing ", f)

z <- fread(file.path(DD, "enrichment_fig3f_zscore.csv"))
fd <- fread(file.path(DD, "enrichment_fig3f_permFDR_specific.csv"))
setnames(z, 1, "celltype"); setnames(fd, 1, "celltype")

zl <- melt(z, id.vars = "celltype", variable.name = "set", value.name = "Z")
fl <- melt(fd, id.vars = "celltype", variable.name = "set", value.name = "FDR")
dt <- merge(zl, fl, by = c("celltype", "set"))
SET_ORDER <- c("M6A_UP_full", "m6A_UP_CDS", "m6A_UP_3UTR", "m6A_UP_5UTR")
SET_LAB <- c(M6A_UP_full = "m6A↑ all (ref)", m6A_UP_CDS = "m6A↑-CDS",
             m6A_UP_3UTR = "m6A↑-3'UTR", m6A_UP_5UTR = "m6A↑-5'UTR")
dt[, set := factor(set, levels = SET_ORDER, labels = SET_LAB[SET_ORDER])]
ct_order <- z$celltype[order(-z$M6A_UP_full)]
dt[, celltype := factor(celltype, levels = ct_order)]
dt[, star := ifelse(FDR < 0.01, "**", ifelse(FDR < 0.05, "*", ""))]
zcap <- max(abs(dt$Z), na.rm = TRUE)

p <- ggplot(dt, aes(x = celltype, y = set, fill = Z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f%s", Z, star)), size = 2.1, family = "Arial") +
  scale_fill_gradient2(low = "#3B6EA5", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-zcap, zcap), name = "Z") +
  labs(x = NULL, y = NULL,
       title = "f  Region-stratified m6A-gain targets × cell-type enrichment",
       subtitle = "Validated m6A↑ set split by transcript region of HC UP sites (canonical Gene_Struc); GSE274819, 1,000 dual permutations",
       caption = paste0("Sets: build_Fig3F_region_sets.R (CDS 852 / 3'UTR 753 / 5'UTR 83 genes; 89 mixed-region excluded).\n",
                        "* FDR<0.05, ** FDR<0.01 (specificity permutation). Scoring mirrors Fig2J pipeline (mean(colSums)/k, scale).")) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "grey30"),
        plot.caption = element_text(size = 5.5, color = "grey50", hjust = 0, lineheight = 1.15),
        legend.position = "right", legend.key.height = unit(4, "mm"))

ggsave(file.path(OUT_DIR, "Fig3F_Region_CellType_Validated.pdf"), p,
       width = 180/25.4, height = 84/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3F_Region_CellType_Validated.png"), p,
       width = 180/25.4, height = 84/25.4, dpi = 600)
cat("Fig3F v2 saved\n")

# ---- transcript-visible summary for manuscript rewrite ----
cat("\n=== top |Z| per set ===\n")
for (s in SET_ORDER) {
  d <- dt[set == SET_LAB[s]][order(-abs(Z))][1:3]
  cat(sprintf("%-14s: %s\n", SET_LAB[s],
              paste(sprintf("%s Z=%+.2f%s", d$celltype, d$Z, d$star), collapse = "; ")))
}
