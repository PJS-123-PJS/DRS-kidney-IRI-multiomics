# Fig2 Composite v5 — VALIDATED 13-panel assembly (A–M continuous)
# =====================================================================
# 2026-07-23 re-sync per continuous-lettering rule:
#   v4 referenced retired filenames (Ridge_v3, Panorama_v5, J_Nature...) and
#   had letter gaps (no C/F). v5 = validated roster, A–M no gaps:
#   A Panorama_v6 | B Ridge_v5_5mod | C HC_trajectory | D ModVerdict arbitration
#   E UMAP(validated) | F PartialCor_v3_Light | G Enzyme_Timecourse
#   H Writer_Heatmap | I 9Quadrant | J CellType_Enrichment_v6
#   K Composition_Drift_v2 | L CoTargeting_v4_Validated | M StarGenes
# Implemented with patchwork (no magick dependency); letters via tag_levels.
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(png)
})

OUTDIR  <- "../outputs"
FIG3OUT <- "../../fig3_8mod_spatial/output"
SINGOUT <- "../output"            # legacy singular output dir (Fig2M)
COMPDIR <- "../outputs"
dir.create(COMPDIR, showWarnings = FALSE, recursive = TRUE)

panels <- list(
  A = file.path(OUTDIR,  "Fig2A_Panorama_v6.png"),
  B = file.path(OUTDIR,  "Fig2B_Ridge_v5_5mod.png"),
  C = file.path(OUTDIR,  "Fig2C_HC_trajectory.png"),
  D = file.path(FIG3OUT, "ModVerdict_Paired_Arbitration.png"),
  E = file.path(OUTDIR,  "Fig2E_UMAP.png"),
  F = file.path(OUTDIR,  "Fig2F_PartialCor_v3_Light.png"),
  G = file.path(OUTDIR,  "Fig2G_Enzyme_Timecourse.png"),
  H = file.path(OUTDIR,  "Fig2H_Writer_Heatmap.png"),
  I = file.path(OUTDIR,  "Fig2I_9Quadrant.png"),
  J = file.path(OUTDIR,  "Fig2J_CellType_Enrichment_v6_validated.png"),
  K = file.path(OUTDIR,  "Fig2K_Composition_Drift_v2.png"),
  L = file.path(FIG3OUT, "Fig2L_v4_CoTargeting_Validated.png"),
  M = file.path(SINGOUT, "Fig2M_StarGenes.png")
)
missing <- names(panels)[!file.exists(unlist(panels))]
if (length(missing)) stop("missing panel files: ", paste(missing, collapse=", "))

# aspect ratio (height/width) of each png
aspect <- sapply(unlist(panels), function(f) {
  d <- dim(readPNG(f, info = TRUE)); d[1] / d[2]
})

# 4-column design; row heights proportional to panel aspect x span
design <- "
AAAA
BBCC
DDDD
EEFF
GGGG
HHII
JJKK
LLMM
"
row_h <- c(
  4 * aspect["A"],
  2 * max(aspect["B"], aspect["C"]),
  4 * aspect["D"],
  2 * max(aspect["E"], aspect["F"]),
  4 * aspect["G"],
  2 * max(aspect["H"], aspect["I"]),
  2 * max(aspect["J"], aspect["K"]),
  2 * max(aspect["L"], aspect["M"])
)

plots <- lapply(unlist(panels), function(f) {
  ggplot() + annotation_raster(readPNG(f), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_void() + coord_fixed(ratio = 1, expand = FALSE) +
    theme(plot.margin = margin(2, 2, 2, 2))
})
# annotation_raster stretches; preserve per-panel aspect via coord_free is not
# possible in a fixed grid — row heights above carry the aspect correction.

composite <- wrap_plots(plots, design = design, heights = row_h) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = "Arial"),
        plot.tag.position = c(0.005, 0.995))

W_IN <- 10
H_IN <- W_IN * sum(row_h) / 4
ggsave(file.path(COMPDIR, "Fig2_Composite_v5.png"), composite,
       width = W_IN, height = H_IN, dpi = 400, limitsize = FALSE)
ggsave(file.path(COMPDIR, "Fig2_Composite_v5.pdf"), composite,
       width = W_IN, height = H_IN, device = cairo_pdf, limitsize = FALSE)
cat(sprintf("Composite v5 saved (A-M validated roster): %.1f x %.1f in\n", W_IN, H_IN))
