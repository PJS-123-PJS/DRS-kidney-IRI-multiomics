# Fig7 Composite v8.0 — PT-S3 Centered Multi-Omics (Nature theme)
# ============================================================================
# 5 panels: A(PT-S3 Circos) | B(Temporal Decoupling) | C(Hub Network) |
#            D(SAM-m6A Metabolome) | E(DRS Single-Platform Advantage)
# ============================================================================

library(png); library(grid)

source("Fig6_config.R")

# ── Section colors ──
sec_col <- c(
  "Overview"     = "#D6404C",
  "Temporal"     = "#E8753A",
  "Network"      = "#3D7EA6",
  "Metabolome"   = "#2E7D6F",
  "Advantage"    = "#7B5294"
)

# ── Panel definitions ──
panels <- list(
  list(id = "a", file = "Fig7A_4Omics_Circos_Journal",       section = "Overview",    row = 1, col = 1, cspan = 2),
  list(id = "b", file = "Fig7E_TimeDecoupling_Core",         section = "Temporal",    row = 2, col = 1, cspan = 1),
  list(id = "c", file = "Fig7F_HubGene_Network",             section = "Network",     row = 2, col = 2, cspan = 1),
  list(id = "d", file = "Fig7J_Metabolome_Epitranscriptome", section = "Metabolome",  row = 3, col = 1, cspan = 1),
  list(id = "e", file = "Fig7E_DRS_Advantage",               section = "Advantage",   row = 3, col = 2, cspan = 1)
)

# ── Layout (Nature full-page: 180mm × 220mm) ──
total_w <- 7.1; total_h <- 9.0
left_margin <- 0.55; right_margin <- 0.15
top_margin <- 0.30; bottom_margin <- 0.15
gap_x <- 0.10; gap_y <- 0.08
ncols <- 2
col_w <- (total_w - left_margin - right_margin - gap_x) / ncols
row_h <- c(2.2, 3.0, 3.0)  # Row 1: wide circos, Rows 2-3: paired panels

compute_rect <- function(row, col, cspan = 1) {
  y_above <- sum(row_h[1:(row-1)]) + (row - 1) * gap_y
  y_top <- top_margin + y_above
  y_bot <- y_top + row_h[row]
  x_left <- left_margin + (col - 1) * (col_w + gap_x)
  w <- col_w * cspan + gap_x * (cspan - 1)
  c(x = x_left, y = total_h - y_top, w = w, h = row_h[row])
}

sections <- list(
  list(label = "PT-S3 Multi-Omics Overview", rows = 1:1, color = sec_col["Overview"]),
  list(label = "Temporal + Network",          rows = 2:2, color = sec_col["Temporal"]),
  list(label = "Metabolome + DRS Advantage",  rows = 3:3, color = sec_col["Advantage"])
)

# ── Read panels ──
cat("Loading panels...\n")
panel_grobs <- list()
for (p in panels) {
  fpath <- file.path(OUTDIR, paste0(p$file, ".png"))
  if (!file.exists(fpath)) {
    cat(sprintf("  MISSING: %s\n", fpath))
    next
  }
  img <- readPNG(fpath)
  panel_grobs[[p$id]] <- img
  cat(sprintf("  %s: %dx%d\n", p$id, dim(img)[2], dim(img)[1]))
}

# If Panel E (DRS Advantage) is missing, generate placeholder
if (!"e" %in% names(panel_grobs)) {
  cat("  Creating Panel E placeholder...\n")
  # Generate DRS Advantage panel
  png(file.path(OUTDIR, "Fig7E_DRS_Advantage.png"), width = 800, height = 600, res = 150)
  par(mar = c(4, 8, 3, 2), family = "Arial")

  layer_sizes <- c("DRS Expression" = 19385, "m6A" = 530, "m5C" = 7,
                   "PSU" = 10, "Um" = 118, "AS" = 124, "PolyA" = 6955,
                   "Novel Tx" = 12083, "ATAC" = 3411)
  layer_platform <- c("DRS Expression" = "DRS", "m6A" = "DRS", "m5C" = "DRS",
                      "PSU" = "DRS", "Um" = "DRS", "AS" = "DRS", "PolyA" = "DRS",
                      "Novel Tx" = "DRS", "ATAC" = "ATAC-seq")

  bp <- barplot(rev(log10(layer_sizes)), horiz = TRUE, las = 1,
                col = ifelse(rev(layer_platform) == "DRS", "#D6404C", "#3D7EA6"),
                border = NA, xlab = "log10(Gene Count)",
                main = "DRS Single-Platform Advantage\nGenes Captured per Omics Layer")
  text(rev(log10(layer_sizes)) + 0.1, bp, labels = rev(layer_sizes),
       cex = 0.8, adj = 0)
  legend("bottomright", legend = c("DRS (single run)", "ATAC-seq (separate)"),
         fill = c("#D6404C", "#3D7EA6"), cex = 0.8, border = NA)
  dev.off()

  img <- readPNG(file.path(OUTDIR, "Fig7E_DRS_Advantage.png"))
  panel_grobs[["e"]] <- img
  cat("  Panel E created.\n")
}

# ── Draw composite ──
draw_composite <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  # Section labels
  for (sec in sections) {
    r_start <- sec$rows[1]; r_end <- sec$rows[length(sec$rows)]
    rect_s <- compute_rect(r_start, 1, 2)
    rect_e <- compute_rect(r_end, 1, 2)
    y_mid <- (rect_s["y"] + rect_e["y"] + rect_e["h"]) / 2

    pushViewport(viewport(x = unit(0.22, "inches"), y = unit(y_mid, "inches"),
                          width = unit(0.44, "inches"), height = unit(1, "npc"),
                          just = "center"))
    grid.text(sec$label, rot = 90, gp = gpar(fontsize = 5.5, fontface = "bold",
              col = sec$color, fontfamily = "Arial"))
    popViewport()
  }

  # Panels
  for (p in panels) {
    if (!p$id %in% names(panel_grobs)) next
    rect <- compute_rect(p$row, p$col, p$cspan)
    img <- panel_grobs[[p$id]]

    pushViewport(viewport(x = unit(rect["x"], "inches"), y = unit(rect["y"], "inches"),
                          width = unit(rect["w"], "inches"), height = unit(rect["h"], "inches"),
                          just = c(0, 0)))
    grid.rect(gp = gpar(fill = NA, col = "grey85", lwd = 0.3))
    grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
    grid.text(p$id, x = unit(0.025, "npc"), y = unit(0.975, "npc"),
              just = c(0, 1), gp = gpar(fontsize = 8.5, fontface = "bold",
              col = sec_col[p$section], fontfamily = "Arial"))
    popViewport()
  }

  # Title
  grid.text("PT-S3 Multi-Omics Convergence: DRS Captures Expression, Modifications, Splicing & PolyA in One Run",
    x = unit(total_w / 2, "inches"), y = unit(0.08, "inches"),
    just = c(0.5, 0), gp = gpar(fontsize = 8, fontface = "bold", col = "grey20", fontfamily = "Arial"))
}

# ── Export ──
cat("\nExporting...\n")
cairo_pdf(file.path(OUTDIR, "Fig7_Composite_PTS3.pdf"), width = total_w, height = total_h, family = "Arial")
draw_composite(); dev.off()
cat("  PDF saved\n")

png(file.path(OUTDIR, "Fig7_Composite_PTS3.png"), width = total_w, height = total_h, units = "in", res = 600)
draw_composite(); dev.off()
cat("  PNG saved\n")

cat(sprintf("\nFig7 PT-S3 composite — %.1f × %.1f inches — 5 panels (A-E)\n", total_w, total_h))
