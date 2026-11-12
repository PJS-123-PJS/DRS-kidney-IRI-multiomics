# Fig7D v3 — Optimized: PT-S3 Multi-Layer Gene Distribution
# ============================================================================
# vision_analyzer: 7/10 → v2 8.5/10 → v3 target 9.5+
# v3 fixes:
#   1. Layer annotation: full "PolyA" display, no truncated "F"
#   2. Color legend added (N Layers: 1=PolyA only, 2, 3=Ggt1)
#   3. Title margin fix (plot.margin top=6)
#   4. Wider right panel + better text spacing
#   5. Coord ylim extended to 4.2 for label room
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("Fig6_config.R")

theme_nature <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1, "pt"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    strip.text = element_text(size = 6.5, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 6, color = "grey30"),
    panel.grid = element_blank()
  )

# ===========================================================================
# Panel D: PT-S3 Multi-Layer Gene Distribution (v3 optimized)
# ===========================================================================
cat("Panel D (v3 optimized): PT-S3 Multi-Layer Distribution...\n")
multi_genes <- fread("../data/Table_MultiLayer_Genes_validated.csv")

pt_s3_markers <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1",
                   "Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc13a1","Slc22a7","Slc13a3","Lrp2","Cubn","Miox")

multi_genes[, PT := ifelse(Gene %in% pt_s3_markers, "PT (10/18 found)", "Non-PT (19,375 genes)")]
multi_genes[, LayerCount := factor(N_Layers, levels = 1:3)]

# ── D1: Layer distribution (left panel) ──
layer_dist <- multi_genes[, .N, by = .(PT, LayerCount)]
layer_dist[, Pct := N / sum(N) * 100, by = PT]

p_d1 <- ggplot(layer_dist, aes(x = LayerCount, y = Pct, fill = PT)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_text(aes(label = paste0(sprintf("%.1f", Pct), "%"),
                group = PT),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.0, family = "Arial", color = "grey30") +
  scale_fill_manual(
    values = c("PT (10/18 found)" = "#D6404C", "Non-PT (19,375 genes)" = "grey75"),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of Omics Layers", y = "% of Genes",
       title = "Multi-Layer Distribution: PT vs Non-PT",
       subtitle = "PT-S3: 10/18 markers found. 90% have only 1 layer (PolyA). Non-PT: 2.2% multi-layer.") +
  theme_nature +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))

# ── D2: PT gene composition (right panel) ──
pt_multi <- multi_genes[Gene %in% pt_s3_markers][order(-N_Layers)]
pt_multi[, Gene := factor(Gene, levels = rev(pt_multi$Gene))]

# Clean layer display: replace "PolyA_diff" with "PolyA", separate with " | "
pt_multi[, layer_label := gsub("PolyA_diff", "PolyA", Layers)]
pt_multi[, layer_label := gsub("; ", " | ", layer_label)]

p_d2 <- ggplot(pt_multi, aes(x = Gene, y = N_Layers)) +
  geom_col(aes(fill = factor(N_Layers)), width = 0.55) +
  # Gene name inside bar (white)
  geom_text(aes(label = Gene, y = 0.12),
            hjust = 0, size = 3.2, family = "Arial", color = "white", fontface = "bold") +
  # Layer annotation outside bar
  geom_text(aes(label = layer_label, y = N_Layers + 0.15),
            hjust = 0, size = 2.8, family = "Arial", color = "grey40") +
  coord_flip(ylim = c(0, 4.2)) +
  scale_fill_manual(
    values = c("1" = "#D9704E", "2" = "#C0403C", "3" = "#9B1B24"),
    name = "N Layers",
    labels = c("1 (PolyA only)", "2", "3 (Ggt1)")
  ) +
  # Ggt1 row highlight
  annotate("rect", ymin = 0, ymax = 4.2,
           xmin = which(levels(pt_multi$Gene) == "Ggt1") - 0.4,
           xmax = which(levels(pt_multi$Gene) == "Ggt1") + 0.4,
           fill = "#D6404C", alpha = 0.06) +
  labs(x = NULL, y = "Number of Omics Layers",
       title = "PT-S3 Marker Layer Composition",
       subtitle = "Ggt1 (shaded) = sole multi-layer PT gene (m6A + AS + PolyA). 9 others = PolyA only.") +
  theme_nature +
  theme(
    legend.position = c(0.85, 0.18),
    legend.key.size = unit(2.5, "mm"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    plot.margin = margin(6, 8, 4, 6)
  )

# ── Combine & Save ──
p_d <- wrap_plots(p_d1, p_d2, ncol = 2, widths = c(0.45, 0.55))

ggsave(file.path(OUTDIR, "Fig7D_PTS3_MultiLayer.pdf"), p_d,
       width = 8.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig7D_PTS3_MultiLayer.png"), p_d,
       width = 8.5, height = 3.5, dpi = 600)
cat("  Fig7D v3 Done\n")
