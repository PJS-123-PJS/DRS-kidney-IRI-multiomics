# Fig7B v3 — Optimized: Layer Size + Jaccard Overlap
# ============================================================================
# v3 fixes (vision_analyzer 7.5→9.0):
#   1. Per-layer colors (NOT generic "Modification" red) — each mod gets its own 5-mod color
#   2. Small-layer callouts: arrow+hollow for m5C(7), PSU(10)
#   3. Bigger heatmap text: 3.0pt, better legend
#   4. Legend moved to top
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
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

# Per-layer colors (from 5-mod palette + extensions)
LAYER_COLORS <- c(
  "m6A↓"        = "#D6404C",  # m6A red
  "m5C↓"        = "#E8753A",  # m5C orange
  "PSU↑"        = "#3D7EA6",  # PSU blue
  "Um↑"         = "#C49A3C",  # Um gold
  "AS_diff"     = "#5B8FA8",  # AS light-blue
  "PolyA_diff"  = "#2E7D6F",  # PolyA green
  "Novel_Tx"    = "grey60"
)

# Ordered category groups (for subtle grouping in bar chart)
layer_categories <- data.table(
  Layer = c("m6A↓", "m5C↓", "PSU↑", "Um↑", "AS_diff", "PolyA_diff", "Novel_Tx"),
  Category = c("Modification", "Modification", "Modification", "Modification",
               "AS", "PolyA", "Novel Transcript")
)

# ===========================================================================
# Panel B: Layer Size + Jaccard Overlap (v3 optimized)
# ===========================================================================
cat("Panel B (v3 optimized): Layer Size + Jaccard Overlap...\n")
layer_sizes <- fread("../data/Table_Layer_Sizes.csv")
layer_jaccard <- fread("../data/Table_Layer_Jaccard.csv")

# ── B1: Log-scale bar chart with per-layer colors ──
layer_order <- c("m6A↓", "m5C↓", "PSU↑", "Um↑", "AS_diff", "PolyA_diff", "Novel_Tx")
layer_sizes[, Layer := factor(Layer, levels = rev(layer_order))]

# Build labels
total_genes_all <- sum(layer_sizes$N_Genes)
layer_sizes[, pct := N_Genes / total_genes_all * 100]
layer_sizes[, bar_label := paste0(comma(N_Genes), " (", sprintf("%.1f", pct), "%)")]
layer_sizes[, is_small := N_Genes < 15]  # m5C (7) and PSU (10)

# Callout data for small layers
small_layers <- layer_sizes[is_small == TRUE]
small_layers[, callout_x := 15]  # fixed position for callout arrow end

p_b1 <- ggplot(layer_sizes, aes(x = N_Genes, y = Layer)) +
  geom_col(aes(fill = Layer), width = 0.55) +
  # Labels for bars >= 15 genes (inside bar area on log scale)
  geom_text(data = layer_sizes[is_small == FALSE],
            aes(label = bar_label, x = N_Genes * 1.08),
            hjust = 0, size = 2.8, family = "Arial", color = "grey30") +
  # Arrow + label for very small bars (m5C, PSU)
  geom_segment(data = small_layers,
               aes(x = N_Genes, xend = 15, y = Layer, yend = Layer),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = "grey30", linewidth = 0.7) +
  geom_text(data = small_layers,
            aes(label = bar_label, x = 18),
            hjust = 0, size = 3.2, family = "Arial", color = "grey20", fontface = "bold") +
  scale_fill_manual(values = LAYER_COLORS, guide = "none") +
  scale_x_log10(
    labels = comma,
    breaks = c(1, 10, 100, 1000, 10000, 20000),
    expand = expansion(mult = c(0.05, 0.45))
  ) +
  annotation_logticks(sides = "b", linewidth = 0.2, alpha = 0.4,
                      short = unit(0.5, "mm"), mid = unit(0.8, "mm"),
                      long = unit(1.2, "mm")) +
  # Category grouping brackets (subtle, right-aligned)
  annotate("text", x = 20000, y = 6.5, label = "Modification",
           size = 2.5, family = "Arial", color = "#D6404C", fontface = "italic", hjust = 1) +
  annotate("text", x = 20000, y = 2.5, label = "AS",
           size = 2.5, family = "Arial", color = "#5B8FA8", fontface = "italic", hjust = 1) +
  annotate("text", x = 20000, y = 1.5, label = "PolyA",
           size = 2.5, family = "Arial", color = "#2E7D6F", fontface = "italic", hjust = 1) +
  labs(x = "Gene Count (log10 scale)", y = NULL,
       title = "Omics Layer Gene Counts",
       subtitle = sprintf("Total: %s genes across 7 layers. Arrows highlight m5C(7) and PSU(10) — sparse but functionally validated.",
                          comma(total_genes_all))) +
  theme_nature +
  theme(
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.2),
    plot.margin = margin(4, 8, 8, 4)
  )

# ── B2: Jaccard heatmap (6×6, clustered, larger text) ──
jac_mat <- as.matrix(layer_jaccard[, -"Layer", with = FALSE])
rownames(jac_mat) <- layer_jaccard$Layer
jac_mat <- jac_mat[rownames(jac_mat) != "Novel_Tx", colnames(jac_mat) != "Novel_Tx"]

# Cluster by similarity
hc <- hclust(as.dist(1 - jac_mat))
clust_order <- hc$order
jac_clust <- jac_mat[clust_order, clust_order]

# Convert to long
jac_long <- as.data.table(as.table(jac_clust))
setnames(jac_long, c("Layer", "Layer2", "Jaccard"))
jac_long[, Layer := factor(Layer, levels = rev(rownames(jac_clust)))]
jac_long[, Layer2 := factor(Layer2, levels = colnames(jac_clust))]
jac_long[, is_diag := Layer == Layer2]
jac_long[Jaccard == 0, Jaccard := NA]

p_b2 <- ggplot(jac_long, aes(x = Layer2, y = Layer)) +
  geom_tile(aes(fill = Jaccard), color = "white", linewidth = 0.7) +
  geom_tile(data = jac_long[is_diag == TRUE],
            fill = "grey88", color = "white", linewidth = 0.7) +
  # Show all non-zero values (larger font)
  geom_text(data = jac_long[!is.na(Jaccard) & is_diag == FALSE],
            aes(label = sprintf("%.3f", Jaccard)),
            size = 3.0, family = "Arial", color = "grey15",
            fontface = ifelse(jac_long[!is.na(Jaccard) & is_diag == FALSE]$Jaccard >= 0.02,
                              "bold", "plain")) +
  geom_text(data = jac_long[is_diag == TRUE],
            aes(label = "─"), size = 3.5, family = "Arial", color = "grey60") +
  scale_fill_gradientn(
    colors = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    values = c(0, 0.05, 0.2, 0.5, 1),
    name = "Jaccard\nIndex",
    na.value = "white",
    limits = c(0, 0.04),
    guide = guide_colorbar(barwidth = unit(3, "mm"), barheight = unit(12, "mm"),
                           title.position = "top", title.hjust = 0.5,
                           label.theme = element_text(size = 5.5))
  ) +
  labs(x = NULL, y = NULL,
       title = "Layer Overlap (Jaccard Index)",
       subtitle = "6 layers clustered by similarity. Novel_Tx excluded (structural zero). m6A↓ = central hub (5 overlaps).") +
  theme_nature +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
    axis.text.y = element_text(size = 6.5),
    legend.key.height = unit(10, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.position = "right"
  )

# ── Combine & Save ──
p_b <- wrap_plots(p_b1, p_b2, ncol = 2, widths = c(0.42, 0.58))

ggsave(file.path(OUTDIR, "Fig7B_LayerOverlap.pdf"), p_b,
       width = 8.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig7B_LayerOverlap.png"), p_b,
       width = 8.5, height = 3.5, dpi = 600)
cat("  Fig7B v3 Done\n")
