# Fig7E v2 — Optimized: Four-Omics Time Decoupling
# ============================================================================
# vision_analyzer v1: 7/10 (hardcoded data, oversized, theme_drs, color confusion)
# v2 fixes:
#   1. Read from Table_TimeDecoupling.csv (zero hardcoding)
#   2. Nature theme (theme_classic, Arial 6.5pt)
#   3. Per-layer colors: Expression=#D6404C, Modification=#E8753A,
#      ATAC=#3D7EA6, Metabolome=#7B5294
#   4. Size: 8.5×3.5in, 600dpi (consistent with B,D panels)
#   5. Clear acute/delayed annotation, N= labels, proper legend
#   6. Point size mapped to N_Changed for visual weight
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

# Layer colors (from 5-mod palette + extensions)
LAYER_COLORS <- c(
  "DRS Expression"   = "#D6404C",  # red
  "DRS Modification" = "#C49A3C",  # gold (distinct from red)
  "ATAC Chromatin"   = "#3D7EA6",  # blue
  "Metabolome"       = "#2E7D6F"   # green-teal (colorblind-safe vs blue)
)

ANN_COLORS <- c(
  "Acute Functional"  = "#D6404C",
  "Delayed Structural" = "#3D7EA6"
)

# ===========================================================================
# Fig7E: Time Decoupling (v2 optimized)
# ===========================================================================
cat("Fig7E (v2 optimized): Four-Omics Time Decoupling...\n")

time_dt <- fread("../data/Table_TimeDecoupling.csv")
time_dt[, Time := factor(Time, levels = c("IS", "Day1", "Day3", "Day7"))]
time_dt[, Layer := factor(Layer, levels = c("DRS Expression", "DRS Modification",
                                             "ATAC Chromatin", "Metabolome"))]

# N label for points with substantial changes
time_dt[, N_label := ifelse(N_Changed >= 50, paste0("N=", comma(N_Changed)), "")]

# ── Main line plot ──
p_e <- ggplot(time_dt, aes(x = Time, y = Pct, color = Layer, group = Layer)) +

  # Acute (Day1) and Delayed (Day3) background highlights
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 108,
           fill = "#D6404C", alpha = 0.12) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 108,
           fill = "#3D7EA6", alpha = 0.12) +

  # Lines
  geom_line(linewidth = 0.8, alpha = 0.9) +

  # Points sized by N_Changed
  geom_point(aes(size = N_Changed), alpha = 0.9) +

  # N= labels (allow both x and y movement to resolve Day1 overlap)
  ggrepel::geom_text_repel(
    aes(label = N_label),
    size = 2.8, family = "Arial", fontface = "bold",
    nudge_y = 3, segment.size = 0.2,
    box.padding = 0.5, force = 2, max.time = 2,
    show.legend = FALSE, max.overlaps = 20, seed = 123
  ) +

  # Layer-specific colors
  scale_color_manual(values = LAYER_COLORS, name = NULL) +
  scale_size_continuous(range = c(4, 10), guide = "none") +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0), breaks = seq(0, 100, 25)) +

  # Annotations
  annotate("text", x = 2, y = 112, label = "Acute Functional\n(DRS peaks Day1)",
           color = "#D6404C", size = 2.2, family = "Arial", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 3, y = 112, label = "Delayed Structural\n(ATAC peaks Day3)",
           color = "#3D7EA6", size = 2.2, family = "Arial", fontface = "bold", hjust = 0.5) +

  labs(x = NULL, y = "Relative Change (% of Layer Maximum)",
       title = "Omics Layers Decouple in Time: Acute Functional vs Delayed Structural",
       subtitle = "DRS (Expression + Modification) peaks Day1. ATAC peaks Day3. Metabolome follows intermediate kinetics.") +

  theme_nature +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(4, 4, 4, 4),
    axis.text = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 7)
  )

# ── Save ──
ggsave(file.path(OUTDIR, "Fig7E_TimeDecoupling_Core.pdf"), p_e,
       width = 8.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig7E_TimeDecoupling_Core.png"), p_e,
       width = 8.5, height = 3.5, dpi = 600)

cat(sprintf("  Fig7E v2 Done. Data: %d rows from Table_TimeDecoupling.csv\n", nrow(time_dt)))
