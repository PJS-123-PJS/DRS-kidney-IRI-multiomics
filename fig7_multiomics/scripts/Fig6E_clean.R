#!/usr/bin/env Rscript
# Fig7E v4 — Clean: Boxplot (Modification) + subtle line trend
# Removed: dual axis, ggrepel labels, long subtitle, secondary annotations
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig7_multiomics/output"

# Boxplot stats (from server panelE_trans_frac, 89K transcripts)
bp_stats <- fread(text = "Time,N,ymin,lower,middle,upper,ymax
IS,11352,-0.07,0.33,0.46,0.59,0.99
Day1,33649,-2.45,-0.48,0.50,0.83,2.80
Day3,23931,-1.82,-0.30,0.48,0.72,2.24
Day7,20925,-1.76,-0.28,0.48,0.70,2.18")
bp_stats[, Time := factor(Time, levels = c("IS", "Day1", "Day3", "Day7"))]

# Line data (3 other layers, % of max → mapped to log2FC range for overlay)
line_dt <- data.table(
  Time = factor(rep(c("IS", "Day1", "Day3", "Day7"), 3),
                levels = c("IS", "Day1", "Day3", "Day7")),
  Layer = rep(c("DRS Expression", "ATAC Chromatin", "Metabolome"), each = 4),
  Pct = c(3.2, 100.0, 93.7, 64.4, 27.2, 40.6, 100.0, 54.0, 5.9, 100.0, 52.9, 23.5)
)
# Map Pct (0-100) → log2FC range (-2 to 3) for overlay on same axis
line_dt[, y_scaled := Pct / 100 * 5 - 2]

LAYER_COLORS <- c(
  "DRS Modification" = "#E8753A",
  "DRS Expression"   = "#D6404C",
  "ATAC Chromatin"   = "#3D7EA6",
  "Metabolome"       = "#2E7D6F"
)

theme_clean <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1.2, "pt"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, color = "grey40"),
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(6, 6, 4, 6)
  )

p <- ggplot() +

  # Subtle background: Acute (Day1) + Delayed (Day3)
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -3, ymax = 3.5,
           fill = "#D6404C", alpha = 0.06) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -3, ymax = 3.5,
           fill = "#3D7EA6", alpha = 0.06) +

  # ── BOXPLOT: DRS Modification ──
  geom_rect(data = bp_stats,
            aes(xmin = as.numeric(Time) - 0.22, xmax = as.numeric(Time) + 0.22,
                ymin = lower, ymax = upper),
            fill = alpha("#E8753A", 0.35), color = "#E8753A", linewidth = 0.5) +
  geom_segment(data = bp_stats,
               aes(x = as.numeric(Time) - 0.22, xend = as.numeric(Time) + 0.22,
                   y = middle, yend = middle),
               color = "#E8753A", linewidth = 0.9) +
  # Whiskers
  geom_segment(data = bp_stats,
               aes(x = as.numeric(Time), xend = as.numeric(Time),
                   y = ymin, yend = lower), color = "#E8753A", linewidth = 0.4) +
  geom_segment(data = bp_stats,
               aes(x = as.numeric(Time), xend = as.numeric(Time),
                   y = upper, yend = ymax), color = "#E8753A", linewidth = 0.4) +

  # n= count above each box
  geom_text(data = bp_stats, aes(x = Time, y = 3.0,
            label = paste0("n=", comma(N))),
            size = 2.8, family = "Arial", color = "#E8753A", vjust = 0) +

  # ── SUBTLE TREND LINES (no points, no labels) ──
  geom_line(data = line_dt,
            aes(x = Time, y = y_scaled, color = Layer, group = Layer),
            linewidth = 0.6, alpha = 0.6, linetype = "dotted") +

  # ── Scales ──
  scale_color_manual(values = LAYER_COLORS,
                     breaks = names(LAYER_COLORS)) +
  scale_y_continuous(
    name = "DRS Modification log2FC",
    limits = c(-3, 3.5),
    breaks = seq(-3, 3, 1)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +

  # ── Minimal annotations ──
  annotate("text", x = 2.0, y = 3.3,
           label = "Acute", color = "#D6404C", size = 2.8,
           family = "Arial", fontface = "bold", hjust = 0.5) +
  annotate("text", x = 3.0, y = 3.3,
           label = "Delayed", color = "#3D7EA6", size = 2.8,
           family = "Arial", fontface = "bold", hjust = 0.5) +

  labs(x = NULL, y = "DRS Modification log2FC",
       title = "Temporal Decoupling of Omics Layers",
       subtitle = "Boxplot: DRS Modification (8 mods, 89K transcripts). Dotted: trend for Expression, ATAC, Metabolome (% max scaled).") +

  theme_clean +
  theme(legend.position = "bottom")

# ── Save ──
ggsave(file.path(OUT_DIR, "Fig7E_TimeDecoupling_Core.pdf"), p,
       width = 7.5, height = 3.2, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig7E_TimeDecoupling_Core.png"), p,
       width = 7.5, height = 3.2, dpi = 600)
cat("Fig7E v4 clean saved\n")
