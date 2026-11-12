#!/usr/bin/env Rscript
# Fig7E v3 — Local: Mixed Boxplot + Line using server-computed boxplot stats
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig7_multiomics/output"

bp_stats <- fread(text = "Time,Layer,N,ymin,lower,middle,upper,ymax,mean_val
IS,DRS Modification,11352,-0.0707,0.3259,0.4574,0.5902,0.9868,0.3826
Day1,DRS Modification,33649,-2.4470,-0.4790,0.5046,0.8330,2.8009,0.2673
Day3,DRS Modification,23931,-1.8195,-0.2970,0.4824,0.7180,2.2405,0.3159
Day7,DRS Modification,20925,-1.7608,-0.2821,0.4770,0.7037,2.1824,0.3127")
bp_stats[, Time := factor(Time, levels = c("IS", "Day1", "Day3", "Day7"))]

line_dt <- data.table(
  Time = factor(rep(c("IS", "Day1", "Day3", "Day7"), 3),
                levels = c("IS", "Day1", "Day3", "Day7")),
  Layer = rep(c("DRS Expression", "ATAC Chromatin", "Metabolome"), each = 4),
  Pct = c(3.2, 100.0, 93.7, 64.4,
          27.2, 40.6, 100.0, 54.0,
          5.9, 100.0, 52.9, 23.5),
  N_Changed = c(31, 958, 898, 617, 142, 212, 522, 282, 5, 85, 45, 20)
)
line_dt[, N_label := ifelse(N_Changed >= 50, paste0("n=", comma(N_Changed)), "")]

LAYER_COLORS <- c(
  "DRS Expression"   = "#D6404C",
  "DRS Modification" = "#E8753A",
  "ATAC Chromatin"   = "#3D7EA6",
  "Metabolome"       = "#2E7D6F"
)

pct_to_log2fc <- function(pct) (pct / 100) * 5 - 2
log2fc_to_pct <- function(lfc) (lfc + 2) / 5 * 100

theme_nature <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1, "pt"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 6, color = "grey30"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 7, color = "black"),
    axis.title.y.left = element_text(size = 7, color = "#E8753A"),
    axis.title.y.right = element_text(size = 7, color = "grey40"),
    plot.margin = margin(4, 4, 4, 4)
  )

p <- ggplot() +

  # Background shading
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -3, ymax = 4.5,
           fill = "#D6404C", alpha = 0.08) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -3, ymax = 4.5,
           fill = "#3D7EA6", alpha = 0.08) +

  # BOXPLOT: DRS Modification (manual from summary stats)
  # Box body (lower to upper quartile)
  geom_rect(
    data = bp_stats,
    aes(xmin = as.numeric(Time) - 0.2, xmax = as.numeric(Time) + 0.2,
        ymin = lower, ymax = upper),
    fill = alpha("#E8753A", 0.3), color = "#E8753A", linewidth = 0.4
  ) +
  # Median line
  geom_segment(
    data = bp_stats,
    aes(x = as.numeric(Time) - 0.2, xend = as.numeric(Time) + 0.2,
        y = middle, yend = middle),
    color = "#E8753A", linewidth = 0.8
  ) +
  # Lower whisker
  geom_segment(
    data = bp_stats,
    aes(x = as.numeric(Time), xend = as.numeric(Time),
        y = ymin, yend = lower),
    color = "#E8753A", linewidth = 0.4
  ) +
  # Upper whisker
  geom_segment(
    data = bp_stats,
    aes(x = as.numeric(Time), xend = as.numeric(Time),
        y = upper, yend = ymax),
    color = "#E8753A", linewidth = 0.4
  ) +

  # Median line + point
  geom_line(data = bp_stats, aes(x = Time, y = middle, group = 1,
            color = "DRS Modification"), linewidth = 1.2, alpha = 0.9) +
  geom_point(data = bp_stats, aes(x = Time, y = middle,
             color = "DRS Modification"), size = 3.5, shape = 21,
             fill = "white", stroke = 1.5) +

  # n= labels for boxplot
  geom_text(data = bp_stats, aes(x = Time, y = 3.6,
            label = paste0("n=", comma(N), "\ntranscripts")),
            size = 2.5, family = "Arial", color = "#E8753A",
            vjust = 0, lineheight = 0.85) +

  # LINES: Expression, ATAC, Metabolome
  geom_line(data = line_dt,
            aes(x = Time, y = pct_to_log2fc(Pct), color = Layer, group = Layer),
            linewidth = 1.0, alpha = 0.9, linetype = "dashed") +
  geom_point(data = line_dt,
             aes(x = Time, y = pct_to_log2fc(Pct), color = Layer,
                 size = N_Changed), alpha = 0.8) +

  # N labels for lines
  ggrepel::geom_text_repel(
    data = line_dt[N_label != ""],
    aes(x = Time, y = pct_to_log2fc(Pct), label = N_label, color = Layer),
    size = 2.5, family = "Arial", fontface = "bold",
    box.padding = 0.5, force = 2, max.overlaps = 20, seed = 42,
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = LAYER_COLORS, name = NULL,
    breaks = c("DRS Modification", "DRS Expression", "ATAC Chromatin", "Metabolome"),
    guide = guide_legend(override.aes = list(
      linetype = c("solid", "dashed", "dashed", "dashed"),
      linewidth = c(1.2, 1.0, 1.0, 1.0),
      shape = c(21, 16, 16, 16),
      fill = c("#E8753A", NA, NA, NA)
    ))
  ) +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  scale_y_continuous(
    name = "DRS Modification log2FC (boxplot)",
    limits = c(-3, 4.5),
    sec.axis = sec_axis(
      transform = ~ log2fc_to_pct(.),
      name = "Relative Change (% of Layer Maximum)",
      breaks = seq(0, 100, 25)
    )
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +

  annotate("text", x = 2, y = 4.3,
           label = "Acute Functional\nDRS peaks Day1",
           color = "#D6404C", size = 2.3, family = "Arial",
           fontface = "bold", hjust = 0.5) +
  annotate("text", x = 3, y = 4.3,
           label = "Delayed Structural\nATAC peaks Day3",
           color = "#3D7EA6", size = 2.3, family = "Arial",
           fontface = "bold", hjust = 0.5) +

  labs(x = NULL,
       title = "Omics Layers Decouple in Time: Acute Functional vs Delayed Structural",
       subtitle = "Boxplot: DRS Modification per-transcript log2FC (89K transcripts, whiskers = 1.5xIQR). Dashed lines: % of layer maximum (each layer normalized to its own peak = 100%). Expression = gene-level, Modification = site-level.") +
  theme_nature +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))

ggsave(file.path(OUT_DIR, "Fig7E_Boxplot_Line.pdf"), p,
       width = 8.5, height = 3.8, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig7E_Boxplot_Line.png"), p,
       width = 8.5, height = 3.8, dpi = 600)
cat("Fig7E v3 boxplot+line saved\n")
