#!/usr/bin/env Rscript
# Fig7A v9 — PT-S3 Multi-Omics: Clean faceted bars
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

OUT_DIR <- "../output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Fig7A v9 ===\n")

# ── Data ──
layer_sizes <- fread("../data/Table_Layer_Sizes.csv")
multi_genes <- fread("../data/Table_MultiLayer_Genes_validated.csv")
metab_dyn <- fread("../output/Fig7G_Metabolite_Dynamic.csv")

pt_markers <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1",
                "Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc13a1","Slc22a7","Slc13a3","Lrp2","Cubn","Miox")
n_pt <- multi_genes[Gene %in% pt_markers, .N]

# Build per-ring
rings <- list()

# A: Expression
dA <- multi_genes[, .N, by = N_Layers][order(N_Layers)]
rings[["A  Expression"]] <- data.table(
  item = paste0(dA$N_Layers, "-layer genes"), N = dA$N,
  color = c("#F4A582", "#D6604D", "#B2182B")[seq_len(nrow(dA))],
  total = sum(dA$N)
)

# B: Modifications
dB <- layer_sizes[Category == "Modification"][1:4]
rings[["B  Modifications"]] <- data.table(
  item = gsub("[↓↑]", "", dB$Layer), N = dB$N_Genes,
  color = c("#D6404C", "#E8753A", "#3D7EA6", "#C49A3C")[1:4],
  total = sum(dB$N_Genes)
)

# C: AS + PolyA
dC <- layer_sizes[Category %in% c("AS", "PolyA")]
rings[["C  AS + PolyA"]] <- data.table(
  item = gsub("_diff", "", dC$Layer), N = dC$N_Genes,
  color = c("#5B8FA8", "#2E7D6F")[seq_len(nrow(dC))],
  total = sum(dC$N_Genes)
)

# D: Metabolome
dV <- metab_dyn[!is.na(log2FC), .N, by = Category][order(-N)]
met_cols <- c("Methyl Cycle"="#D6404C", "Redox"="#3D7EA6", "Energy"="#E8753A")
rings[["D  Metabolome"]] <- data.table(
  item = dV$Category, N = dV$N,
  color = met_cols[dV$Category],
  total = sum(dV$N)
)

# Combine
all_dt <- rbindlist(lapply(names(rings), function(nm) {
  rings[[nm]][, ring := nm]
}), fill = TRUE)
all_dt[, ring := factor(ring, levels = names(rings))]
all_dt[, pct := round(N / total * 100, 1), by = ring]

# Total annotations
totals <- all_dt[, .(total = total[1], maxN = max(N)), by = ring]

# ── Theme ──
theme_bar <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(size = 6.5, color = "black"),
    axis.title.x = element_text(size = 7),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 6.5, color = "grey40", hjust = 0.5, margin = margin(b = 4)),
    plot.caption = element_text(size = 6, color = "grey50", hjust = 0.5),
    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 12, 2, 4)
  )

# ── Plot ──
p <- ggplot(all_dt, aes(x = N, y = reorder(item, N))) +
  geom_col(aes(fill = I(color)), width = 0.5) +
  # Count label on bar
  geom_text(aes(label = paste0(comma(N), " (", pct, "%)"),
                x = ifelse(N < 20, N * 3, N * 1.06)),
            hjust = 0, size = 2.6, family = "Arial", color = "grey30") +
  # Total per ring
  geom_text(data = totals,
            aes(x = maxN * 8, y = -Inf, label = paste0("Total: ", comma(total))),
            hjust = 1, vjust = -1, size = 2.6, family = "Arial", color = "grey50") +
  scale_x_log10(
    labels = comma,
    breaks = c(1, 10, 100, 1000, 10000),
    limits = c(1, 25000),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  annotation_logticks(sides = "b", linewidth = 0.2, alpha = 0.3,
                      short = unit(0.4, "mm"), mid = unit(0.6, "mm"),
                      long = unit(1, "mm")) +
  facet_wrap(~ ring, ncol = 1, scales = "free_y") +
  labs(x = "Number of Features (log10 scale)", y = NULL,
       title = "PT-S3 Multi-Omics Convergence in a Single DRS Run",
       subtitle = sprintf("4 omics layers captured simultaneously. %d/%d PT-S3 markers detected across layers. ", n_pt, length(pt_markers)),
       caption = "m6A=N6-methyladenosine, m5C=5-methylcytosine, PSU=pseudouridine, Um=2'-O-methyluridine, AS=alternative splicing, PolyA=polyadenylation. All panels share identical log10 x-axis range for cross-layer comparison. Metabolome features are sparse (n=10) relative to transcriptomic layers.") +
  theme_bar

ggsave(file.path(OUT_DIR, "Fig7A_Linear_4Omics.pdf"), p,
       width = 7, height = 5.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig7A_Linear_4Omics.png"), p,
       width = 7, height = 5.5, dpi = 600)
cat("  Fig7A v9 saved (7×5.5in)\n")
