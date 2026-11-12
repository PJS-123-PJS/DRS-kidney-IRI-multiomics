# =============================================================================
# GSE298899 Stereo-seq projection figure (supplement candidate / Fig 3 support)
# 2026-09-28 — input: GSE298899_stereoseq_projection_{sample_means,stats}.csv
# A: deltaZ heatmap (7 direction sets x key compartments, IRI14d - Sham)
# B: per-sample Z in PT compartments (n=3/group dots; animal-level spirit)
# Canonical palette: Sham #4A4A4A / IRI #D6404C. Arial 7pt, cairo_pdf.
# Output: output/FigSX_StereoSeq_Orthogonal_Projection.{pdf,png}
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

dd <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
od <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"

sm <- fread(file.path(dd, "GSE298899_stereoseq_projection_sample_means.csv"))
st <- fread(file.path(dd, "GSE298899_stereoseq_projection_stats.csv"))

set_levels <- c("m6A_UP","m6A_DOWN_resid","m5C_DOWN","PSU_UP","Inosine_UP",
                "Nm_UP","Um_UP","Nm_DOWN")
set_labels <- c(m6A_UP="m6A↑", m6A_DOWN_resid="m6A↓ residual", m5C_DOWN="m5C↓",
                PSU_UP="PSU↑", Inosine_UP="inosine↑", Nm_UP="Nm↑",
                Um_UP="Um↑", Nm_DOWN="Nm↓")
comp_levels <- c("PT-S1/S2","PT-S3","FR-PT","Injured tubule","TAL","DCT","CNT",
                 "collecting duct","glomeruli","interstitial cells")

# ---- A: deltaZ heatmap -------------------------------------------------------
ha <- st[figure_idents %in% comp_levels]
ha[, compartment := factor(figure_idents, levels = comp_levels)]
ha[, set := factor(set, levels = set_levels, labels = set_labels[set_levels])]
ha[, lab := sprintf("%+.2f%s", deltaZ,
                    ifelse(welch_p < 0.05, "*", ""))]

pA <- ggplot(ha, aes(compartment, set, fill = deltaZ)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 1.9, family = "Arial") +
  scale_fill_gradient2(low = "#2E7D6F", mid = "white", high = "#D6404C",
                       midpoint = 0, name = "ΔZ (IRI−Sham)") +
  labs(x = NULL, y = NULL,
       title = "Stereo-seq orthogonal check (GSE298899, 14 d, n = 3 per group)") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(size = 7.5, face = "bold"),
        legend.key.height = unit(0.28, "cm"))

# ---- B: per-sample dots in PT compartments -----------------------------------
pb_d <- sm[figure_idents %in% c("PT-S1/S2","PT-S3","FR-PT") &
           set %in% c("m6A_UP","Um_UP","Nm_UP","Inosine_UP","m5C_DOWN","m6A_DOWN_resid")]
pb_d[, set := factor(set, levels = c("m6A_UP","Um_UP","Nm_UP","Inosine_UP",
                                     "m5C_DOWN","m6A_DOWN_resid"),
                     labels = set_labels[c("m6A_UP","Um_UP","Nm_UP","Inosine_UP",
                                           "m5C_DOWN","m6A_DOWN_resid")])]
pb_d[, figure_idents := factor(figure_idents, levels = c("PT-S1/S2","PT-S3","FR-PT"))]
pb_d[, group := factor(group, levels = c("Sham","IRI_14d"))]
pal <- c(Sham = "#4A4A4A", IRI_14d = "#D6404C")

pB <- ggplot(pb_d, aes(group, z, color = group)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey75", linewidth = 0.3) +
  geom_point(size = 1.5, position = position_jitter(width = 0.12, seed = 7)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "black",
               linewidth = 0.3, fatten = 1.2) +
  facet_grid(figure_idents ~ set) +
  scale_color_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Direction-set score (Z, pooled across bins)") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 6.3),
        axis.title = element_text(size = 7))

library(patchwork)
p <- pA / pB + plot_layout(heights = c(1, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))
ggsave(file.path(od, "FigSX_StereoSeq_Orthogonal_Projection.pdf"), p,
       width = 7.0, height = 7.6, device = cairo_pdf)
ggsave(file.path(od, "FigSX_StereoSeq_Orthogonal_Projection.png"), p,
       width = 7.0, height = 7.6, dpi = 600)
cat("written: FigSX_StereoSeq_Orthogonal_Projection.{pdf,png}\n")
