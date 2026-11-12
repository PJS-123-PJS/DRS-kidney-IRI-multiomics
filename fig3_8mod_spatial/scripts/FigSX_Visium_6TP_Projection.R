# =============================================================================
# Fig S-candidate: SIX-timepoint Visium projection (GSE269622, direct
# module-score caliber) — extension of Fig 3F time coverage.
# 2026-09-28
# A: PT-compartment trajectories of 7 validated direction sets across
#    Sham/4h/12h/2d/14d/42d; scRNA-weighted anchors of the Fig 3F pipeline
#    (Sham +0.445, 12h->Early +1.156, 14d->Late -0.060) overlaid as open
#    triangles. 2026-10-06: y-axis label neutralized to "PT-compartment
#    score (Z)" (two calibers share the axis); anchor tag renamed from
#    internal jargon "canonical anchor" -> "scRNA-weighted anchor (Fig 3F)".
# B: CT x mod heatmaps for the three NEW timepoints (4h/2d/42d) not covered
#    by main Fig 3F.
# Arial 7pt, cairo_pdf. Input: data/visium_6tp_{PT_summary,ctmod_z}.csv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

dd <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
od <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"

pt  <- fread(file.path(dd, "visium_6tp_PT_summary.csv"))
ctm <- fread(file.path(dd, "visium_6tp_ctmod_z.csv"))

MODS <- c("m6A_UP","m6A_DOWN_resid","m5C_DOWN","PSU_UP","Inosine_UP",
          "Nm_UP","Nm_DOWN")
mod_lab <- c(m6A_UP = "m6A↑", m6A_DOWN_resid = "m6A↓ residual",
             m5C_DOWN = "m5C↓", PSU_UP = "PSU↑", Inosine_UP = "inosine↑",
             Nm_UP = "Nm↑", Nm_DOWN = "Nm↓")
mod_col <- c(m6A_UP = "#D6404C", m6A_DOWN_resid = "#8A8A8A",
             m5C_DOWN = "#2E7D6F", PSU_UP = "#E8753A",
             Inosine_UP = "#799acb", Nm_UP = "#9d86bd", Nm_DOWN = "#ca95c1")
tp_levels <- c("Sham","Hour4","Hour12","Day2","Day14","Week6")
tp_lab    <- c(Sham="Sham", Hour4="4h", Hour12="12h", Day2="2d",
               Day14="14d", Week6="42d")

# ---- A: PT trajectories ------------------------------------------------------
long <- melt(pt, id.vars = "timepoint",
             measure.vars = paste0(MODS, "_PT_Z"),
             variable.name = "mod", value.name = "z")
long[, mod := sub("_PT_Z$", "", mod)]
long[, timepoint := factor(timepoint, levels = tp_levels)]

anchors <- data.table(
  timepoint = factor(c("Sham","Hour12","Day14"), levels = tp_levels),
  z = c(0.44546, 1.155633, -0.060388),
  src = "canonical scRNA-weight (m6A↑)")

pA <- ggplot(long, aes(timepoint, z, group = mod)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey80", linewidth = 0.3) +
  geom_line(aes(color = mod, linewidth = mod)) +
  geom_point(aes(color = mod), size = 1.3) +
  geom_point(data = anchors, aes(timepoint, z), inherit.aes = FALSE,
             shape = 2, size = 2.2, stroke = 0.9, color = "black", fill = NA) +
  scale_color_manual(values = mod_col, labels = mod_lab, name = NULL) +
  scale_linewidth_manual(values = c("m6A_UP" = 0.9, "m6A_DOWN_resid" = 0.9,
                                    rep(0.45, 5)), guide = "none") +
  scale_x_discrete(labels = tp_lab) +
  labs(x = NULL, y = "PT-compartment score (Z)",
       title = "Six-timepoint PT trajectory (Visium GSE269622)") +
  annotate("text", x = 0.72, y = 1.18, label = "△ scRNA-weighted anchor (Fig 3F)",
           size = 2.1, family = "Arial", hjust = 0) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        legend.position = c(0.80, 0.72),
        legend.background = element_blank(),
        legend.key.width = unit(0.32, "cm"),
        legend.key.height = unit(0.26, "cm"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 7.5, face = "bold"))

# ---- B: CT x mod heatmaps for NEW timepoints ---------------------------------
newtp <- c("Hour4","Day2","Week6")
ct_levels <- c("PT","TAL","DCT","CNT_CD","ICB","Podo","Endo","Macrophage",
               "T_NK","B_cell","Neutrophil","Proliferating")
hb <- ctm[timepoint %in% newtp]
hb[, timepoint := factor(timepoint, levels = tp_levels, labels = tp_lab)]
hb[, celltype := factor(celltype, levels = rev(ct_levels))]
hb[, mod := factor(mod, levels = MODS, labels = mod_lab[MODS])]
hb[, lab := sprintf("%+.2f", meanZ)]

pB <- ggplot(hb, aes(mod, celltype, fill = meanZ)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = lab), size = 1.55, family = "Arial") +
  scale_fill_gradient2(low = "#2E7D6F", mid = "white", high = "#D6404C",
                       midpoint = 0, name = "mean Z") +
  facet_wrap(~ timepoint, nrow = 1) +
  labs(x = NULL, y = NULL,
       title = "Cell-type x direction-set enrichment at timepoints not in Fig 3F") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7.5, face = "bold"),
        legend.key.height = unit(0.28, "cm"))

p <- pA / pB + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))

ggsave(file.path(od, "FigSX_Visium_6TP_Projection.pdf"), p,
       width = 7.0, height = 6.4, device = cairo_pdf)
ggsave(file.path(od, "FigSX_Visium_6TP_Projection.png"), p,
       width = 7.0, height = 6.4, dpi = 600)
cat("written: FigSX_Visium_6TP_Projection.{pdf,png}\n")

# ---- console summary ----------------------------------------------------------
cat("\n=== PT trajectory (direct caliber) ===\n")
print(dcast(long, timepoint ~ mod, value.var = "z"))
cat("\nanchors (canonical scRNA-weight): Sham +0.445 / Hour12->Early +1.156 / Day14->Late -0.060\n")
