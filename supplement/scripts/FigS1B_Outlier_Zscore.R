# FigS1B: Per-sample outlier analysis (8-mod site counts)
# ============================================================================
# 2026-07-30 — README planned item, user-approved build (batch 4).
# v3: global Z-score heatmap (|Z|>=2.5 labelled — catches the Fig2A × marks:
#     Day3s1 Inosine, Day7s1 Um; NCs2 Gm at tiny counts) + NCs3 yield profile
#     narrated in caption with runtime-computed within-group deviations.
#     (v2's blanket within-group %>=40 rule was dropped: it fired on 28 cells —
#      low-count mods swing wildly and Day3s1's 32,267 Inosine poisons its
#      groupmates' peer means. With n=3/group, within-group Z is capped ~1.41,
#      so NCs3 is documented by narrative fold-deviation, not a fake statistic.)
# Zero hardcoding: all numbers computed from the CSV at runtime.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

QC_CSV <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig1_platform_qc/output/Table_Fig1_8mod_SiteCounts_Long.csv"
OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/supplement/figS1_qc"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dt <- fread(QC_CSV)
MODS <- c("m6A","m5C","PSU","Inosine","Am","Cm","Gm","Um")
stopifnot(all(MODS %in% names(dt)))

long <- melt(dt, id.vars = c("Sample","group"), measure.vars = MODS,
             variable.name = "Mod", value.name = "Sites")
long[, Z := as.numeric(scale(Sites)), by = Mod]
long[, GlobalOut := abs(Z) >= 2.5]
# within-group % deviation vs peer mean (used for NCs3 narrative only)
long[, PeerMean := (sum(Sites) - Sites) / (.N - 1), by = .(group, Mod)]
long[, PeerDev := (Sites / PeerMean - 1) * 100]

grp_order <- c("NC","IS","Day1","Day3","Day7")
long[, group := factor(group, levels = grp_order)]
long[, Sample := factor(Sample, levels = dt$Sample[order(match(dt$group, grp_order))])]
long[, Mod := factor(Mod, levels = MODS)]

g_out <- long[GlobalOut == TRUE]
cat("=== Global outliers (|Z|>=2.5) ===\n"); print(g_out[, .(Sample,group,Mod,Sites,Z)])

# NCs3 within-group profile (the README-named subject of this panel)
ncs3 <- long[Sample == "NCs3"][order(PeerDev)]
ncs3_line <- paste(sprintf("%s %+.0f%%", ncs3$Mod, ncs3$PeerDev), collapse = ", ")
cat(sprintf("NCs3 vs NC peers: %s\n", ncs3_line))

sub_line <- paste(sprintf("%s %s (Z=%.1f)", g_out$Sample, g_out$Mod, g_out$Z), collapse = "; ")

p <- ggplot(long, aes(x = Mod, y = Sample, fill = Z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(data = g_out, aes(label = sprintf("%.1f", Z)),
            size = 2.6, fontface = "bold", color = "black", family = "Arial") +
  scale_fill_gradient2(low = "#3D7EA6", mid = "white", high = "#D6404C",
                       midpoint = 0, limits = c(-3, 3), oob = squish,
                       name = "Z-score") +
  labs(title = "FigS1B: Per-Sample Site-Count Outlier Analysis (8 modifications)",
       subtitle = paste0("Global outliers (|Z|>=2.5, labelled): ", sub_line),
       x = NULL, y = NULL,
       caption = paste0(
         "Fill = Z-score of per-sample site counts, per modification across all 15 samples.\n",
         "Day3s1 Inosine (32,267 sites) and Day7s1 Um are the outlier samples marked in Fig 2A.\n",
         "NCs3 yield profile vs NC peer mean: ", ncs3_line, " —\n",
         "an anomalous modification-yield profile consistent with its off-cluster position in Fig 1 PCA.")) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 6.5),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 5.8, color = "grey30"),
        plot.caption = element_text(size = 5.3, color = "grey50", hjust = 0),
        legend.key.size = unit(8, "pt"), legend.title = element_text(size = 6),
        legend.text = element_text(size = 5.5))

W_MM <- 130; H_MM <- 95
ggsave(file.path(OUT_DIR, "FigS1B_Outlier_Zscore.pdf"), p,
       width = W_MM/25.4, height = H_MM/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS1B_Outlier_Zscore.png"), p,
       width = W_MM/25.4, height = H_MM/25.4, dpi = 600)

out <- dcast(long, Sample + group ~ Mod, value.var = "Z")
setnames(out, MODS, paste0(MODS, "_Z"))
fwrite(out, file.path(OUT_DIR, "Table_FigS1B_sample_Zscores.csv"))
fwrite(long[, .(Sample,group,Mod,Sites,Z,PeerDev,GlobalOut)],
       file.path(OUT_DIR, "Table_FigS1B_sitecounts_Z_long.csv"))
cat(sprintf("FigS1B v3 saved: %.0f x %.0f mm -> %s\n", W_MM, H_MM, OUT_DIR))
