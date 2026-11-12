# Fig3 Panel B — 3'UTR Preference Timecourse (Nature-figure standard)
# ============================================================================
# v2 (2026-07-29) — post-audit rebuild (checkpoint v8):
#   * reads fig3b_5mod_dedup_rows.rds — NC/Sham = UNION across NC_vs_Day1/3/7
#     (v1 pooled site×comparison rows, inflating Sham counts ×1.14–1.61);
#     Day groups identical to frac_5mod_cache_v2.rds (fidelity-gated).
#   * reliability tiers now PER-TIMEPOINT n (was per-mod total, which graded
#     all 5 mods "robust" and made the fade/line-weight encoding a no-op).
#   * removed false "INOSINE: 0 sites in 3'UTR" subtitle (v1-cache leftover;
#     v2 canonical data has 106 INOSINE 3'UTR sites, 57.0%→85.7% 3'UTR).
#   * n= labels only on low-n points (n<20) — full labeling was unreadable.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales) })

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Nature theme ──
theme_nature <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1, "pt"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.key.size = unit(8, "pt"),
    strip.text = element_text(size = 6.5, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 6, color = "grey30"),
    panel.grid = element_blank(),
    panel.spacing = unit(1.5, "mm")
  )

MOD_COLORS <- c("M6A"="#D6404C","M5C"="#E8753A","PSU"="#3D7EA6","INOSINE"="#7B5294","Nm"="#2E7D6F")
REGION_FILL <- c("5'UTR"="#F2E8D0","CDS"="#B8D8B8","3'UTR"="#6090C0")
MOD_ORDER <- c("M6A","M5C","PSU","INOSINE","Nm")
TIME_ORDER <- c("Sham","Day1","Day3","Day7")
TIER_LV <- c(">=20 sites (robust)", "10-19 sites (marginal)", "<10 sites (unreliable)")

tier_of <- function(n) {
  fifelse(n >= 20, TIER_LV[1], fifelse(n >= 10, TIER_LV[2], TIER_LV[3]))
}

# ── Load deduplicated rows (cache v3 = same gated content, project-wide name) ──
frac <- readRDS(file.path(OUT_DIR, "frac_5mod_cache_v3.rds")); setDT(frac)
frac[, Group := factor(Group, levels = c("NC","Day1","Day3","Day7"), labels = TIME_ORDER)]
frac[, ModLabel := factor(Modification, levels = MOD_ORDER)]
frac[, Group := factor(Group, levels = TIME_ORDER)]

reg_counts <- frac[, .N, by = .(ModLabel, Group, Gene_Struc)]
reg_total <- frac[, .(Total = .N), by = .(ModLabel, Group)]
reg_prop <- merge(reg_counts, reg_total, by = c("ModLabel","Group"))
reg_prop[, Proportion := N / Total]

# per-timepoint reliability (B2: by 3'UTR n; B3: by bar total n)
utr3_data <- reg_prop[Gene_Struc == "three_prime_UTR"]
utr3_data[, Reliability := factor(tier_of(N), levels = TIER_LV)]
reg_total[, Reliability := factor(tier_of(Total), levels = TIER_LV)]

# ── Panel B₁: unique site count bars ──
# labels: M6A/PSU/Nm above bars (never horizontally adjacent); M5C/INOSINE
# white inside bar top — avoids all adjacent-label collisions on log scale
p_count <- ggplot(reg_total, aes(x = Group, y = Total, fill = ModLabel)) +
  geom_col(position = position_dodge(0.8), width = 0.65) +
  geom_text(data = reg_total[ModLabel %in% c("M6A","PSU","Nm")],
            aes(label = Total, color = ModLabel),
            position = position_dodge(0.8), vjust = -0.3, size = 2, family = "Arial") +
  geom_text(data = reg_total[ModLabel %in% c("M5C","INOSINE")],
            aes(label = Total),
            position = position_dodge(0.8), vjust = 1.15, size = 1.8,
            family = "Arial", color = "white") +
  scale_fill_manual(values = MOD_COLORS, name = NULL) +
  scale_color_manual(values = MOD_COLORS, guide = "none") +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0.05, 0.15))) +
  annotation_logticks(sides = "l", alpha = 0.25, linewidth = 0.2) +
  labs(x = NULL, y = "Unique sites (log10)",
       title = "Differential Site Count",
       subtitle = "unique sites; filtered: fraction>0.7, |ES|>=0.1, p<0.05") +
  theme_nature + theme(legend.position = "right", legend.margin = margin(l = -5))

# ── Panel B₂: 3'UTR trend (per-timepoint reliability) ──
# n= labels on low-n points only; label goes below the point if another mod's
# point sits just above (else it would visually attach to the wrong series)
lab_dt <- utr3_data[N < 20]
lab_above <- lab_dt[, {
  y0 <- Proportion * 100
  crowded <- utr3_data[ModLabel != .BY$ModLabel & Group == .BY$Group &
                       Proportion * 100 > y0 + 1 & Proportion * 100 < y0 + 9]
  .(above = nrow(crowded) == 0)
}, by = .(ModLabel, Group)]
lab_dt <- merge(lab_dt, lab_above, by = c("ModLabel","Group"))

p_utr3 <- ggplot(utr3_data, aes(x = Group, y = Proportion * 100,
                                 color = ModLabel, group = ModLabel)) +
  geom_line(aes(alpha = Reliability, linewidth = Reliability)) +
  geom_point(aes(size = N, alpha = Reliability)) +
  geom_text(data = lab_dt[above == TRUE],
            aes(label = paste0("n=", N)),
            vjust = -1.3, size = 2, family = "Arial", show.legend = FALSE) +
  geom_text(data = lab_dt[above == FALSE],
            aes(label = paste0("n=", N)),
            vjust = 1.8, size = 2, family = "Arial", show.legend = FALSE) +
  scale_color_manual(values = MOD_COLORS, name = NULL) +
  scale_alpha_manual(values = c(">=20 sites (robust)"=0.9, "10-19 sites (marginal)"=0.45, "<10 sites (unreliable)"=0.18), guide = "none") +
  scale_linewidth_manual(values = c(">=20 sites (robust)"=0.5, "10-19 sites (marginal)"=0.3, "<10 sites (unreliable)"=0.15), guide = "none") +
  scale_size_continuous(range = c(0.8, 3.5), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "3'UTR sites (%)",
       title = "3'UTR Preference", subtitle = "faded = low-n timepoint (3'UTR n<20)") +
  theme_nature + theme(legend.position = "right", legend.margin = margin(l = -5))

# ── Panel B₃: Regional stacked bars (per-bar reliability) ──
reg_prop[, Region := factor(Gene_Struc,
  levels = c("three_prime_UTR","CDS","five_prime_UTR"),
  labels = c("3'UTR","CDS","5'UTR"))]
reg_prop <- merge(reg_prop, reg_total[, .(ModLabel, Group, Reliability)],
                  by = c("ModLabel","Group"), all.x = TRUE)

p_stack <- ggplot(reg_prop, aes(x = Group, y = Proportion, fill = Region, alpha = Reliability)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.15) +
  facet_wrap(~ ModLabel, ncol = 5) +
  scale_fill_manual(values = REGION_FILL, name = NULL) +
  scale_alpha_manual(values = c(">=20 sites (robust)"=1, "10-19 sites (marginal)"=0.6, "<10 sites (unreliable)"=0.3), guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Proportion",
       title = "Regional Distribution (5'UTR / CDS / 3'UTR)",
       subtitle = "faded = low-n timepoint (total n<20)") +
  theme_nature + theme(legend.position = "bottom", legend.margin = margin(t = -2))

# ── Assemble ──
top_row <- wrap_plots(p_count, p_utr3, ncol = 2, widths = c(0.48, 0.52))
combined <- wrap_plots(top_row, p_stack, ncol = 1, heights = c(0.38, 0.62)) +
  plot_annotation(
    title = "Fig3B: Site Count Dynamics & 3'UTR Preference Across IRI Time Course",
    subtitle = "Upper: unique site count (log) + 3'UTR proportion trend. Lower: regional distribution stacked bars.",
    caption = paste0("n=", format(nrow(frac), big.mark=","), " unique differential sites ",
                     "(Sham = union across NC-vs-Day1/3/7; no double-counting). Nm=Am+Cm+Gm+Um.\n",
                     "Reliability tiers by per-timepoint n: robust>=20, marginal=10-19, unreliable<10 (faded)."),
    theme = theme(plot.title = element_text(size = 8.5, face = "bold", family = "Arial"),
                  plot.subtitle = element_text(size = 6.5, color = "grey30", family = "Arial"),
                  plot.caption = element_text(size = 5.5, color = "grey50", hjust = 0, family = "Arial"))
  )

W_MM <- 183; H_MM <- 120
ggsave(file.path(OUT_DIR, "Fig3B_3UTR_Timecourse.pdf"), combined, width = W_MM/25.4, height = H_MM/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3B_3UTR_Timecourse.png"), combined, width = W_MM/25.4, height = H_MM/25.4, dpi = 600)
cat(sprintf("Fig3B v2 (dedup) saved: %.0f × %.0f mm\n", W_MM, H_MM))
