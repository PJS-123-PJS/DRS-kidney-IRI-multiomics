# =============================================================================
# Fig2B v2: HC site-count trajectory + count-change vs threshold (new path)
# 2026-07-22 — replaces Fig2B_Suppressive_DoubleBar ("Both decreased" was
# FALSE for m6A). Left: HC70 trajectory (m6A stable -6.5%, m5C real -34%).
# Right: Sham->Day1 %count-change vs fraction gate from Part C bins — the
# m6A "loss" deflates with stringency, the m5C loss persists.
# Mock-testable via env FIG2B_DATA_DIR pointing at a dir with the 4 tsvs.
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
source("Fig2_config.R")

dd <- Sys.getenv("FIG2B_DATA_DIR", "../../fig3_8mod_spatial/data")
out_dir <- "../outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(Sham = "#E64B35", IS = "#F39B7F", Day1 = "#4DBBD5",
                  Day3 = "#00A087", Day7 = "#3C5488")
mod_colors <- c(m6A = "#D6404C", m5C = "#E8753A")
grp_lv <- c("Sham", "IS", "Day1", "Day3", "Day7")

read_hc <- function(f, mod) {
  d <- fread(file.path(dd, f))
  d[, mod := mod]
  d[, group := factor(group, levels = grp_lv)]
  d
}
hc <- rbind(read_hc("m6A_highconf_counts_NEWPATH.tsv", "m6A"),
            read_hc("m5c_highconf_counts_NEWPATH.tsv", "m5C"))

# ---- left: HC70 trajectory, indexed to Sham (solves 10x scale gap) -------------
hc[, idx := 100 * N_hc70 / mean(N_hc70[group == "Sham"]), by = mod]
traj <- hc[, .(mean_idx = mean(idx), sd_idx = sd(idx),
               mean_n = mean(N_hc70)), by = .(mod, group)]
raw_lab <- traj[group %in% c("Sham", "Day1"),
                .(txt = sprintf("%.1fK→%.1fK",
                                mean_n[group == "Sham"] / 1000,
                                mean_n[group == "Day1"] / 1000)), by = mod]

p_l <- ggplot(traj, aes(x = group, y = mean_idx, color = mod, group = mod)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(data = hc, aes(x = group, y = idx, color = mod, group = mod),
             position = position_jitter(width = 0.12, height = 0, seed = 3),
             alpha = 0.4, size = 1.1, show.legend = FALSE) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.4) +
  geom_errorbar(aes(ymin = mean_idx - sd_idx, ymax = mean_idx + sd_idx),
                width = 0.18, linewidth = 0.6) +
  geom_text(data = traj[group == "Day1"], inherit.aes = FALSE,
            aes(x = ifelse(mod == "m6A", 2.76, 3.24),
                y = mean_idx + ifelse(mod == "m6A", 7.5, -7.5),
                color = mod, label = sprintf("%+.1f%%", mean_idx - 100),
                hjust = ifelse(mod == "m6A", 1, 0)),
            size = 3.4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = mod_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.1))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "HC sites (frac ≥ 0.7),\nSham = 100%",
       title = "C1  m6A stable, m5C genuinely lost at high confidence",
       subtitle = sprintf("HC70 site counts, indexed to Sham; mean ± SD (n=3)\nraw: m6A %s, m5C %s (Sham→Day1)",
                          raw_lab$txt[raw_lab$mod == "m6A"], raw_lab$txt[raw_lab$mod == "m5C"]),
       color = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8.5),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.2, color = "grey30"),
        plot.margin = margin(5.5, 44, 5.5, 5.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = c(1, 0))

# ---- right: %count-change vs fraction gate (from Part C bins) ------------------
read_bins <- function(f, mod) {
  b <- fread(file.path(dd, f))
  b[, .(N = sum(N)), by = .(sample, group, bin)][, mod := mod]
}
bins <- rbind(read_bins("m6a_threshold_motif_counts_NEWPATH.tsv", "m6A"),
              read_bins("m5c_threshold_bins_NEWPATH.tsv", "m5C"))
# cumulative N >= threshold per sample
setorder(bins, mod, sample, -bin)
bins[, cumN := cumsum(N), by = .(mod, sample)]
thr_dt <- bins[, .(cumN = mean(cumN)), by = .(mod, group, thr = bin / 10)]
wide <- dcast(thr_dt, mod + thr ~ group, value.var = "cumN")
wide[, pctD1 := 100 * (Day1 - Sham) / Sham]
wide <- wide[thr >= 0.1 & thr <= 0.9]

p_r <- ggplot(wide, aes(x = thr, y = pctD1, color = mod)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4, shape = 21, fill = "white", stroke = 1.3) +
  geom_text(data = wide[thr == 0.9], inherit.aes = FALSE,
            aes(x = 0.92, y = pctD1, label = mod, color = mod),
            size = 3.4, fontface = "bold", hjust = 0, show.legend = FALSE) +
  scale_color_manual(values = mod_colors) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  coord_cartesian(clip = "off") +
  labs(x = "Fraction gate (≥)", y = "Sham → Day1 count change (%)",
       title = "C2  The m6A 'loss' deflates with stringency",
       subtitle = "counts above each fraction gate (group means); m5C loss persists") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8.5),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.2, color = "grey30"),
        plot.margin = margin(5.5, 34, 5.5, 5.5),
        legend.position = "none")

p_out <- p_l + p_r + plot_layout(widths = c(1.15, 1))
ggsave(file.path(out_dir, "Fig2C_HC_trajectory.pdf"), p_out,
       width = 10, height = 4.1, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig2C_HC_trajectory.png"), p_out,
       width = 10, height = 4.1, dpi = 600)
cat("saved Fig2B_v2_HC_trajectory.{pdf,png}\n")
print(traj[, .(mod, group, mean_n = round(mean_n), mean_idx = round(mean_idx, 1))])
