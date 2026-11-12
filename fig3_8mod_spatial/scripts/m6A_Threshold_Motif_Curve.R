# =============================================================================
# m6A Threshold-Motif Curve (NAR reliability panel, 2026-07-21)
# Question: does DRACH% rise with the modification-fraction (stoichiometry)
# threshold? If high-fraction sites are ~90% DRACH, the "23.5% DRACH" of the
# lenient call set reflects low-stoichiometry noise, not false calling.
# Input: NEWPATH binned counts from canonical all_m6a.rds
#   (data/m6a_threshold_motif_counts_NEWPATH.tsv, 2026-07-22; schema bin/N/DRACH/sample/group,
#   bin k = fraction [k/10,(k+1)/10) ).
#   2026-08-19: switched from old-path m6a_threshold_motif_counts.tsv (superseded;
#   old render showed 21%/89%/bg7.5% vs canonical 34.1%/94.4%/bg9.9% — archived
#   supplement/figS1_qc/_archive_stale_20260819/). Background = bin 0 ([0,0.1)).
#   2026-09-23: palette unified to project canonical five-group colors
#   (Sham #4A4A4A / IS #8A8A8A / Day1 #D6404C / Day3 #E8753A / Day7 #2E7D6F,
#   same as wet-lab figures), replacing the earlier NPG mapping.
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

data_dir <- "../data"
out_dir  <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Canonical five-group palette (2026-09-23 unification; wet-lab figures同款)
group_map <- c(NCs1="Sham",NCs2="Sham",NCs3="Sham",
               Is1="IS",Is2="IS",Is3="IS",
               Day1s1="Day1",Day1s2="Day1",Day1s3="Day1",
               Day3s1="Day3",Day3s2="Day3",Day3s3="Day3",
               Day7s1="Day7",Day7s2="Day7",Day7s3="Day7")
group_levels <- c("Sham","IS","Day1","Day3","Day7")
group_colors <- c(Sham="#4A4A4A", IS="#8A8A8A", Day1="#D6404C",
                  Day3="#E8753A", Day7="#2E7D6F")

cnt <- fread(file.path(data_dir, "m6a_threshold_motif_counts_NEWPATH.tsv"))
cnt[, Group := factor(group, levels = group_levels)]
stopifnot(!any(is.na(cnt$Group)))

# ---- cumulative (>= t) stats from 0.1-wide bins -----------------------------
thresholds <- seq(0.1, 0.9, by = 0.1)
cum_dt <- rbindlist(lapply(thresholds, function(t) {
  sub <- cnt[bin >= t * 10 - 1e-9]
  sub[, .(N = sum(N), DRACH = sum(DRACH)), by = Group][, Threshold := t]
}))
cum_dt[, DRACH_pct := DRACH / N * 100]
cum_dt[, Group := factor(Group, levels = group_levels)]

# pooled (all 15 samples) for annotations
pooled <- cum_dt[, .(N = sum(N), DRACH = sum(DRACH)), by = Threshold]
pooled[, DRACH_pct := DRACH / N * 100]

# background: bin [0,0.1), >99% unmodified A positions
bg <- cnt[bin == 0, .(N = sum(N), DRACH = sum(DRACH))]
bg_pct <- bg$DRACH / bg$N * 100
cat(sprintf("Background (Fraction<0.1): DRACH = %.1f%% (n=%s)\n",
            bg_pct, format(bg$N, big.mark = ",")))
print(pooled[, .(Threshold, N, DRACH_pct = round(DRACH_pct, 1))])

key_t <- c(0.1, 0.5, 0.7)
# round() guards against seq() floating-point representation mismatch
key_rows <- pooled[round(Threshold, 1) %in% key_t & N > 0]
key_rows[, lab := sprintf("≥%.1f: %.0f%%", round(Threshold, 1), DRACH_pct)]
# leftmost label would clip the y-axis -> nudge right+below; others above center
key_rows[, nudge_x := fifelse(round(Threshold, 1) == 0.1, 0.045, 0)]
key_rows[, nudge_y := fifelse(round(Threshold, 1) == 0.1, -9, 6)]
key_rows[, hj      := fifelse(round(Threshold, 1) == 0.1, 0, 0.5)]

subtitle_txt <- sprintf("Lenient call set (≥0.1) %.0f%% DRACH; high-confidence (≥0.7) %.0f%%; noise floor (Fraction<0.1) %.1f%%",
                        pooled[round(Threshold, 1) == 0.1]$DRACH_pct,
                        pooled[round(Threshold, 1) == 0.7]$DRACH_pct, bg_pct)

# ---- sub-panel A: DRACH% vs threshold ---------------------------------------
p_a <- ggplot(cum_dt, aes(x = Threshold, y = DRACH_pct, color = Group)) +
  geom_hline(yintercept = bg_pct, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  annotate("text", x = 0.55, y = bg_pct + 4,
           label = sprintf("noise floor (Fraction<0.1) %.1f%%", bg_pct), color = "grey40",
           size = 2.8, hjust = 0, family = "Arial") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  annotate("point", x = key_rows$Threshold, y = key_rows$DRACH_pct,
           shape = 21, size = 3.2, stroke = 1.1, fill = "white", color = "black") +
  geom_text(data = key_rows,
            aes(x = Threshold + nudge_x, y = DRACH_pct + nudge_y, label = lab, hjust = hj),
            inherit.aes = FALSE, size = 2.9, fontface = "bold", family = "Arial") +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = thresholds) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "m6A fraction threshold (≥)", y = "DRACH motif (%)",
       title = "Stoichiometry stratification rescues motif specificity",
       subtitle = subtitle_txt,
       color = NULL) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.4, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = c(0.25, 0.72),
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 7.5),
        legend.background = element_blank())

# ---- sub-panel B: sites retained vs threshold (pooled, log y) ----------------
p_b <- ggplot(pooled, aes(x = Threshold, y = N)) +
  geom_col(fill = "#8A8A8A", alpha = 0.85, width = 0.07) +
  geom_text(aes(label = scales::comma(N)), vjust = -0.4, size = 2.4,
            family = "Arial", color = "grey30") +
  scale_x_continuous(breaks = thresholds) +
  scale_y_log10(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "m6A fraction threshold (≥)", y = "Sites retained (log10)",
       title = "High-confidence yield") +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.4, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        plot.title = element_text(size = 10, face = "bold"))

p_full <- p_a + p_b + plot_layout(widths = c(0.62, 0.38))

ggsave(file.path(out_dir, "m6A_Threshold_Motif_Curve.pdf"), p_full,
       width = 9.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(out_dir, "m6A_Threshold_Motif_Curve.png"), p_full,
       width = 9.5, height = 4.2, dpi = 600)

cat("Saved: m6A_Threshold_Motif_Curve.{pdf,png}\n")
