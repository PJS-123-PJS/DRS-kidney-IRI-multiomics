# =============================================================================
# 5-mod paired-arbitration verdict figure (replaces all "wave" narrative)
# 2026-07-22
# Panel A: paired stoichiometry delta (Sham->Day1) per caller, BOTH anchors
#   (Day1-HC50 filled circle, Sham-HC50 open triangle) with wilcox CI.
#   Self-diagnosing: same-sign = real; symmetric straddle = stable;
#   opposite with |SH|>|D1| = dominant loss. Verdict label DERIVED by rule.
# Panel B: anchor-free HC70 site-count trajectory (indexed to Sham), all
#   callers, two facets (single-base callers / Nm bases).
# Mock-testable: MODVERDICT_DATA_DIR pointing to mock dir with same files.
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

dd <- Sys.getenv("MODVERDICT_DATA_DIR", "../data")
out_dir <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load -------------------------------------------------------------------
st <- fread(file.path(dd, "modverdict_paired_stats.tsv"))
st <- st[anchor %in% c("D1-HC50", "SH-HC50") & n >= 10]

caller_lv <- c("m6A", "m5C", "PSU", "inosine",
               "Nm-Um", "Nm-Cm", "Nm-Am", "Nm-Gm")
class_map <- c("m6A" = "m6A", "m6A DRACH+" = "m6A",
               "m5C" = "m5C", "m5C gate+" = "m5C",
               "PSU" = "PSU", "inosine" = "inosine",
               "Nm-Um" = "Nm (2'-O-Me)", "Nm-Am" = "Nm (2'-O-Me)",
               "Nm-Cm" = "Nm (2'-O-Me)", "Nm-Gm" = "Nm (2'-O-Me)")
cls_colors <- c("m6A" = "#D6404C", "m5C" = "#E8753A", "PSU" = "#3C5488",
                "inosine" = "#00A087", "Nm (2'-O-Me)" = "#6A3D9A")

st[, caller := factor(mod, levels = caller_lv)]
st[, cls := class_map[as.character(mod)]]
st[, anchor_f := factor(anchor, levels = c("SH-HC50", "D1-HC50"))]

# ---- Panel B data first: anchor-free HC70 counts (needed as verdict witness) -
cnt_files <- c("m6A" = "m6A_highconf_counts_NEWPATH.tsv",
               "m5C" = "m5c_highconf_counts_NEWPATH.tsv",
               "PSU" = "psu_highconf_counts_NEWPATH.tsv",
               "inosine" = "inosine_highconf_counts_NEWPATH.tsv",
               "Nm-Um" = "um_highconf_counts_NEWPATH.tsv",
               "Nm-Am" = "am_highconf_counts_NEWPATH.tsv",
               "Nm-Cm" = "cm_highconf_counts_NEWPATH.tsv",
               "Nm-Gm" = "gm_highconf_counts_NEWPATH.tsv")
cnt <- rbindlist(lapply(names(cnt_files), function(cl) {
  f <- file.path(dd, cnt_files[[cl]])
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d[, caller := cl]
  d
}))
grp_lv <- c("Sham", "IS", "Day1", "Day3", "Day7")
cnt[, group := factor(group, levels = grp_lv)]
# caller each row belongs to (gate subsets inherit their caller's counts)
cnt_change <- cnt[, .(pct_d1 = 100 * (mean(N_hc70[group == "Day1"]) /
                                      mean(N_hc70[group == "Sham"]) - 1)), by = caller]
caller_of <- c("m6A" = "m6A", "m6A DRACH+" = "m6A",
               "m5C" = "m5C", "m5C gate+" = "m5C",
               "PSU" = "PSU", "inosine" = "inosine",
               "Nm-Um" = "Nm-Um", "Nm-Am" = "Nm-Am",
               "Nm-Cm" = "Nm-Cm", "Nm-Gm" = "Nm-Gm")

# ---- rule-derived verdict per caller ----------------------------------------
# Opposite-sign anchors = regression-to-mean suspected. Third witness: the
# anchor-free HC70 count change. Small anchor effects (<5pp) with big count
# swings = sub-threshold-band churn (threshold-crossing inflation), i.e.
# stoichiometrically stable. Large anchor effects: count sign arbitrates.
verdict_dt <- st[, {
  d1 <- .SD[anchor == "D1-HC50"]; sh <- .SD[anchor == "SH-HC50"]
  ccnt <- cnt_change[caller == caller_of[mod]]$pct_d1
  d1_up <- d1$ci_lo_pp > 0; d1_dn <- d1$ci_hi_pp < 0
  sh_up <- sh$ci_lo_pp > 0; sh_dn <- sh$ci_hi_pp < 0
  v <- if (d1_up && sh_up) "Real gain"
  else if (d1_dn && sh_dn) "Real loss"
  else if ((d1_up && sh_dn) || (d1_dn && sh_up)) {
    mx <- max(abs(d1$est_pp), abs(sh$est_pp))
    if (mx < 5) "Stable"
    else if (length(ccnt) == 0 || abs(ccnt) < 15) "Stable (symmetric noise)"
    else if (ccnt < 0) ifelse(mx >= 10, "Dominant loss", "Modest loss (band)")
    else ifelse(mx >= 10, "Dominant gain", "Modest gain (band)")
  }
  else if (d1_up) "Real gain"
  else if (d1_dn) "Real loss"
  else if (sh_dn) "Loss (Sham anchor)"
  else if (sh_up) "Gain (Sham anchor)"
  else "Stable"
  .(verdict = v, count_pct_d1 = ccnt,
    lab = sprintf("%s\nHC70 counts %s%.0f%%", v, ifelse(ccnt >= 0, "+", ""), ccnt))
}, by = .(mod, cls)]
verdict_dt[, caller := factor(mod, levels = caller_lv)]
print(verdict_dt)

st <- merge(st, verdict_dt[, .(mod, verdict)], by = "mod")

p_a <- ggplot(st, aes(y = caller, color = cls)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo_pp, xmax = ci_hi_pp, group = anchor_f),
                 position = position_dodge(width = 0.55), height = 0.25, linewidth = 0.7) +
  geom_point(aes(x = est_pp, shape = anchor_f, fill = cls, group = anchor_f),
             position = position_dodge(width = 0.55), size = 2.8, stroke = 1.1) +
  geom_text(aes(x = est_pp, label = sprintf("n=%d", n), group = anchor_f),
            position = position_dodge(width = 0.55), vjust = 1.9,
            size = 2.3, color = "grey45", show.legend = FALSE) +
  geom_text(data = verdict_dt, inherit.aes = FALSE,
            aes(x = Inf, y = caller, label = lab, color = cls),
            hjust = -0.06, size = 2.7, fontface = "bold", lineheight = 0.95,
            show.legend = FALSE) +
  scale_shape_manual(values = c("SH-HC50" = 24, "D1-HC50" = 21),
                     labels = c("SH-HC50" = "Sham-anchored HC50",
                                "D1-HC50" = "Day1-anchored HC50")) +
  scale_color_manual(values = cls_colors, guide = "none") +
  scale_fill_manual(values = cls_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.30))) +
  coord_cartesian(clip = "off") +
  labs(x = "Paired stoichiometry change, Sham → Day1 (percentage points)",
       y = NULL,
       title = "A  Dual-anchor paired arbitration per modification",
       subtitle = "fully paired sites (3/3 reps each condition); points = pseudomedian Δ, whiskers = 95% CI (Wilcoxon)\nsame-sign anchors → real change; symmetric straddle of 0 → stable; opposite signs → HC70 count change arbitrates",
       shape = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8.5),
        plot.title = element_text(size = 10.5, face = "bold"),
        plot.subtitle = element_text(size = 7.0, color = "grey30"),
        plot.margin = margin(5.5, 96, 5.5, 5.5),
        legend.position = "top",
        legend.justification = c(0, 0),
        legend.text = element_text(size = 8))

# ---- Panel B: HC70 count trajectory (anchor-free) ----------------------------
cnt[, cls := ifelse(caller %in% c("m6A", "m5C"), caller,
             ifelse(caller %in% c("PSU", "inosine"), caller, "Nm (2'-O-Me)"))]
# caller-outlier flag: within caller x group, max/median > 10x = caller glitch
# (e.g. Day3s1 inosine 32,341 vs siblings 312/87). Excluded from mean, plotted
# as an x so the point is reported, not hidden.
cnt[, outlier := N_hc70 / median(N_hc70) > 10, by = .(caller, group)]
if (any(cnt$outlier)) {
  cat("caller-outlier samples excluded from trajectory means:\n")
  print(cnt[outlier == TRUE, .(caller, sample, group, N_hc70)])
}
cnt[, facet := ifelse(caller %in% c("m6A", "m5C", "PSU", "inosine"),
                      "Single-base callers", "Nm bases")]
cnt_ok <- cnt[outlier == FALSE]
cnt[, idx := 100 * N_hc70 / mean(N_hc70[group == "Sham" & outlier == FALSE]), by = caller]
cnt_ok[, idx := 100 * N_hc70 / mean(N_hc70[group == "Sham"]), by = caller]
traj <- cnt_ok[, .(mean_idx = mean(idx), sd_idx = sd(idx)), by = .(caller, cls, facet, group)]
traj[, facet := factor(facet, levels = c("Single-base callers", "Nm bases"))]
cnt[, facet := factor(facet, levels = c("Single-base callers", "Nm bases"))]
cnt_ok[, facet := factor(facet, levels = c("Single-base callers", "Nm bases"))]

caller_colors <- c("m6A" = "#D6404C", "m5C" = "#E8753A", "PSU" = "#3C5488",
                   "inosine" = "#00A087", "Nm-Um" = "#6A3D9A", "Nm-Cm" = "#9E77B5",
                   "Nm-Am" = "#C2A9D6", "Nm-Gm" = "#4B2673")

# cap axis from non-outlier range; off-scale outliers named in subtitle
y_cap <- ceiling(max(c(cnt_ok$idx, traj$mean_idx + traj$sd_idx), na.rm = TRUE) * 1.15 / 50) * 50
outl <- cnt[outlier == TRUE]
out_note <- if (nrow(outl) == 0) "mean ± SD (n=3), indexed to Sham" else
  paste(strwrap(paste0("mean ± SD (n=3), indexed to Sham; × = caller-outlier sample excluded from mean (off-scale): ",
         paste(sprintf("%s %s (idx %.0f%%)", outl$caller, outl$sample, outl$idx), collapse = ", ")),
        width = 78), collapse = "\n")

p_b <- ggplot(traj, aes(x = group, y = mean_idx, color = caller, group = caller)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_line(linewidth = 1.0) +
  geom_point(data = cnt_ok, aes(x = group, y = idx, color = caller, group = caller),
             position = position_jitter(width = 0.1, height = 0, seed = 3),
             alpha = 0.35, size = 0.9, show.legend = FALSE) +
  geom_point(data = cnt[outlier == TRUE],
             aes(x = group, y = idx, color = caller, group = caller),
             shape = 4, size = 2.2, stroke = 1.2, show.legend = FALSE) +
  geom_point(size = 2.4, shape = 21, fill = "white", stroke = 1.2) +
  geom_errorbar(aes(ymin = mean_idx - sd_idx, ymax = mean_idx + sd_idx),
                width = 0.15, linewidth = 0.5) +
  facet_wrap(~facet, ncol = 1) +
  scale_color_manual(values = caller_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.1))) +
  coord_cartesian(ylim = c(0, y_cap)) +
  labs(x = NULL, y = "HC70 sites (Sham = 100%)",
       title = "B  Anchor-free net change:\nHC site counts over time",
       subtitle = out_note,
       color = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.2, color = "grey30"),
        strip.text = element_text(size = 8.5, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))

p_out <- p_a + p_b + plot_layout(widths = c(1.45, 1))
ggsave(file.path(out_dir, "ModVerdict_Paired_Arbitration.pdf"), p_out,
       width = 12.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(out_dir, "ModVerdict_Paired_Arbitration.png"), p_out,
       width = 12.5, height = 5.6, dpi = 600)
cat("saved ModVerdict_Paired_Arbitration.{pdf,png}\n")
