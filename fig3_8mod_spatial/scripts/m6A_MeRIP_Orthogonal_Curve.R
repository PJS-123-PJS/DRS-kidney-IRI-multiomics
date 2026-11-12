# =============================================================================
# MeRIP orthogonal validation curve v2 (NAR reliability package, piece 2)
# 2026-07-22
# Input : data/public_GSE224289/merip_persite_records.tsv  (per-site window
#         query: sample, group, set(tier|control), frac, log2fc_mean,
#         input_mean, supported)
# Output: output/m6A_MeRIP_Orthogonal_Curve.{pdf,png}
# Two analysis modes side by side:
#   raw        = all queried DRS sites
#   detectable = only sites whose ±100bp window has mean input RPKM >= 1
#                (i.e. transcripts the MeRIP experiment could have detected)
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

data_dir <- "../data/public_GSE224289"
out_dir  <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(Sham = "#E64B35", Day1 = "#4DBBD5")
tier_lv <- c("0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-0.9", ">=0.9", "control")

dt <- fread(file.path(data_dir, "merip_persite_records.tsv"))
dt[, set := factor(set, levels = tier_lv)]
dt[, group := factor(group, levels = c("Sham", "Day1"))]
dt[, mode := ifelse(input_mean >= 1, "detectable", "low_input")]

# ---- summarise: per (analysis, group, sample, tier) ---------------------------
z95 <- 1.96
mk_summ <- function(d) {
  d[, .(n = .N, sup = sum(supported)), by = .(group, sample, set)]
}
summ <- rbind(
  mk_summ(dt)[, analysis := "All queried sites"],
  mk_summ(dt[input_mean >= 1])[, analysis := "MeRIP-detectable (input>=1)"])
summ[, pct := 100 * sup / n]
summ[, phat := sup / n]
summ[, `:=`(center = (phat + z95^2/(2*n)) / (1 + z95^2/n),
            half = z95 * sqrt(phat*(1 - phat)/n + z95^2/(4*n^2)) / (1 + z95^2/n))]
summ[, `:=`(lo = 100 * (center - half), hi = 100 * (center + half))]

pooled <- summ[, .(n = sum(n), sup = sum(sup)), by = .(analysis, group, set)]
pooled[, pct := 100 * sup / n]
ctrl_pool <- pooled[set == "control"]
main_pool <- pooled[set != "control"]
main_ps   <- summ[set != "control"]

mode_lv <- c("All queried sites", "MeRIP-detectable (input>=1)")
main_pool[, analysis := factor(analysis, levels = mode_lv)]
main_ps[, analysis := factor(analysis, levels = mode_lv)]
ctrl_pool[, analysis := factor(analysis, levels = mode_lv)]

# ---- A: support curves --------------------------------------------------------
p_a <- ggplot() +
  geom_hline(data = ctrl_pool, aes(yintercept = pct, color = group),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_line(data = main_ps,
            aes(x = set, y = pct, color = group, group = interaction(group, sample)),
            alpha = 0.22, linewidth = 0.5, show.legend = FALSE) +
  geom_line(data = main_pool, aes(x = set, y = pct, color = group, group = group),
            linewidth = 1.2) +
  geom_point(data = main_pool, aes(x = set, y = pct, color = group, group = group),
             size = 2.6, shape = 21, fill = "white", stroke = 1.3) +
  facet_wrap(~analysis) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  labs(x = "DRS m6A fraction tier", y = "MeRIP-supported sites (%)",
       title = "A  Public MeRIP-seq concordance rises with DRS stoichiometry",
       subtitle = paste0("window ±100 bp around each DRS A; supported = IP/input≥2 in ≥2/3 WT replicates (GSE224289 WT, IRI-24h, GRCm38)\n",
                         "faint = per-sample (n=3), bold = pooled; dashed = low-fraction control (0.02–0.1)"),
       color = "DRS group") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.0, color = "grey30"),
        strip.text = element_text(size = 8.5, face = "bold"),
        legend.position = c(0.15, 0.82),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

# ---- B: HC vs control bars + Fisher, both analyses ----------------------------
hc_id <- c("0.7-0.9", ">=0.9")
bar_dt <- rbind(
  summ[set %in% hc_id, .(n = sum(n), sup = sum(sup), set = "HC (frac>=0.7)"),
       by = .(analysis, group)],
  summ[set == "control", .(n = sum(n), sup = sum(sup), set = "control (0.02-0.1)"),
       by = .(analysis, group)])
bar_dt[, pct := 100 * sup / n]
bar_dt[, set := factor(set, levels = c("control (0.02-0.1)", "HC (frac>=0.7)"))]
bar_dt[, analysis := factor(analysis, levels = mode_lv)]

fisher_dt <- bar_dt[, {
  h <- .SD[set == "HC (frac>=0.7)"]
  c0 <- .SD[set == "control (0.02-0.1)"]
  ft <- fisher.test(matrix(c(h$sup, h$n - h$sup, c0$sup, c0$n - c0$sup),
                           nrow = 2, byrow = TRUE))
  .(OR = unname(ft$estimate), p = ft$p.value)
}, by = .(analysis, group)]
fisher_dt[, label := sprintf("OR=%.1f\np=%.1e", OR, p)]
# place each group's Fisher label above its own HC bar; stagger vertically
# (the two bars are near-identical height -> same-height labels would collide)
fisher_dt <- merge(fisher_dt,
                   bar_dt[set == "HC (frac>=0.7)", .(analysis, group, pct, set)],
                   by = c("analysis", "group"))
fisher_dt[, y_lab := pct + ifelse(group == "Sham", 2.4, 5.6)]

p_b <- ggplot(bar_dt, aes(x = set, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", pct), group = group),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.6, fontface = "bold") +
  geom_text(data = fisher_dt, inherit.aes = FALSE,
            aes(x = set, y = y_lab, label = label, group = group, color = group),
            position = position_dodge(width = 0.7), vjust = 0,
            size = 2.4, fontface = "bold", lineheight = 0.9, show.legend = FALSE) +
  facet_wrap(~analysis) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.34))) +
  labs(x = NULL, y = "MeRIP-supported (%)",
       title = "B  HC DRS sites vs low-fraction control",
       subtitle = "Fisher exact test, pooled n=3 replicates per group",
       fill = "DRS group") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.0, color = "grey30"),
        strip.text = element_text(size = 8.5, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

p_out <- p_a + p_b + plot_layout(widths = c(1.3, 1))
ggsave(file.path(out_dir, "m6A_MeRIP_Orthogonal_Curve.pdf"), p_out,
       width = 11.5, height = 4.4, device = cairo_pdf)
ggsave(file.path(out_dir, "m6A_MeRIP_Orthogonal_Curve.png"), p_out,
       width = 11.5, height = 4.4, dpi = 600)
cat("saved m6A_MeRIP_Orthogonal_Curve.{pdf,png}\n")
print(fisher_dt)
