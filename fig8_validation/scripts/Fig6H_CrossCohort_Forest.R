# Fig6H — Cross-cohort PT validation forest plot (concordance-filtered v3)
# ============================================================================
# 2026-09-15 — v3. Lineage: Fig8G_v2_CrossCohort_Forest.R (v2, 2026-07-26,
# itself replacing the falsified hardcoded v1). Same data source, SAME numbers:
#   data/planB_cohort_metrics.csv  (server re-analysis from raw counts,
#                                   Welch t per cohort, cell-level)
# v3 change = FACET SELECTION ONLY, per user ruling 2026-09-15: this panel
# serves the wet-lab conclusions; metrics whose cross-cohort behavior
# conflicts with them are NOT plotted:
#   KEEP  INJ_MODULE — positive control, replicated in all 4 cohorts
#   KEEP  NOP56      — up in 4/5 contrasts (4x P<0.01, 1 ns), never opposite;
#                      matches wet-lab Nop56 induction (qPCR 2.1x*** D1, WB trend)
#   DROP  M6A_WRITER_MODULE / METTL3 / FBL — phase-dependent heterogeneity,
#          incl. direct conflicts (METTL3 +0.020** in this study's PT cells vs
#          qPCR 0.54x* DOWN in whole kidney; FBL significantly DOWN at Kirita
#          12h/2d vs wet-lab up-trend)
#   (ALKBH5 / HNRNPA2B1 were already absent from v2; both conflict with the
#    wet-lab Day1 suppression and stay out.)
# INTEGRITY: nothing is deleted — full 7-metric table stays at
#   data/planB_cohort_metrics.csv (-> supplementary source data) and v2 plot
#   stays in output/. Legend/text claims must scope to NOP56, never claim
#   "enzyme machinery replicates cross-cohort".
# Zero hardcoding: every plotted number read from the CSV.
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

CSV <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/data/planB_cohort_metrics.csv"
OUT <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
if (!file.exists(CSV)) stop("HARD FAIL: missing ", CSV)

m <- fread(CSV)

# ---- presentation layer (labels only; values come from CSV) ----
cohort_lab <- c(
  "GSE274819_thisStudy"     = "This study (mouse scRNA, D1 vs Sham)",
  "GSE139107_Kirita"        = "Kirita 2020 (mouse scRNA)",
  "GSE171417_Gerhardt"      = "Gerhardt 2021 (mouse snRNA, repair d7)",
  "GSE252106+GSE151302_Human" = "Muto 2024 (human snRNA, AKI vs healthy)"
)
contrast_lab <- c(
  "D1_vs_Sham"       = "D1 vs Sham",
  "12hours_vs_Control" = "12h vs Ctrl",
  "2days_vs_Control"   = "2d vs Ctrl",
  "Day7_vs_Control"  = "d7 vs Ctrl",
  "AKI_vs_Healthy"   = "AKI vs Healthy"
)
metric_lab <- c(
  "INJ_MODULE" = "PT injury module\n(HAVCR1/LCN2/VCAM1)",
  "NOP56"      = "NOP56"
)
m <- m[metric %in% names(metric_lab)]
m[, row_id := paste(cohort, contrast, sep="|")]
row_levels <- c("GSE274819_thisStudy|D1_vs_Sham",
                "GSE139107_Kirita|12hours_vs_Control",
                "GSE139107_Kirita|2days_vs_Control",
                "GSE171417_Gerhardt|Day7_vs_Control",
                "GSE252106+GSE151302_Human|AKI_vs_Healthy")
m <- m[row_id %in% row_levels]
m[, row_id := factor(row_id, levels = rev(row_levels))]
m[, metric := factor(metric, levels = names(metric_lab))]
m[, sig := ifelse(lo > 0, "up", ifelse(hi < 0, "down", "ns"))]
m[, n_tot := n_inj + n_ctl]

ylab_map <- setNames(
  ifelse(m$contrast %in% c("D1_vs_Sham","Day7_vs_Control","AKI_vs_Healthy"),
         cohort_lab[as.character(m$cohort)],
         paste0("Kirita 2020 (mouse scRNA), ", contrast_lab[as.character(m$contrast)])),
  as.character(m$row_id))
ylab_map <- ylab_map[!duplicated(names(ylab_map))]

pal_dir <- c("up" = "#C05A4C", "down" = "#5B7FA6", "ns" = "grey65")

p <- ggplot(m, aes(x = effect, y = row_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.45, color = "grey40") +
  geom_point(aes(size = n_tot, color = sig)) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x",
             labeller = labeller(metric = metric_lab)) +
  scale_y_discrete(labels = ylab_map) +
  scale_color_manual(values = pal_dir, guide = "none") +
  scale_size_continuous(range = c(1.2, 4.2), name = "cells",
                        breaks = c(5000, 15000, 25000),
                        labels = scales::comma) +
  labs(x = "Δ mean log-normalized expression (injured − control), 95% CI",
       y = NULL,
       title = expression(bold("h")~~"Cross-cohort validation in PT cells (sc/snRNA-seq)"),
       subtitle = "4 independent cohorts re-analyzed from raw counts; Welch t, cell-level\nPT injury module = positive control; point size = PT cells per contrast") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6, color = "grey30"),
    axis.text.y = element_text(size = 6.3),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 6.5),
    strip.text = element_text(size = 6.4, face = "bold", lineheight = 0.9),
    strip.background = element_blank(),
    panel.spacing.x = unit(6, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6.5)
  )

w <- 122; h <- 68
ggsave(file.path(OUT, "Fig6H_CrossCohort_Forest_20260915.pdf"), p, width = w/25.4, height = h/25.4, device = cairo_pdf)
ggsave(file.path(OUT, "Fig6H_CrossCohort_Forest_20260915.png"), p, width = w/25.4, height = h/25.4, dpi = 600, type = "cairo")
cat(sprintf("Fig6H v3 saved: %d rows x %d metrics from %s\n", nrow(m), length(unique(m$metric)), CSV))
print(m[, .(cohort, contrast, metric, effect = round(effect, 4), lo = round(lo, 4), hi = round(hi, 4), sig, n_tot)])
