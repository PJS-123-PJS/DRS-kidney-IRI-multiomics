# Fig8G v2 — Cross-cohort PT validation forest plot (REPLACES qualitative grid)
# ============================================================================
# 2026-07-26 — old G (hardcoded 4x3 word-grid) archived: cells were string
# literals with zero file reads; audit showed GSE184171/GSE252106 claims had
# NO analysis behind them, and "GSE184171 (Gerhardt 2021)" was a WRONG
# accession (true Gerhardt 2021 = GSE171417). v2 plots REAL re-analyses:
#   data/planB_cohort_metrics.csv  (server: planB_mouse_metrics.csv +
#                                   planB_human_metrics.csv, Welch t per cohort)
# Design: INJ_MODULE = positive control (PT injury axis must replicate);
# enzyme metrics shown HONESTLY incl. heterogeneous/phase-dependent effects.
# Zero hardcoding: every number read from the CSV.
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
  "INJ_MODULE"        = "PT injury module\n(HAVCR1/LCN2/VCAM1)",
  "NOP56"             = "NOP56",
  "M6A_WRITER_MODULE" = "m6A writers\n(METTL3/14/WTAP)",
  "METTL3"            = "METTL3",
  "FBL"               = "FBL"
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
         paste0("   ", contrast_lab[as.character(m$contrast)])),
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
       title = expression(bold("g")~~"Cross-cohort validation in PT cells (sc/snRNA-seq)"),
       subtitle = "4 independent cohorts re-analyzed from raw counts; injury module = positive control; Welch t, cell-level") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.5, face = "bold"),
    plot.subtitle = element_text(size = 6.3, color = "grey30"),
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

w <- 178; h <- 68
ggsave(file.path(OUT, "Fig8G_v2_CrossCohort_Forest.pdf"), p, width = w/25.4, height = h/25.4, device = cairo_pdf)
ggsave(file.path(OUT, "Fig8G_v2_CrossCohort_Forest.png"), p, width = w/25.4, height = h/25.4, dpi = 600)
cat(sprintf("Fig8G v2 saved: %d rows x %d metrics from %s\n", nrow(m), length(unique(m$metric)), CSV))
print(m[, .(cohort, contrast, metric, effect = round(effect, 4), lo = round(lo, 4), hi = round(hi, 4), sig, n_tot)])
