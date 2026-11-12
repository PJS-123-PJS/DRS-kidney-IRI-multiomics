# =============================================================================
# FigS2H: Paired stoichiometry density of Sham high-confidence m6A sites
# 2026-09-23 — canonical re-render on the NewPath per-site context.
#
# Supersedes the 2026-07-22 render (generator archived 2026-08-19): that render
# paired sites by depth>=20 in all 6 samples (n=8,108 DRACH+ / 1,031 DRACH-).
# Pairing here = detected 3/3 in BOTH conditions (n_det_NC==3 & n_det_Day1==3),
# identical to the ModVerdict dual-anchor arbitration and manuscript text:
#   DRACH+ HC70: n=5,174, median delta = +2.84pp (wilcox p = 7.1e-226)
#   DRACH- HC70: n=379,   median delta = -30.7pp (reference collapsing to ~0)
# Fidelity gate below hard-asserts these n values against
# data/m6A_paired_stoichiometry_NEWPATH.tsv (Part B of m6A_NewPath_Reliability.R).
#
# Palette unified 2026-09-23: DRACH+ #D6404C / DRACH- #8A8A8A (project canonical).
# Input: data/site_ctx_nc_day1_NEWPATH.tsv.gz (NewPath is pre-filtered depth>=20)
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

data_dir <- "../data"
out_dir  <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ctx <- fread(file.path(data_dir, "site_ctx_nc_day1_NEWPATH.tsv.gz"))
paired <- ctx[n_det_NC == 3 & n_det_Day1 == 3]
hc <- paired[mean_frac_Sham >= 0.7]
hc[, delta := mean_frac_Day1 - mean_frac_Sham]
hc[, cls := fifelse(is_drach, "DRACH+ (certain m6A)", "DRACH- (reference)")]

stats <- hc[, .(n = .N, med = median(delta)), by = cls]
print(stats)

# ---- fidelity gate: must reproduce the canonical Part B paired stats --------
n_plus  <- stats[grepl("^DRACH\\+", cls), n]
n_minus <- stats[grepl("^DRACH-",  cls), n]
stopifnot(n_plus == 5174L, n_minus == 379L)
med_plus  <- stats[grepl("^DRACH\\+", cls), med]
med_minus <- stats[grepl("^DRACH-",  cls), med]
stopifnot(abs(med_plus - 0.02836667) < 1e-4, abs(med_minus - (-0.3074667)) < 1e-4)
cat("fidelity gate PASS: n=5,174/379, medians +2.84pp/-30.7pp == m6A_paired_stoichiometry_NEWPATH.tsv\n")

pal <- c("DRACH+ (certain m6A)" = "#D6404C", "DRACH- (reference)" = "#8A8A8A")

p <- ggplot(hc, aes(x = delta, fill = cls, color = cls)) +
  geom_density(alpha = 0.30, linewidth = 0.9) +
  geom_vline(data = stats, aes(xintercept = med, color = cls),
             linetype = "dotted", linewidth = 0.8, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_color_manual(values = pal, guide = "none") +
  labs(x = expression(Delta * " m6A fraction (Day1 − Sham), paired sites"),
       y = "Density",
       title = "Paired stoichiometry of Sham high-confidence m6A sites",
       subtitle = sprintf("sites with mean Sham fraction ≥ 0.7, detected 3/3 in both conditions (n = %s DRACH+, %s DRACH−); dotted = median",
                          format(n_plus, big.mark = ","), format(n_minus, big.mark = ","))) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.4, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = c(0.28, 0.85),
        legend.background = element_blank(),
        legend.text = element_text(size = 8))

ggsave(file.path(out_dir, "m6A_Paired_Stoichiometry.pdf"), p,
       width = 6.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(out_dir, "m6A_Paired_Stoichiometry.png"), p,
       width = 6.5, height = 4.2, dpi = 600)
cat("Saved: m6A_Paired_Stoichiometry.{pdf,png}\n")
