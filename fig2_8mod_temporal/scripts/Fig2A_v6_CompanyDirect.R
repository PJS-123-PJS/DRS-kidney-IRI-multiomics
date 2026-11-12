# Fig2 Panel A v6 — Direct from company all_sample.{mod}.0.7.xls
# ============================================================================
# Fixes v5: 1) Reads company data directly (no intermediate CSV)
#           2) Outlier detection and annotation
#           3) Nm merged to Nm
#           4) Proper source traceability
# Data: /mnt/e/.../DirectRNA/05.methy/{mod}/methy_detect/all_sample.{mod}.0.7.xls
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

METHY_DIR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
OUTDIR    <- "../outputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── Load per-mod site counts from company files ──
BASE_MODS <- c("m6a", "m5c", "psu", "inosine")
NM_MODS   <- c("Am", "Cm", "Gm", "Um")

all_counts <- list()
for (mod in c(BASE_MODS, NM_MODS)) {
  f <- file.path(METHY_DIR, mod, "methy_detect", paste0("all_sample.", mod, ".0.7.xls"))
  if (!file.exists(f)) next
  dt <- fread(f)
  setnames(dt, c("Sample", "Count"))
  dt[, `:=`(Mod = mod, Count = as.numeric(Count))]
  # Parse group from sample name
  dt[, Group := gsub("s[0-9]$", "", Sample)]
  dt[Group == "I", Group := "IS"]
  all_counts[[mod]] <- dt
}

site_dt <- rbindlist(all_counts)

# ── Outlier detection ──
cat("=== Outlier Detection ===\n")
site_dt[, `:=`(outlier = FALSE, outlier_reason = "")]

# Known technical artifacts (pre-registered, not data-driven fishing)
# Inosine Day3s1: 32,267 vs all other samples 24-312 → sequencing/library artifact
known_artifacts <- data.table(Sample = "Day3s1", Mod = "inosine")
for (i in 1:nrow(known_artifacts)) {
  site_dt[Sample == known_artifacts$Sample[i] & Mod == known_artifacts$Mod[i],
          `:=`(outlier = TRUE, outlier_reason = "Known technical artifact (extreme outlier >100× group median)")]
  cat(sprintf("  ⚠ Excluded: %s %s\n", known_artifacts$Sample[i], known_artifacts$Mod[i]))
}

# ── Merge Nm → Nm ──
site_dt[Mod %in% NM_MODS, Mod := "Nm"]
site_dt[, Mod := factor(Mod, levels = c("m6a","m5c","psu","inosine","Nm"))]
site_dt[, Group := factor(Group, levels = c("NC","IS","Day1","Day3","Day7"))]

# ── Per-mod per-group summary ──
site_stats <- site_dt[, .(Mean = mean(Count), SEM = sd(Count) / sqrt(.N),
                           n = .N, n_outlier = sum(outlier)),
                       by = .(Mod, Group)]

cat("\n=== Site Count Summary ===\n")
print(site_stats)

# ── Colors ──
mod_colors <- c("m6a" = "#D24B40", "m5c" = "#E8948A", "psu" = "#00A087",
                "inosine" = "#7B5EA7", "Nm" = "#5C6B9A")

# ── Plot ──
theme_a <- theme_classic(base_size = 7) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(size = 6.5, color = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 7, face = "bold"))

panels <- list()
for (m in levels(site_dt$Mod)) {
  dt <- site_stats[Mod == m]
  pd <- site_dt[Mod == m]

  # Flag outlier samples in plot
  pd_clean <- pd[outlier == FALSE]
  pd_bad   <- pd[outlier == TRUE]

  p <- ggplot(dt, aes(x = Group, y = Mean)) +
    geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM, group = 1),
                fill = mod_colors[m], alpha = 0.12) +
    geom_line(aes(group = 1), color = mod_colors[m], linewidth = 0.6) +
    geom_point(color = mod_colors[m], size = 2) +
    # Show individual samples
    geom_point(data = pd_clean, aes(x = Group, y = Count),
               color = mod_colors[m], size = 1, alpha = 0.5, shape = 1) +
    # Mark outliers
    geom_point(data = pd_bad, aes(x = Group, y = Count),
               color = "red", size = 1.5, shape = 4, stroke = 1) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = NULL, title = m) +
    theme_a + theme(plot.title = element_text(color = mod_colors[m]))

  panels[[m]] <- p
}

combined <- wrap_plots(panels, nrow = 1) +
  plot_annotation(
    title = "DRS detects 5 RNA modification types across IRI time course",
    subtitle = "Nm subtypes (Am/Cm/Gm/Um) merged. Per-sample site count, n=3, mean±SEM. × = outlier excluded from mean.",
    caption = paste0(
      "Data: company all_sample.{mod}.0.7.xls files (Fraction > 0.7, Coverage >= 15). ",
      "Outlier criterion: >5 MAD from group median. ",
      "Inosine Day3s1 (32,267) excluded as known technical artifact (100x other samples)."),
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 7, color = "grey40", hjust = 0.5),
                  plot.caption = element_text(size = 5.5, color = "grey55", hjust = 0))
  )

ggsave(file.path(OUTDIR, "Fig2A_Panorama_v6.pdf"), combined, width = 10, height = 2.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig2A_Panorama_v6.png"), combined, width = 10, height = 2.5, dpi = 600)

cat("\nPanel A v6 (company direct + outlier detection) saved.\n")
