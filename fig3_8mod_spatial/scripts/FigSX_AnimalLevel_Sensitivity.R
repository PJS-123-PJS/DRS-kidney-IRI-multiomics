# =============================================================================
# Per-animal sensitivity summary figure (supplement candidate)
# 2026-09-27 — one point per animal: median per-site delta (vs other-group
# mean) across each caller's verdict-defining anchor set.
# Day1 animals should be concordant and separated from Sham animals.
# Canonical palette: Sham #4A4A4A / Day1 #D6404C. Arial 7pt, cairo.
# Output: output/FigSX_AnimalLevel_Sensitivity.{pdf,png}
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); 
})

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
out_dir  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"

a7  <- fread(file.path(data_dir, "animal_level_sensitivity_7callers.tsv"))
am6 <- fread(file.path(data_dir, "animal_level_sensitivity_m6a.tsv"))
all <- rbind(a7, am6, fill = TRUE)

# verdict-defining anchor per caller (the anchor quoted in the main text)
key_sets <- list(
  "m6A DRACH+" = c("SH-HC70", "Real gain (mild)"),
  "m5C"        = c("SH-HC50", "Dominant loss"),
  "PSU"        = c("SH-HC70", "Stable"),
  "inosine"    = c("D1-HC50", "Dominant gain"),
  "Nm-Um"      = c("D1-HC50", "Dominant gain (largest)"),
  "Nm-Am"      = c("D1-HC50", "Real gain"),
  "Nm-Cm"      = c("D1-HC50", "Real gain"),
  "Nm-Gm"      = c("SH-HC50", "Dominant loss")
)

rows <- list()
for (k in names(key_sets)) {
  anch <- key_sets[[k]][1]; verdict <- key_sets[[k]][2]
  r <- all[mod == k & anchor == anch]
  if (!nrow(r)) next
  d1 <- as.numeric(strsplit(r$Day1_animal_median_deltas_pp, ";")[[1]])
  sh <- as.numeric(strsplit(r$Sham_animal_median_deltas_pp, ";")[[1]])
  rows[[k]] <- data.table(
    caller = paste0(k, "\n", verdict, " (n=", r$n_sites, ")"),
    group = rep(c("Day1", "Sham"), each = 3),
    animal = rep(c("s1", "s2", "s3"), 2),
    delta = c(d1, sh))
}
dt <- rbindlist(rows)
dt[, group := factor(group, levels = c("Sham", "Day1"))]

pal <- c("Sham" = "#4A4A4A", "Day1" = "#D6404C")

p <- ggplot(dt, aes(group, delta, color = group)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey70", linewidth = 0.3) +
  geom_point(size = 1.6, position = position_jitter(width = 0.12, seed = 7)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45,
               color = "black", linewidth = 0.3, fatten = 1.2) +
  facet_wrap(~caller, nrow = 2, scales = "free_y") +
  scale_color_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Per-animal median Δ stoichiometry (pp)",
       title = "Animal-level sensitivity of the arbitration verdicts",
       subtitle = "One point per animal (n = 3 per group); median per-site Sham→Day1 change across each verdict-defining anchor set; bars = group medians") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 6.5),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "grey30"),
        axis.title = element_text(size = 7))

ggsave(file.path(out_dir, "FigSX_AnimalLevel_Sensitivity.pdf"), p,
       width = 6.8, height = 4.4, device = cairo_pdf)
ggsave(file.path(out_dir, "FigSX_AnimalLevel_Sensitivity.png"), p,
       width = 6.8, height = 4.4, dpi = 600)
cat("written:", file.path(out_dir, "FigSX_AnimalLevel_Sensitivity.{pdf,png}"), "\n")
