# Fig6D Nm Western blot densitometry — fills the Fig6D quantification placeholder
# Input : output/Fig6D_Nm_WB_densitometry_20260915.csv (python densitometry, 2026-09-15)
# Output: Fig6D_Nm_WB_densitometry_20260915.{pdf,png}
# Caveats on the plot: single blot (n=1, descriptive only); strong bands partially
# saturated -> Day3/Day7 folds are lower bounds. Style: canonical group colors, Arial 7pt.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

OUT <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
dt <- fread(file.path(OUT, "Fig6D_Nm_WB_densitometry_20260915.csv"))
dt <- dt[protein != "beta-actin"]
dt[, lane := factor(lane, levels = c("Sham","IS","Day1","Day3","Day7"))]
dt[, protein := factor(protein, levels = c("Fibrillarin","Nop56"))]

GCOL <- c(Sham="#4A4A4A", IS="#8A8A8A", Day1="#D6404C", Day3="#E8753A", Day7="#2E7D6F")

p <- ggplot(dt, aes(lane, rel, fill = lane)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#B0B0B0", linewidth = 0.3) +
  geom_col(width = 0.62, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", rel)), vjust = -0.4, size = 1.9, family = "Arial") +
  facet_wrap(~ protein, nrow = 1) +
  scale_fill_manual(values = GCOL) +
  scale_y_continuous(limits = c(0, 2.8), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Relative band intensity\n(target / β-actin, Sham = 1)",
       title = "Nm machinery Western blot densitometry (Fig 6D)",
       subtitle = "single blot, n = 1 — descriptive only\nstrong bands partially saturated (Day3/Day7 = lower bounds)") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(), legend.position = "none",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7.5),
        plot.subtitle = element_text(size = 5.8, color = "grey40"))

ggsave(file.path(OUT, "Fig6D_Nm_WB_densitometry_20260915.pdf"), p,
       width = 3.5, height = 1.75, device = cairo_pdf)
ggsave(file.path(OUT, "Fig6D_Nm_WB_densitometry_20260915.png"), p,
       width = 3.5, height = 1.75, dpi = 600, type = "cairo")
cat("DONE\n")
