# Fig6C m6A WB densitometry — 3 biological replicates, publication version
# SUPERSEDES the single-blot version (Fig6C_m6A_WB_densitometry_20260915.*, archive):
# run1 re-quantified with main-band windows for ALKBH5 (consistency with run2/3).
# Input : Fig6C_m6A_WB_densitometry_3reps_long_20260915.csv
# Stats : one-way ANOVA per protein + Tukey HSD (vs Sham reported)
# Output: Fig6C_m6A_WB_densitometry_3reps_20260915.{pdf,png} + _stats.csv
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

OUT <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
dt <- fread(file.path(OUT, "Fig6C_m6A_WB_densitometry_3reps_long_20260915.csv"))
dt[, lane := factor(lane, levels = c("Sham","IS","Day1","Day3","Day7"))]
dt[, protein := factor(protein, levels = c("Mettl3","Hnrnpa2b1","Alkbh5"))]

st <- rbindlist(lapply(split(dt, dt$protein), function(d){
  fit <- aov(rel ~ lane, data = d)
  tk <- TukeyHSD(fit)$lane
  keep <- grep("Sham$", rownames(tk))
  data.table(protein = as.character(d$protein[1]),
             comparison = rownames(tk)[keep],
             diff = tk[keep, "diff"], p_adj = tk[keep, "p adj"],
             anova_p = summary(fit)[[1]]["lane", "Pr(>F)"])
}))
st[, group := factor(sub("-Sham", "", comparison), levels = c("IS","Day1","Day3","Day7"))]
st[, star := fcase(p_adj < 0.001, "***", p_adj < 0.01, "**", p_adj < 0.05, "*", default = "")]
fwrite(st, file.path(OUT, "Fig6C_m6A_WB_densitometry_3reps_stats_20260915.csv"))
print(st)

sm <- dt[, .(m = mean(rel), s = sd(rel)), by = .(protein, lane)]
starpos <- merge(st[star != ""], sm, by.x = c("protein","group"), by.y = c("protein","lane"))
rp <- dt[, .(rmax = max(rel)), by = .(protein, lane)]
starpos <- merge(starpos, rp, by.x = c("protein","group"), by.y = c("protein","lane"))
starpos[, y := pmax(m + s, rmax) + 0.10]
# merge coerces protein to character -> facet order lost; re-factor every layer frame
PLEV <- c("Mettl3","Hnrnpa2b1","Alkbh5")
starpos[, protein := factor(protein, levels = PLEV)]
sm[, protein := factor(protein, levels = PLEV)]

GCOL <- c(Sham="#4A4A4A", IS="#8A8A8A", Day1="#D6404C", Day3="#E8753A", Day7="#2E7D6F")
p <- ggplot(dt, aes(lane)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#B0B0B0", linewidth = 0.3) +
  geom_col(data = sm, aes(y = m, fill = lane), width = 0.62, alpha = 0.9) +
  geom_errorbar(data = sm, aes(y = m, ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.22, linewidth = 0.32, color = "grey20") +
  geom_point(aes(y = rel, color = lane), size = 0.9,
             position = position_jitter(width = 0.10, height = 0, seed = 42)) +
  geom_text(data = starpos, aes(x = group, y = y, label = star),
            family = "Arial", size = 2.6, vjust = 0.3) +
  facet_wrap(~ protein, nrow = 1) +
  scale_fill_manual(values = GCOL) + scale_color_manual(values = GCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "Relative band intensity\n(target / β-actin, Sham = 1)",
       title = "m6A machinery Western blot (Fig 6C) — 3 biological replicates",
       subtitle = "bar = mean ± SD of n = 3 independent blots; ANOVA + Tukey HSD vs Sham\nALKBH5 main band quantified; strong bands partially saturated = lower bounds") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(), legend.position = "none",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7.5),
        plot.subtitle = element_text(size = 5.8, color = "grey40"))
ggsave(file.path(OUT, "Fig6C_m6A_WB_densitometry_3reps_20260915.pdf"), p,
       width = 4.7, height = 2.0, device = cairo_pdf)
ggsave(file.path(OUT, "Fig6C_m6A_WB_densitometry_3reps_20260915.png"), p,
       width = 4.7, height = 2.0, dpi = 600, type = "cairo")
cat("DONE\n")
