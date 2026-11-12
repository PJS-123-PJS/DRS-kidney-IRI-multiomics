# Fig8B wet-lab qPCR visualization — 10 modification enzymes, 4 groups (NC/Day1/Day3/Day7)
# Input : Prism_Grouped格式_qPCR数据.xlsx (sheet 长格式_原始RQ; RQ = 2^-ΔΔCt, n=3)
# Output: Fig8B_qPCR_Wetlab_10Enzymes.{pdf,png} + stats/summary CSV
# Style : per-gene free-y mini-bar array (mirrors meRIP 3-layer final v4), DRS canonical
#         group colors (NC#4A4A4A / D1#D6404C / D3#E8753A / D7#2E7D6F), Arial 7pt.
#         Gene order = Day1 mean log2(RQ) descending — mirrors finalized Fig8B heatmap
#         ("ordered by Day1 response"). Stats: Welch t-test on log2(RQ) vs NC.
# do NOT use ggbreak/stat_summary(mean_sdl) — precompute mean±SD (Hmisc-free).

suppressPackageStartupMessages({
  library(data.table); library(readxl); library(ggplot2)
})

XLSX <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/Prism_Grouped格式_qPCR数据.xlsx"
OUT  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
font_family <- "Arial"
GROUPS  <- c("NC","DAY1","DAY3","DAY7")
G_LABEL <- c(NC="NC", DAY1="Day1", DAY3="Day3", DAY7="Day7")
G_COL   <- c(NC="#4A4A4A", DAY1="#D6404C", DAY3="#E8753A", DAY7="#2E7D6F")

dt <- as.data.table(read_excel(XLSX, sheet = "长格式_原始RQ"))
stopifnot(all(c("Gene","Replicate","Group","RQ") %in% names(dt)))
dt[, Group := factor(Group, levels = GROUPS)]
dt[, log2RQ := log2(RQ)]
cat("genes:", length(unique(dt$Gene)), "| rows:", nrow(dt), "\n")

# gene order: Day1 mean log2RQ descending (mirror Fig8B heatmap)
ord <- dt[Group == "DAY1", .(m = mean(log2RQ)), by = Gene][order(-m)]
dt[, Gene := factor(Gene, levels = ord$Gene)]
cat("order:", paste(ord$Gene, collapse = " > "), "\n")

# summary + Welch t on log2RQ vs NC
sm <- dt[, .(n = .N, m = mean(RQ), s = sd(RQ), mlog2 = mean(log2RQ)), by = .(Gene, Group)]
st <- rbindlist(lapply(split(dt, dt$Gene), function(d){
  nc <- d[Group == "NC", log2RQ]
  rbindlist(lapply(GROUPS[-1], function(g){
    tt <- t.test(d[Group == g, log2RQ], nc)   # Welch (default)
    data.table(Gene = as.character(d$Gene[1]), Group = g,
               t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
  }))
}))
st[, star := fifelse(p < 0.0001, "****", fifelse(p < 0.001, "***",
              fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", ""))))]
st[, Gene := factor(Gene, levels = levels(dt$Gene))]
st[, Group := factor(Group, levels = GROUPS)]
fwrite(st, file.path(OUT, "Table_Fig8B_qPCR_wetlab_stats_20260914.csv"))
sm_out <- merge(sm, st[, .(Gene, Group, p, star)], by = c("Gene","Group"), all.x = TRUE)
setorder(sm_out, Gene, Group)
fwrite(sm_out, file.path(OUT, "Table_Fig8B_qPCR_wetlab_summary_20260914.csv"))

# star positions: above max(mean+SD, max replicate)
rp <- dt[, .(rmax = max(RQ)), by = .(Gene, Group)]
starpos <- merge(st[star != ""], sm, by = c("Gene","Group"))
starpos <- merge(starpos, rp, by = c("Gene","Group"))
starpos[, y := pmax(m + s, rmax)]
starpos[, Gene := factor(Gene, levels = levels(dt$Gene))]

p <- ggplot(dt, aes(Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#B0B0B0", linewidth = 0.3) +
  geom_col(data = sm, aes(y = m, fill = Group), width = 0.62, alpha = 0.88, color = NA) +
  geom_errorbar(data = sm, aes(ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.22, linewidth = 0.32, color = "grey20") +
  geom_point(aes(y = RQ, color = Group), size = 0.9,
             position = position_jitter(width = 0.10, height = 0, seed = 42)) +
  geom_text(data = starpos, aes(y = y, label = star),
            family = font_family, size = 2.5, vjust = 0.35) +
  facet_wrap(~ Gene, nrow = 2, scales = "free_y") +
  scale_x_discrete(labels = G_LABEL) +
  scale_fill_manual(values = G_COL, labels = G_LABEL) +
  scale_color_manual(values = G_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  labs(x = NULL, fill = NULL,
       y = expression("Relative mRNA expression (2"^{-Delta*Delta*Ct}*")"),
       title = "qRT-PCR validation: modification-enzyme expression across IRI time course",
       subtitle = "n = 3 biological replicates; bar = mean ± SD; Welch t-test on log2(RQ) vs NC; genes ordered by Day1 response") +
  theme_bw(base_size = 7, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold.italic", size = 7),
    axis.text.x = element_text(size = 6),
    legend.position = "top",
    legend.key.size = unit(2.6, "mm"),
    legend.text = element_text(size = 6.5),
    plot.title = element_text(face = "bold", size = 8.5),
    plot.subtitle = element_text(size = 6.5, color = "grey30"),
    panel.border = element_rect(linewidth = 0.4)
  )

ggsave(file.path(OUT, "Fig8B_qPCR_Wetlab_10Enzymes.pdf"), p,
       width = 180/25.4, height = 78/25.4, device = cairo_pdf)
ggsave(file.path(OUT, "Fig8B_qPCR_Wetlab_10Enzymes.png"), p,
       width = 180/25.4, height = 78/25.4, dpi = 600, type = "cairo")
cat("DONE\n")
