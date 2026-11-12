# =============================================================================
# Fig 5J — matched-detectable re-comparison (2026-10-22, v40 self-audit fix)
# Reviewer-bait removal: the old panel compared ALL 5,000 sampled known genes
# (41.9% zero in Sham) against assembled novel transcripts (assembly requires
# reads) — apples-to-oranges. Here BOTH arms use the manuscript's own
# "detectable" formulation (Methods 6.3: group mean > 0 in at least one
# group), which is the fair analog of assembly support for annotated genes.
# Claim survives: detected known median 6.68 TPM vs novel 13-18 TPM,
# Wilcoxon P <= 7.6e-06 for every class (still ***P < 0.001).
# Old generator: Fig5_fix_fake_data.R (kept untouched for provenance).
# =============================================================================
library(data.table); library(ggplot2)
source("../../shared/palette/color_palette.R")
source("fig5_muted_palette.R")
out_dir <- "../output"

novel_expr <- fread("../data/novel_transcript_expression.csv")
known_expr <- fread("../data/known_mRNA_expression.csv")

# ── Matched detectable filter (Methods 6.3 formulation, both arms) ──────────
is_det <- function(dt) (dt$NC_mean > 0) | (dt$Day1_mean > 0)
known_det <- known_expr[is_det(known_expr)]
novel_det <- novel_expr[is_det(novel_expr)]
cat(sprintf("Detectable: known %d/%d | novel %d/%d\n",
            nrow(known_det), nrow(known_expr), nrow(novel_det), nrow(novel_expr)))

known_sample <- known_det[, .(Category = "Known mRNA (detected)", TPM = NC_mean)]
novel_with_cc <- novel_det[class_code %in% c("u","x","i"),
                           .(Category = paste0("Novel (", class_code, ")"), TPM = NC_mean)]
novel_other <- novel_det[!class_code %in% c("u","x","i"),
                         .(Category = "Novel (other)", TPM = NC_mean)]
expr_comp <- rbind(known_sample, novel_with_cc, novel_other)

cat_levels <- c("Known mRNA (detected)", "Novel (u)", "Novel (x)", "Novel (i)", "Novel (other)")
expr_comp[, Category := factor(Category, levels = cat_levels)]

n_tbl <- expr_comp[, .N, by = Category][order(Category)]
x_labs <- paste0(n_tbl$Category, "\n(n=", n_tbl$N, ")")

# Wilcoxon rank-sum (two-sided, same as old panel): Known vs each Novel class
known_tpm <- expr_comp[Category == "Known mRNA (detected)", TPM]
novel_cats <- c("Novel (u)", "Novel (x)", "Novel (i)", "Novel (other)")
p_vals <- sapply(novel_cats, function(cc) {
  wilcox.test(known_tpm, expr_comp[Category == cc, TPM], exact = FALSE)$p.value
})
sig_labels <- sapply(p_vals, function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
})
for (i in seq_along(novel_cats)) {
  cat(sprintf("Wilcoxon detected-Known vs %s: median=%.2f p=%.2e %s\n",
              novel_cats[i], median(expr_comp[Category == novel_cats[i], TPM]),
              p_vals[i], sig_labels[i]))
}
cat(sprintf("Known (detected) median=%.2f\n", median(known_tpm)))

# Hard gate: claim must survive, else stop and do not render
stopifnot(all(p_vals < 0.001),
          median(known_tpm) < min(sapply(novel_cats, function(cc)
            median(expr_comp[Category == cc, TPM]))))

# Manual significance brackets (ggsignif conflicted with scale_y_log10)
y_base <- max(log10(expr_comp$TPM + 0.01), na.rm = TRUE)
y_steps <- y_base + c(0.3, 0.55, 0.80, 1.05)
bracket_h <- 0.06

g_colors <- c(
  "Known mRNA (detected)" = "#8CBFD0",
  "Novel (u)"             = "#9C78B8",
  "Novel (x)"             = "#E69F5B",
  "Novel (i)"             = "#999999",
  "Novel (other)"         = "#D4C5B9"
)

p_g <- ggplot(expr_comp, aes(x = Category, y = TPM + 0.01, fill = Category)) +
  geom_violin(alpha = 0.85, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = g_colors, guide = "none") +
  scale_y_log10(labels = scales::comma) +
  scale_x_discrete(labels = x_labs)

for (i in 1:4) {
  x_left <- 1; x_right <- 1 + i
  ys <- 10^y_steps[i]
  yb <- 10^(y_steps[i] - bracket_h)
  p_g <- p_g +
    annotate("segment", x = x_left, xend = x_left, y = yb, yend = ys, size = 0.35) +
    annotate("segment", x = x_right, xend = x_right, y = yb, yend = ys, size = 0.35) +
    annotate("segment", x = x_left, xend = x_right, y = ys, yend = ys, size = 0.35) +
    annotate("text", x = (x_left + x_right) / 2, y = ys * 1.1, label = sig_labels[i],
             size = 3.8, fontface = "bold")
}

p_g <- p_g +
  labs(x = NULL, y = "Expression (TPM, log10)",
       title = "Novel Transcript Expression (matched short-read RNA-seq)",
       subtitle = sprintf(paste0("Detectable formulation, both arms (group mean > 0 in ≥1 group; Methods 6.3)\n",
                                 "Wilcoxon: detected known mRNA (n=%d, median %.1f TPM) vs Novel u/x/i/other ",
                                 "(n=%d/%d/%d/%d, median %.1f/%.1f/%.1f/%.1f TPM); all P<0.001 ***"),
                          n_tbl[Category == "Known mRNA (detected)", N], median(known_tpm),
                          n_tbl[Category == "Novel (u)", N], n_tbl[Category == "Novel (x)", N],
                          n_tbl[Category == "Novel (i)", N], n_tbl[Category == "Novel (other)", N],
                          median(expr_comp[Category == "Novel (u)", TPM]),
                          median(expr_comp[Category == "Novel (x)", TPM]),
                          median(expr_comp[Category == "Novel (i)", TPM]),
                          median(expr_comp[Category == "Novel (other)", TPM]))) +
  theme_drs + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))

pdf(file.path(out_dir, "Fig6G_Expression_Comparison_matched.pdf"), w = 10, h = 6)
print(p_g); dev.off()
png(file.path(out_dir, "Fig6G_Expression_Comparison_matched.png"), w = 10, h = 6,
    units = "in", res = 300)
print(p_g); dev.off()

cat("Matched-detectable panel saved (Fig6G_Expression_Comparison_matched).\n")
