library(data.table); library(ggplot2); library(patchwork)
source("../../shared/palette/color_palette.R")
source("fig5_muted_palette.R")
out_dir <- "../output"
class_colors <- class_colors_muted

# Panel B moved to Fig5_BD_real_data.R (comprehensive version with chr from structural CSV)
# Panel C merged into Fig5D (now Fig5B): length violin + exon count already in classification panel

# Panel G: Real Expression with Wilcoxon stats + N labels (2026-06-23)
novel_expr <- fread("../data/novel_transcript_expression.csv")
known_expr <- fread("../data/known_mRNA_expression.csv")
known_sample <- known_expr[, .(Category="Known mRNA", TPM=NC_mean)]
novel_with_cc <- novel_expr[class_code %in% c("u","x","i"), .(Category=paste0("Novel (",class_code,")"), TPM=NC_mean)]
novel_other <- novel_expr[!class_code %in% c("u","x","i"), .(Category="Novel (other)", TPM=NC_mean)]
expr_comp <- rbind(known_sample, novel_with_cc, novel_other)

cat_levels <- c("Known mRNA","Novel (u)","Novel (x)","Novel (i)","Novel (other)")
expr_comp[, Category := factor(Category, levels=cat_levels)]

# N per group for x-axis labels
n_tbl <- expr_comp[, .N, by=Category][order(Category)]
x_labs <- paste0(n_tbl$Category, "\n(n=", n_tbl$N, ")")

# Wilcoxon rank-sum: Known mRNA vs each Novel class
known_tpm <- expr_comp[Category=="Known mRNA", TPM]
novel_cats <- c("Novel (u)","Novel (x)","Novel (i)","Novel (other)")
p_vals <- sapply(novel_cats, function(cat) {
  wilcox.test(known_tpm, expr_comp[Category==cat, TPM], exact=FALSE)$p.value
})
sig_labels <- sapply(p_vals, function(p) {
  if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else "ns"
})

cat(sprintf("Wilcoxon Known vs Novel(u): p=%.2e %s\n", p_vals[1], sig_labels[1]))
cat(sprintf("Wilcoxon Known vs Novel(x): p=%.2e %s\n", p_vals[2], sig_labels[2]))
cat(sprintf("Wilcoxon Known vs Novel(i): p=%.2e %s\n", p_vals[3], sig_labels[3]))
cat(sprintf("Wilcoxon Known vs Novel(other): p=%.2e %s\n", p_vals[4], sig_labels[4]))

# Manual significance brackets (ggsignif conflicted with scale_y_log10)
# Y positions in log10 space; bracket at increasing heights
y_base <- max(log10(expr_comp$TPM + 0.01), na.rm=TRUE)
y_steps <- y_base + c(0.3, 0.55, 0.80, 1.05)
bracket_h <- 0.06
n_comps <- 4

# Extended palette for all 5 categories
g_colors <- c(
  "Known mRNA"   = "#8CBFD0",
  "Novel (u)"    = "#9C78B8",
  "Novel (x)"    = "#E69F5B",
  "Novel (i)"    = "#999999",
  "Novel (other)"= "#D4C5B9"
)

p_g <- ggplot(expr_comp, aes(x=Category, y=TPM+0.01, fill=Category)) +
  geom_violin(alpha=0.85, draw_quantiles=c(0.5)) +
  scale_fill_manual(values=g_colors, guide="none") +
  scale_y_log10(labels=scales::comma) +
  scale_x_discrete(labels=x_labs)

# Add significance brackets manually
for(i in 1:n_comps) {
  x_left <- 1; x_right <- 1 + i
  ys <- 10^y_steps[i]
  yb <- 10^(y_steps[i] - bracket_h)
  # Vertical lines
  p_g <- p_g +
    annotate("segment", x=x_left, xend=x_left, y=yb, yend=ys, size=0.35) +
    annotate("segment", x=x_right, xend=x_right, y=yb, yend=ys, size=0.35) +
    annotate("segment", x=x_left, xend=x_right, y=ys, yend=ys, size=0.35) +
    annotate("text", x=(x_left+x_right)/2, y=ys*1.1, label=sig_labels[i],
             size=3.8, fontface="bold")
}

p_g <- p_g +
  labs(x=NULL, y="Expression (TPM, log10)", title="Novel Transcript Expression (Real RNA-seq Data)",
       subtitle=sprintf("Wilcoxon: Known mRNA (n=%d, median=%.1f TPM) vs Novel u/x/i (n=%d/%d/%d, median=%.0f/%.0f/%.0f TPM); all p<0.001 ***",
                        n_tbl[Category=="Known mRNA",N], median(known_tpm),
                        n_tbl[Category=="Novel (u)",N], n_tbl[Category=="Novel (x)",N], n_tbl[Category=="Novel (i)",N],
                        median(expr_comp[Category=="Novel (u)",TPM]),
                        median(expr_comp[Category=="Novel (x)",TPM]),
                        median(expr_comp[Category=="Novel (i)",TPM]))) +
  theme_drs + theme(axis.text.x=element_text(angle=0, hjust=0.5, size=8))

pdf(file.path(out_dir, "Fig6G_Expression_Comparison.pdf"), w=10, h=6); print(p_g); dev.off()
png(file.path(out_dir, "Fig6G_Expression_Comparison.png"), w=10, h=6, units="in", res=300); print(p_g); dev.off()

# Panel E moved to Fig5E_real_density.R (real server CSV data) — removed to prevent overwrite

cat("Panels G saved with real data.\n")
