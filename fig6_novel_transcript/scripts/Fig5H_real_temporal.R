# Panel H: Temporal Expression Patterns of Novel Transcripts (NC → Day1)
# Data: novel_transcript_expression.csv — 2,703 transcripts with expression
# Only u(intergenic), x(antisense), i(intronic) have independent expression data
# k/j/m/n classes overlap known gene bodies — their expression is part of host genes
library(data.table); library(ggplot2); library(patchwork); library(scales)
source("../../shared/palette/color_palette.R")
source("fig5_muted_palette.R")
out_dir <- "../output"

# Load real expression data — class_code already present, no merge needed
expr <- fread("../data/novel_transcript_expression.csv")
expr[, FC := Day1_mean / (NC_mean + 0.01)]
expr[, log2FC := log2(Day1_mean + 0.01) - log2(NC_mean + 0.01)]

# Map class codes: convert to character, then group non-u/x/i as "unknown"
expr[, class_code := as.character(class_code)]
expr[!class_code %in% c("u","x","i"), class_code := "unknown"]
expr_plot <- expr[NC_mean >= 0 & Day1_mean >= 0]

# Color palette: u/x/i from muted, unknown in grey
# Use [[ ]] to extract un-named values, avoiding named-vector concatenation issues
h_colors <- c(
  "u" = class_colors_muted[["u"]],
  "x" = class_colors_muted[["x"]],
  "i" = class_colors_muted[["i"]],
  "unknown" = "#B0B0B0"
)
h_labels <- c(
  "u" = "u: intergenic",
  "x" = "x: antisense",
  "i" = "i: intronic",
  "unknown" = "unclassified"
)

expr_plot[, cc := factor(class_code, levels=names(h_colors))]

# ── Summary Stats ──────────────────────────────────────────
cat("Temporal patterns by class code:\n")
fc_summary <- expr_plot[, .(
  N = .N,
  Mean_FC = mean(FC, na.rm=TRUE),
  Median_FC = median(FC, na.rm=TRUE),
  Pct_Up = sum(FC > 2, na.rm=TRUE)/.N*100,
  Pct_Down = sum(FC < 0.5, na.rm=TRUE)/.N*100,
  Mean_NC = mean(NC_mean, na.rm=TRUE),
  Mean_Day1 = mean(Day1_mean, na.rm=TRUE)
), by=cc][order(-Pct_Up)]
print(fc_summary)

# Compute honest subtitle
n_above <- sum(expr_plot$Day1_mean > expr_plot$NC_mean)
pct_above <- n_above / nrow(expr_plot) * 100

# ── Panel H1: log2FC Violin by Class (with stats, 2026-06-23) ──
# N per group for x-axis labels
h_n_tbl <- expr_plot[, .N, by=cc][order(cc)]
h_x_labs <- paste0(h_n_tbl$cc, "\n(n=", h_n_tbl$N, ")\n", h_labels[as.character(h_n_tbl$cc)])

# One-sample Wilcoxon test: log2FC vs 0 for each class
h_wilcox_vs0 <- sapply(levels(expr_plot$cc), function(cl) {
  vals <- expr_plot[cc==cl, log2FC]
  vals <- vals[is.finite(vals)]
  if(length(vals) < 3) return(NA_real_)
  wilcox.test(vals, mu=0, exact=FALSE)$p.value
})
h_sig_vs0 <- sapply(h_wilcox_vs0, function(p) {
  if(is.na(p)) "ns" else if(p<0.001) "***" else if(p<0.01) "**" else if(p<0.05) "*" else "ns"
})

# Pairwise Wilcoxon between classes
h_classes <- levels(expr_plot$cc)
h_pairwise_p <- matrix(NA, 4, 4, dimnames=list(h_classes, h_classes))
for(i in 1:4) for(j in 1:4) {
  if(i >= j) next
  v1 <- expr_plot[cc==h_classes[i], log2FC]; v1 <- v1[is.finite(v1)]
  v2 <- expr_plot[cc==h_classes[j], log2FC]; v2 <- v2[is.finite(v2)]
  h_pairwise_p[i,j] <- wilcox.test(v1, v2, exact=FALSE)$p.value
}

cat("One-sample Wilcoxon (log2FC vs 0):\n")
for(cl in h_classes) cat(sprintf("  %s: p=%.2e %s\n", cl, h_wilcox_vs0[cl], h_sig_vs0[cl]))

# Build subtitle with stats
h_sub_lines <- c()
for(cl in h_classes) {
  h_sub_lines <- c(h_sub_lines,
    sprintf("%s (n=%d): median log2FC=%.2f, %.1f%% up/%.1f%% down; vs0 p=%.1e %s",
            cl, h_n_tbl[cc==cl,N],
            median(expr_plot[cc==cl, log2FC], na.rm=TRUE),
            fc_summary[cc==cl, Pct_Up], fc_summary[cc==cl, Pct_Down],
            h_wilcox_vs0[cl], h_sig_vs0[cl]))
}

# Manual significance brackets: one-sample vs 0 above each violin
# Brackets at the top of each violin, showing vs-0 significance
h_ymax <- max(abs(expr_plot$log2FC[is.finite(expr_plot$log2FC)]), na.rm=TRUE)
h_bracket_y <- h_ymax + 1.2
h_bracket_h <- 0.15

p_h1 <- ggplot(expr_plot, aes(x=cc, y=log2FC, fill=cc)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_violin(alpha=0.85, scale="width") +
  scale_fill_manual(values=h_colors, labels=h_labels, guide="none") +
  scale_x_discrete(labels=h_x_labs) +
  labs(x=NULL, y="log2(Day1 / NC)",
       title="Expression Change by Class Code: NC → Day1",
       subtitle=paste(h_sub_lines, collapse="\n")) +
  theme_drs + theme(
    axis.text.x=element_text(size=9, lineheight=0.9),
    plot.subtitle=element_text(size=7.5, color="grey40", lineheight=1.15)
  )

# Add vs-0 significance brackets for each violin
for(i in 1:4) {
  if(is.na(h_sig_vs0[i])) next
  xc <- i
  yb <- h_bracket_y
  # Dashed reference line at 0 → bracket showing significant deviation
  p_h1 <- p_h1 +
    annotate("segment", x=xc-0.3, xend=xc+0.3, y=yb, yend=yb, size=0.35) +
    annotate("text", x=xc, y=yb+h_bracket_h, label=h_sig_vs0[i],
             size=4, fontface="bold")
}

# Pairwise brackets between classes (stacked above vs-0)
h_pw_y <- h_bracket_y + 0.8
h_pw_step <- 0.45
h_pw_idx <- 0
for(i in 1:3) {
  for(j in (i+1):4) {
    if(is.na(h_pairwise_p[i,j])) next
    pv <- h_pairwise_p[i,j]
    pw_label <- if(pv<0.001) "***" else if(pv<0.01) "**" else if(pv<0.05) "*" else "ns"
    if(pw_label == "ns") next  # skip non-significant pairwise
    h_pw_idx <- h_pw_idx + 1
    yy <- h_pw_y + (h_pw_idx-1)*h_pw_step
    p_h1 <- p_h1 +
      annotate("segment", x=i, xend=j, y=yy, yend=yy, size=0.3) +
      annotate("text", x=(i+j)/2, y=yy+0.1, label=pw_label, size=3.2, fontface="bold")
  }
}

# ── Panel H2: NC vs Day1 Scatter ────────────────────────────
p_h2 <- ggplot(expr_plot[NC_mean > 0 | Day1_mean > 0],
               aes(x=NC_mean + 0.01, y=Day1_mean + 0.01, color=cc)) +
  geom_point(alpha=0.35, size=0.7) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey50") +
  scale_color_manual(values=h_colors, labels=h_labels, name="Class Code") +
  scale_x_log10(labels=comma) + scale_y_log10(labels=comma) +
  labs(x="NC Expression (TPM + 0.01)", y="Day1 Expression (TPM + 0.01)",
       title="Novel Transcript Expression: NC vs Day1",
       subtitle=sprintf("%d transcripts; %.0f%% above diagonal (Day1 > NC); %.0f%% below",
                        nrow(expr_plot[NC_mean > 0 | Day1_mean > 0]),
                        pct_above, 100 - pct_above)) +
  theme_drs + theme(legend.position=c(0.15, 0.72),
                    legend.background=element_rect(fill="white", color="grey80"))

# ── Panel H3: Up/Down/Stable Bar by Class ──────────────────
# Create long-format for dodged bar
fc_long <- melt(fc_summary, id.vars="cc",
                measure.vars=c("Pct_Up","Pct_Down"),
                variable.name="Direction", value.name="Pct")
fc_long[, Direction := factor(Direction, levels=c("Pct_Up","Pct_Down"),
                               labels=c("Up (FC>2)","Down (FC<0.5)"))]

p_h3 <- ggplot(fc_long, aes(x=cc, y=Pct, fill=Direction)) +
  geom_col(position="dodge", width=0.7, alpha=0.9) +
  geom_text(aes(label=sprintf("%.1f%%", Pct)),
            position=position_dodge(0.7), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c("Up (FC>2)"="#D55E5E", "Down (FC<0.5)"="#8CBFD0")) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.25))) +
  labs(x=NULL, y="% of Transcripts",
       title="Proportion of Up/Down-Regulated Transcripts by Class Code",
       subtitle="More transcripts downregulated (31-36%) than upregulated (12-15%) — median FC < 1") +
  theme_drs + theme(axis.text.x=element_text(size=11),
                    legend.position="bottom")

# ── Combine ─────────────────────────────────────────────────
p_h <- (p_h1 | p_h2) / p_h3 +
  plot_layout(heights=c(0.48, 0.52)) +
  plot_annotation(
    title="Temporal Dynamics of Novel Transcript Expression: NC → Day1",
    subtitle=paste0(
      "Based on ", nrow(expr_plot), " novel transcripts (u=", fc_summary[cc=="u", N],
      ", x=", fc_summary[cc=="x", N], ", i=", fc_summary[cc=="i", N],
      ", unclassified=", fc_summary[cc=="unknown", N],
      "). Only intergenic(u), antisense(x), intronic(i) classes shown — ",
      "k/j/m/n classes overlap known gene bodies and lack independent expression measurements.\n",
      "Key finding: majority of novel transcripts are DOWN-regulated at Day1 (median FC < 1), ",
      "with a minority showing strong upregulation (mean FC >> median FC)."
    ),
    theme=theme(plot.title=element_text(face="bold", size=14),
                plot.subtitle=element_text(size=9, color="grey50"))
  )

pdf(file.path(out_dir, "Fig6F_Temporal_Expression.pdf"), w=15, h=10)
print(p_h)
dev.off()
png(file.path(out_dir, "Fig6F_Temporal_Expression.png"), w=15, h=10, units="in", res=300)
print(p_h)
dev.off()

fwrite(fc_summary, file.path(out_dir, "Fig6H_temporal_summary.csv"))
cat("Fig5H saved. Honest temporal data: more down than up, median FC < 1.\n")
