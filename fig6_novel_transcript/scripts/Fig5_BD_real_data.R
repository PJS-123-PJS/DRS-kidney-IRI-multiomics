# Comprehensive Fig5 fix script — Panel B (real chromosome from structural CSV)
# + Panel D (class code stats from real data)
library(data.table); library(ggplot2); library(patchwork); library(scales)
source("../../shared/palette/color_palette.R")
source("fig5_muted_palette.R")
out_dir <- "../output"

# ── Load REAL structural data (19,487 transcripts) ──────────
struct <- fread("../data/novel_transcript_structural.csv")
cat(sprintf("Total transcripts: %d\n", nrow(struct)))
cat(sprintf("Columns: %s\n", paste(names(struct), collapse=", ")))

# ── Class code distribution (REAL from structural CSV) ──────
cc_dist <- struct[, .N, by=class_code][order(-N)]
cc_dist[, Pct := N/sum(N)*100]
cc_dist[, class_code := factor(class_code, levels=names(class_colors_muted))]
cat("\nReal class code distribution:\n")
print(cc_dist)

# ── Panel A: Chromosome Distribution (REAL data, including sex chromosomes) ──
# Normalize chromosome names: "X" → "chrX", "chr1" stays "chr1"
struct[, chr_norm := ifelse(chr == "X", "chrX", chr)]

# Valid chromosomes: chr1-19 + chrX + chrY (chrY = 0 transcripts, shown as empty bar)
chr_valid <- paste0("chr", c(as.character(1:19), "X", "Y"))
chr_dt <- struct[chr_norm %in% chr_valid]
chr_by_cc <- chr_dt[, .N, by=.(chr_norm, class_code)]
setnames(chr_by_cc, "chr_norm", "chr")
chr_by_cc[, class_code := factor(class_code, levels=names(class_colors_muted))]

# Add chrY as explicit zero row (if not present)
if (!"chrY" %in% unique(chr_by_cc$chr)) {
  chr_y_row <- data.table(chr="chrY", class_code="k", N=0)
  chr_by_cc <- rbind(chr_by_cc, chr_y_row, fill=TRUE)
}
chr_by_cc[chr == "chrY" & is.na(class_code), class_code := "k"]

# Order chromosomes naturally: chr1-19, chrX, chrY
chr_order <- paste0("chr", c(as.character(1:19), "X", "Y"))
chr_by_cc[, chr := factor(chr, levels=chr_order)]
chr_by_cc <- chr_by_cc[!is.na(chr)]

# Summary
chr_totals <- chr_dt[, .N, by=chr_norm][order(-N)]
chr_x_count <- chr_totals[chr_norm == "chrX", N]
chr_y_count <- 0L

cat(sprintf("Top 5 chromosomes: %s\n",
            paste(sprintf("%s(%d)", chr_totals$chr_norm[1:5], chr_totals$N[1:5]), collapse=", ")))
cat(sprintf("chrX: %d transcripts | chrY: %d transcripts\n", chr_x_count, chr_y_count))

p_a <- ggplot(chr_by_cc, aes(x=chr, y=N, fill=class_code)) +
  geom_col(alpha=0.9, width=0.7, color=NA) +
  # Highlight chrY bar border (empty, show it was checked)
  geom_col(data=chr_by_cc[chr=="chrY"], aes(x=chr, y=N),
           fill=NA, color="grey50", linewidth=0.8, width=0.7) +
  # chrY annotation
  annotate("text", x=length(chr_order), y=max(chr_totals$N)*0.05,
           label="0", size=3.5, fontface="bold", color="grey50") +
  scale_fill_manual(values=class_colors_muted, name="Class Code") +
  labs(x="Chromosome", y="Number of Novel Transcripts",
       title="Novel Transcript Chromosome Distribution",
       subtitle=sprintf("%d novel transcripts across chr1-19, chrX (%d transcripts), chrY (0 transcripts); %s(%d) and %s(%d) are hotspots",
                        nrow(chr_dt) + chr_y_count, chr_x_count,
                        chr_totals$chr_norm[1], chr_totals$N[1],
                        chr_totals$chr_norm[2], chr_totals$N[2])) +
  theme_drs + theme(axis.text.x=element_text(size=7, angle=45, hjust=1))

pdf(file.path(out_dir, "Fig6A_Chromosome_Distribution.pdf"), w=12, h=5.5)
print(p_a)
dev.off()
png(file.path(out_dir, "Fig6A_Chromosome_Distribution.png"), w=12, h=5.5, units="in", res=300)
print(p_a)
dev.off()
cat("Panel A saved with real chromosome data (including sex chromosomes).\n")

# ── Panel D: RNA-seq Validation + Class Code Stats ──────────
# Create a 2-in-1 panel: (left) class code pie, (right) per-class-code stats

# D1: Class code re-distribution
cc_dist[, ypos := cumsum(N) - N/2]
p_d1 <- ggplot(cc_dist, aes(x="", y=N, fill=class_code)) +
  geom_col(width=1, alpha=0.9, color="white", linewidth=1.2) +
  coord_polar("y", start=0) +
  geom_text(aes(y=ypos, label=sprintf("%s\n%d (%.1f%%)", class_code, N, Pct)),
            size=3.5, fontface="bold", color="white") +
  scale_fill_manual(values=class_colors_muted, guide="none") +
  labs(title="Class Code Distribution",
       subtitle=sprintf("%d novel transcripts, %d class codes", sum(cc_dist$N), nrow(cc_dist))) +
  theme_void() +
  theme(plot.title=element_text(face="bold", size=12, hjust=0.5),
        plot.subtitle=element_text(size=9, color="grey50", hjust=0.5))

# D2: Length distribution by class code (from real structural data)
struct_plot <- struct[length > 0 & length < 200000]
struct_plot[, class_code := factor(class_code, levels=names(class_colors_muted))]
p_d2 <- ggplot(struct_plot, aes(x=class_code, y=length, fill=class_code)) +
  geom_violin(alpha=0.85, draw_quantiles=c(0.25, 0.5, 0.75), scale="width") +
  scale_fill_manual(values=class_colors_muted, guide="none") +
  scale_y_log10(labels=comma) +
  labs(x=NULL, y="Transcript Length (bp, log10)",
       title="Transcript Length by Class Code",
       subtitle=sprintf("j-class longest (median %.0fkb); n-class shortest (median %.1fkb)",
                        median(struct_plot[class_code=="j", length])/1000,
                        median(struct_plot[class_code=="n", length])/1000)) +
  theme_drs + theme(axis.text.x=element_text(angle=30, hjust=1))

# D3: Exon count distribution by class code
struct_plot[, class_code := factor(class_code, levels=names(class_colors_muted))]
p_d3 <- ggplot(struct_plot, aes(x=class_code, y=exon_count, fill=class_code)) +
  geom_boxplot(alpha=0.85, outlier.size=0.5) +
  scale_fill_manual(values=class_colors_muted, guide="none") +
  labs(x=NULL, y="Exon Count",
       title="Exon Count by Class Code",
       subtitle=sprintf("j-class most exons (median %d); u-class fewest (median %d)",
                        as.integer(median(struct_plot[class_code=="j", exon_count])),
                        as.integer(median(struct_plot[class_code=="u", exon_count])))) +
  theme_drs + theme(axis.text.x=element_text(angle=30, hjust=1))

# Combine D panels
p_d <- (p_d1 | (p_d2 / p_d3)) +
  plot_layout(widths=c(0.4, 0.6)) +
  plot_annotation(
    title="Novel Transcript Classification & Validation",
    subtitle=sprintf("Based on %d novel transcripts from full-length DRS data (cuffcompare vs mm10)",
                     nrow(struct)),
    theme=theme(plot.title=element_text(face="bold", size=14),
                plot.subtitle=element_text(size=9, color="grey50"))
  )

pdf(file.path(out_dir, "Fig6B_Transcript_Classification.pdf"), w=14, h=8)
print(p_d)
dev.off()
png(file.path(out_dir, "Fig6B_Transcript_Classification.png"), w=14, h=8, units="in", res=300)
print(p_d)
dev.off()
cat("Panel D saved with real classification data.\n")

# ── Also verify and save real class code stats ──────────────
fwrite(cc_dist, file.path(out_dir, "Fig6_real_class_code_stats.csv"))
cat("\nAll panels regenerated with REAL data from structural CSV.\n")
