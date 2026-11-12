# Fig8 Panel H — PolyA PAT real-data panel, ALL-GROUP + gene swap (2026-07-31)
# ============================================================================
# Original Fig8H was a TEXT TEMPLATE (Park7+Gpx1) reading the AGING-PROJECT
# file Module4_polyA_aging_per_gene.tsv — two problems (2026-07-23 reevaluation):
#   (1) Park7 is RETIRED (HC50 verified 0 sites — fabrication-chain gene);
#   (2) cross-project data routing violates the project boundary rule.
# Rebuilt from canonical E:result/06.polyA/length/01.all per-read files:
#   gene swap  Park7 -> Fxyd2 + Kap (both canonical NC_vs_Day1 top shortening
#   hits: Kap dmedian -18..-28 nt FDR=0, Fxyd2 -7..-22 nt FDR<=1.6e-26);
#   Gpx1 retained (real, verified: -13..-20 nt FDR=0).
# Per-sample gene-level median polyA (PAT convention) across ALL of the gene's
# annotated transcripts (19 contigs from isoforms.TPM.matrix.add_info), then
# group mean±SEM (n=3). Old template archived to archive_superseded_20260731/.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

READ_DIR <- "/tmp/pat_reads"
CONTIGS  <- "/tmp/pat_contigs.tsv"
OUT_DIR  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
ARCHIVE  <- file.path(OUT_DIR, "archive_superseded_20260731")
dir.create(ARCHIVE, showWarnings = FALSE, recursive = TRUE)

for (ext in c("pdf","png")) {
  src <- file.path(OUT_DIR, paste0("Fig8H_PAT_Assay_Template.", ext))
  if (file.exists(src)) file.rename(src, file.path(ARCHIVE, basename(src)))
}

GROUPS   <- c("NC", "IS", "Day1", "Day3", "Day7")
G_COLORS <- c("NC" = "#4A4A4A", "IS" = "#8A8A8A", "Day1" = "#D6404C",
              "Day3" = "#E8753A", "Day7" = "#2E7D6F")
SAMPLE_GROUP <- c(NCs1="NC",NCs2="NC",NCs3="NC", Is1="IS",Is2="IS",Is3="IS",
                  Day1s1="Day1",Day1s2="Day1",Day1s3="Day1",
                  Day3s1="Day3",Day3s2="Day3",Day3s3="Day3",
                  Day7s1="Day7",Day7s2="Day7",Day7s3="Day7")

c2g <- fread(CONTIGS, header = FALSE, col.names = c("Gene","contig"))

# ── per-sample gene-level median polyA ──
per_sample <- list()
for (s in names(SAMPLE_GROUP)) {
  f <- file.path(READ_DIR, paste0(s, ".tsv"))
  d <- fread(f, header = FALSE, col.names = c("readname","contig","polya_length"))
  d <- merge(d, c2g, by = "contig")
  ps <- d[, .(median_polyA = median(polya_length), n_reads = .N), by = Gene]
  ps[, Sample := s][, Group := SAMPLE_GROUP[[s]]]
  per_sample[[s]] <- ps
  cat(sprintf("  %s: %s reads\n", s, format(sum(ps$n_reads), big.mark=",")))
}
ps <- rbindlist(per_sample)
ps[, Group := factor(Group, levels = GROUPS)]
ps[, Gene  := factor(Gene, levels = c("Fxyd2","Kap","Gpx1"))]

cat("=== per-sample gene median polyA (nt) ===\n")
print(dcast(ps, Gene + Group ~ Sample, value.var = "median_polyA"))

stat <- ps[, .(m = mean(median_polyA), se = sd(median_polyA)/sqrt(.N),
               reads = sum(n_reads)), by = .(Gene, Group)]

# Day1 delta annotation per gene
dlt <- stat[Group == "Day1", .(Gene, d1 = m)][stat[Group == "NC", .(Gene, nc = m)], on = "Gene"]
dlt[, lab := sprintf("ΔDay1 = %+.1f nt", d1 - nc)]
cat("=== Day1 delta vs NC ===\n"); print(dlt)

p_h <- ggplot() +
  geom_col(data = stat, aes(x = Group, y = m, fill = Group),
           width = 0.62, alpha = 0.85) +
  geom_errorbar(data = stat, aes(x = Group, ymin = m - se, ymax = m + se),
                width = 0.22, linewidth = 0.3, color = "grey30") +
  geom_point(data = ps, aes(x = Group, y = median_polyA),
             size = 1.1, alpha = 0.75, color = "grey20",
             position = position_jitter(width = 0.14, seed = 7)) +
  geom_text(data = dlt, aes(x = 3, y = -Inf, label = lab),
            vjust = -0.8, size = 2.1, fontface = "bold", family = "Arial",
            color = "#D6404C") +
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = G_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.05))) +
  labs(x = NULL, y = "Median polyA tail length (nt)",
       title = expression(bold("h")~~"PolyA Tail Length (PAT): Validated Targets — Full Time Course"),
       subtitle = "DRS per-read polyA, gene-level median per sample (n=3/group, bar=mean±SEM); Park7 replaced by Fxyd2/Kap (canonical NC_vs_Day1 top hits)",
       caption = "Fxyd2/Kap/Gpx1 all FDR<=1.6e-26 (Kap/Gpx1 FDR=0) at Day1. Transient shortening at Day1, rebound toward Day3-Day7.") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.3, color = "grey30"),
        plot.caption = element_text(size = 5.8, color = "grey50", hjust = 0),
        strip.text = element_text(size = 7.5, face = "bold.italic"),
        strip.background = element_blank(),
        panel.spacing = unit(2, "mm"))

ggsave(file.path(OUT_DIR, "Fig8H_PAT_AllGroups.pdf"), p_h,
       width = 183/25.4, height = 80/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig8H_PAT_AllGroups.png"), p_h,
       width = 183/25.4, height = 80/25.4, dpi = 600)

fwrite(ps,   file.path(OUT_DIR, "Table_Fig8H_per_sample_gene_median_polyA.csv"))
fwrite(stat, file.path(OUT_DIR, "Table_Fig8H_group_stats.csv"))
cat("=== Fig8H all-group real-data panel saved ===\n")
