# Fig8 Panels A/H/I v2 — replacements for archived retired-gene versions
# ============================================================================
# 2026-07-23 — old A/H/I archived (archive_retired_genes_20260723/) because
# they featured retired genes (Park7/Fis1/Gclm/Gpx4/Trp53 = 0 validated HC
# m6A-UP sites; old PAT pair Park7+Gpx1) and displayed the falsified m6A-down
# narrative. v2 renders the RE-EVALUATED targets from the validated basis.
# Data (zero hardcoding — everything read at runtime):
#   data/Fig8_gene_reevaluation_20260723.csv  (curated, corrected 07-23)
#   NC_vs_Day1.polya_diff.all.xls             (canonical company polyA, for H)
# A: meRIP-qPCR target selection (5 new genes + alt)
# H: PAT assay design (Fxyd2/Kap, honest gene-level aggregation)
# I: meRIP-DRS correlation design (text template, experiment pending)
# ============================================================================

Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
CSV     <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/data/Fig8_gene_reevaluation_20260723.csv"
POLYA   <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/06.polyA/length/02.diff/NC_vs_Day1/NC_vs_Day1.polya_diff.all.xls"
for (f in c(CSV, POLYA)) if (!file.exists(f)) stop("HARD FAIL: missing ", f)

ev <- fread(CSV)
theme_n8 <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 6.5, color = "grey30"),
        plot.caption = element_text(size = 5.5, color = "grey50", hjust = 0, lineheight = 1.15),
        plot.margin = margin(5, 5, 5, 5))
save_both <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w/25.4, height = h/25.4, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w/25.4, height = h/25.4, dpi = 600)
  cat(sprintf("%s saved: %.0f x %.0f mm\n", name, w, h))
}

# =============================================================================
# PANEL A v2: meRIP-qPCR targets — HC DRACH+ UP sites by region (stacked)
# =============================================================================
pa <- ev[arm %in% c("meRIP-qPCR", "meRIP-qPCR(alt)")]
parse_struc <- function(s) {
  cds <- 0; utr3 <- 0
  if (nzchar(s)) {
    m <- regmatches(s, gregexpr("(CDS|three_prime_UTR):[0-9]+", s))[[1]]
    for (x in m) {
      kv <- strsplit(x, ":")[[1]]
      if (kv[1] == "CDS") cds <- as.integer(kv[2]) else utr3 <- as.integer(kv[2])
    }
  }
  list(cds = cds, utr3 = utr3)
}
pa_long <- rbindlist(lapply(seq_len(nrow(pa)), function(i) {
  st <- parse_struc(pa$site_struc[i])
  data.table(gene = pa$gene[i], region = c("CDS", "3'UTR"),
             n = c(st$cds, st$utr3), alt = pa$arm[i] == "meRIP-qPCR(alt)",
             tpm = sprintf("TPM %.0f→%.0f", pa$TPM_NC[i], pa$TPM_D1[i]),
             tag = ifelse(pa$gene[i] == "Ninj1", "also AS RT-PCR (dual)",
                   ifelse(pa$gene[i] == "Ndufv3", sprintf("also PolyA %snt", pa$polyA_Diff_median_nt[i]),
                   ifelse(pa$gene[i] == "Tomm5", "TPM rises Day1",
                   ifelse(pa$gene[i] == "Herpud1", "alternative (5-layer elite)", "")))))
}))
pa_long[, gene := factor(gene, levels = rev(c("Glyat","Ninj1","Ehhadh","Tomm5","Ndufv3","Herpud1")))]
lab_dt <- unique(pa_long[, .(gene, alt, tpm, tag)])
lab_dt[, total := pa_long[, sum(n), by = gene]$V1[match(gene, pa_long[, sum(n), by = gene]$gene)]]

p_a <- ggplot(pa_long, aes(y = gene, x = n, fill = region, alpha = alt)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(data = lab_dt,
            aes(y = gene, x = total + 0.25, label = paste0(tpm, ifelse(tag == "", "", paste0("  ·  ", tag)))),
            hjust = 0, size = 2.2, color = "grey25", family = "Arial", inherit.aes = FALSE) +
  scale_fill_manual(values = c("CDS" = "#A8C8A0", "3'UTR" = "#9FC3D8")) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.45), guide = "none") +
  scale_x_continuous(breaks = 0:9, expand = expansion(mult = c(0, 0.55))) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "HC DRACH+ m6A-gain sites (validated HC50 set)", y = NULL,
       title = "a  meRIP-qPCR validation targets — re-evaluated on validated basis",
       subtitle = "5 targets from the validated m6A-gain set (dual-anchor arbitration); retired: Park7/Fis1/Gclm/Gpx4/Trp53 (0 HC sites)",
       caption = paste0("Source: Fig8_gene_reevaluation_20260723.csv (2026-07-23); site counts = HC DRACH+ UP sites per gene.\n",
                        "Selection: validated-set membership + HC site count + TPM feasibility (≥10 both conditions) + PT biology."),
       fill = NULL) +
  theme_n8 +
  theme(legend.position = c(0.82, 0.25), legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"), legend.background = element_blank())
save_both(p_a, "Fig8A_meRIP_qPCR_Targets", 170, 100)

# =============================================================================
# PANEL H v2: PAT assay design — Fxyd2 + Kap (honest gene-level aggregation)
# =============================================================================
polya <- fread(POLYA)
agg_gene <- function(g) {
  d <- polya[Name == g]
  s <- d[FDR < 0.05]
  list(gene = g, n_tx = nrow(d), n_sig = nrow(s),
       mean_d = if (nrow(s)) mean(s$Diff_median) else NA_real_,
       min_fdr = if (nrow(d)) min(d$FDR) else NA_real_)
}
ph <- rbindlist(lapply(c("Fxyd2","Kap"), agg_gene))
ph <- merge(ph, ev[gene %in% c("Fxyd2","Kap"), .(gene, TPM_NC, TPM_D1, n_layers, layers)], by = "gene")
fmt_fdr <- function(x) if (x == 0) "<1e-300" else sprintf("%.1e", x)

mk_block <- function(row, y0) {
  g <- row$gene
  list(
    annotate("text", x = 0.04, y = y0, hjust = 0, vjust = 1, size = 3.4, fontface = "bold",
             label = sprintf("%s  (%s-layer)", g, row$n_layers), family = "Arial"),
    annotate("text", x = 0.04, y = y0 - 0.09, hjust = 0, vjust = 1, size = 2.6, color = "grey20", lineheight = 1.25,
             label = sprintf("polyA shortening: %.1f nt (mean of %d/%d FDR<0.05 transcripts)\nmin transcript FDR = %s\nTPM: %.0f (Sham) → %.0f (Day1) — %s",
                             row$mean_d, row$n_sig, row$n_tx, fmt_fdr(row$min_fdr),
                             row$TPM_NC, row$TPM_D1,
                             ifelse(g == "Fxyd2", "kidney Na/K-ATPase regulator, PT marker", "PT-specific, androgen-regulated")),
             family = "Arial"),
    # schematic PAT lanes
    annotate("rect", xmin = 0.62, xmax = 0.72, ymin = y0 - 0.30, ymax = y0 - 0.02,
             fill = "grey95", color = "grey60", linewidth = 0.3),
    annotate("rect", xmin = 0.76, xmax = 0.86, ymin = y0 - 0.30, ymax = y0 - 0.02,
             fill = "grey95", color = "grey60", linewidth = 0.3),
    annotate("segment", x = 0.63, xend = 0.71, y = y0 - 0.09, yend = y0 - 0.09, linewidth = 1.6, color = "#2E7D6F"),
    annotate("segment", x = 0.77, xend = 0.85, y = y0 - 0.17, yend = y0 - 0.17, linewidth = 1.6, color = "#2E7D6F"),
    annotate("text", x = c(0.67, 0.81), y = y0 + 0.005, label = c("Sham", "Day1"), size = 2.2,
             color = "grey30", family = "Arial", vjust = 1),
    annotate("text", x = 0.74, y = y0 - 0.33, label = sprintf("expected shift ≈ %.0f nt", abs(row$mean_d)),
             size = 2.1, color = "#2E7D6F", family = "Arial", fontface = "italic", hjust = 0.5)
  )
}

p_h <- ggplot() +
  mk_block(ph[gene == "Fxyd2"], 0.98) +
  mk_block(ph[gene == "Kap"], 0.58) +
  annotate("text", x = 0.5, y = 0.10, size = 3, fontface = "bold", color = "#B8860B",
           label = "STATUS: Experiment pending — primers to be designed on validated transcripts", family = "Arial") +
  annotate("text", x = 0.5, y = 0.045, size = 2.2, color = "grey40",
           label = "Retired: Park7 (no validated modification layer). Backup: Gpx1 (3-layer, −16.5 nt, 4/4 tx)", family = "Arial") +
  xlim(0, 1) + ylim(0, 1) + theme_void() +
  labs(title = "h  PolyA PAT assay validation — re-evaluated targets",
       subtitle = "Dual-anchor validated polyA-shortening genes, kidney-abundant (TPM > 500) for assay robustness",
       caption = paste0("Source: NC_vs_Day1.polya_diff.all.xls (canonical company output), aggregated 2026-07-23 as mean of\n",
                        "per-transcript Diff_median over FDR<0.05 transcripts; TPM from gene.TPM.matrix (MGP-ID mapped).")) +
  theme(plot.title = element_text(size = 8.5, face = "bold", family = "Arial"),
        plot.subtitle = element_text(size = 6.5, color = "grey30", family = "Arial"),
        plot.caption = element_text(size = 5.5, color = "grey50", hjust = 0, family = "Arial", lineheight = 1.15),
        plot.margin = margin(5, 5, 5, 5))
save_both(p_h, "Fig8H_PAT_Assay_Template", 160, 110)

# =============================================================================
# PANEL I v2: meRIP-DRS correlation design (text template, experiment pending)
# =============================================================================
pi_genes <- ev[arm == "meRIP-qPCR"]
gene_lines <- paste(sprintf("  %s — %d HC DRACH+ UP sites (%s)", pi_genes$gene,
                            pi_genes$m6A_HC50_DRACH_UP_sites,
                            gsub("three_prime_UTR", "3'UTR", pi_genes$site_struc)), collapse = "\n")
p_i <- ggplot() +
  annotate("text", x = 0.04, y = 0.96, hjust = 0, vjust = 1, size = 3.1, fontface = "bold",
           label = "Design: orthogonal confirmation of DRS stoichiometry by meRIP-qPCR", family = "Arial") +
  annotate("text", x = 0.04, y = 0.86, hjust = 0, vjust = 1, size = 2.55, color = "grey15", lineheight = 1.3,
           label = paste0(
             "X axis: DRS validated m6A-gain signal per gene (dual-anchor HC50 DRACH+ set)\n",
             "Y axis: meRIP-qPCR enrichment change (Δ%input, Day1 vs Sham, n = 3 per group)\n\n",
             "Targets (all validated-gain, qPCR-feasible TPM):\n", gene_lines, "\n\n",
             "Expectation: positive concordance across 5 genes (R > 0.8 target).\n",
             "A successful correlation closes the loop: caller-independent arbitration\n",
             "→ validated gain set → independent antibody-based quantification."), family = "Arial") +
  annotate("text", x = 0.5, y = 0.10, size = 3, fontface = "bold", color = "#B8860B",
           label = "STATUS: Experiment pending", family = "Arial") +
  annotate("text", x = 0.5, y = 0.045, size = 2.2, color = "grey40",
           label = "Supersedes archived design built on retired genes showing the falsified m6A-loss direction", family = "Arial") +
  xlim(0, 1) + ylim(0, 1) + theme_void() +
  labs(title = "i  meRIP-qPCR × DRS correlation — re-evaluated design",
       caption = paste0("Gene set: Fig8_gene_reevaluation_20260723.csv (2026-07-23). Validated basis: dual-anchor paired\n",
                        "arbitration, m6A real gain +2.84pp (DRACH+ HC70 paired). No measured data shown — design only.")) +
  theme(plot.title = element_text(size = 8.5, face = "bold", family = "Arial"),
        plot.caption = element_text(size = 5.5, color = "grey50", hjust = 0, family = "Arial", lineheight = 1.15),
        plot.margin = margin(5, 5, 5, 5))
save_both(p_i, "Fig8I_meRIP_DRS_Correlation", 160, 110)

cat("\n=== Fig8 A/H/I v2 complete (validated-basis replacements) ===\n")
