# =============================================================================
# Fig5 Panels D-I: PolyA + AS-PolyA Integration
# Ported from fig4_as_landscape/scripts/Fig4_fix_real_data.R (Panels H,I,J,K,L)
# and fig4_as_landscape/_archive_pre_v7/scripts_old/Fig4_fix_real_data.R.bak (Panel M)
# as part of the Fig4(AS-only)/Fig5(PolyA+PT) split, 2026-07-13.
# Same underlying real data as the original panels — no new analysis, only
# relettering + (Panel H) compression from 3 sub-panels to 2.
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")  # fix: system default "C" locale mangles unicode (arrows) into "..." in cairo_pdf output
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(ggforce)
  library(scales)
})

source("../../shared/palette/color_palette.R")

data_dir   <- "../../fig4_as_landscape/data"
out_dir    <- "../output"
AS_BASE    <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/07.Structure_Analysis/03.AS/02.Diff_AS"
POLYA_BASE <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/06.polyA"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

norm_gene <- function(x) {
  x <- as.character(x)
  ifelse(x == "-" | x == "" | is.na(x), x,
         paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))))
}

save_std5 <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(out_dir, paste0("Fig5", name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0("Fig5", name, ".png")), p, width = w, height = h, dpi = 300)
}
save_phm5 <- function(p, name, w = 8, h = 5) {
  cairo_pdf(file.path(out_dir, paste0("Fig5", name, ".pdf")), width = w, height = h)
  grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
  png(file.path(out_dir, paste0("Fig5", name, ".png")), width = w, height = h, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
}

theme_uni <- theme_bw(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, color = "grey92"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7))

theme_pub <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
        panel.border = element_rect(color = "grey40", linewidth = 0.4),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8, color = "grey20"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        plot.margin = margin(8, 8, 4, 4))

AS_TYPES <- c("SE", "RI", "A3", "A5", "AF", "AL", "MX")

# =============================================================================
# Shared data loading
# =============================================================================
cat("=== Loading real sequencing data for Fig5 D-I ===\n")

as_dir_d1 <- file.path(AS_BASE, "NC_vs_Day1-diff")
as_all <- rbindlist(lapply(AS_TYPES, function(at) {
  f <- file.path(as_dir_d1, paste0(at, "_gene.dpsi.add_info.xls"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, AS_type := at]
  dt
}), fill = TRUE)
as_all[, gene_sym := norm_gene(gene_name)]
as_all[, `:=`(
  dpsi = as.numeric(gsub(",", "", as.character(dpsi))),
  diffAS_pvalue = as.numeric(gsub(",", "", as.character(`diffAS-pvalue`)))
)]
as_sig <- as_all[!is.na(dpsi) & !is.na(diffAS_pvalue) & diffAS_pvalue < 0.05 & abs(dpsi) > 0.05]
cat(sprintf("DRS AS: %d sig events, %d unique genes (NC vs Day1)\n",
            nrow(as_sig), uniqueN(as_sig$gene_sym)))

polya_d1_f <- file.path(POLYA_BASE, "length", "02.diff", "NC_vs_Day1",
                        "NC_vs_Day1.polya_diff.all.xls")
polya_all <- fread(polya_d1_f)

polya_stat <- fread(file.path(POLYA_BASE, "length", "01.all", "ployA-stat.xls"))
polya_stat[, timepoint := gsub("s[0-9]$", "", sample)]
polya_stat[, timepoint := factor(timepoint,
  levels = c("NC", "I", "Day1", "Day3", "Day7"),
  labels = c("NC", "IS", "Day1", "Day3", "Day7"))]
cat(sprintf("PolyA length stat: %d samples\n", nrow(polya_stat)))

polya_merged <- fread(file.path(data_dir, "PolyA_m6A_merged.csv"))

# =============================================================================
# PANEL D: PolyA Length Panorama (was Fig4H)
# =============================================================================
cat("\n=== Panel D: PolyA Length Panorama ===\n")

fig3k_ramp <- colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(101)
fig3k_bar_colors <- c(
  "NC"   = fig3k_ramp[15],
  "IS"   = fig3k_ramp[35],
  "Day1" = fig3k_ramp[85],
  "Day3" = fig3k_ramp[65],
  "Day7" = fig3k_ramp[78]
)
fig3k_line_color <- fig3k_ramp[10]

polya_stat[, sample := factor(sample, levels = polya_stat$sample[order(timepoint)])]

p_d1 <- ggplot(polya_stat, aes(x = sample, y = Q50, fill = timepoint)) +
  geom_col(width = 0.7, alpha = 0.92, color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Q25, ymax = Q75), width = 0.2, linewidth = 0.6, color = "grey40") +
  scale_fill_manual(values = fig3k_bar_colors,
                    name = "Timepoint",
                    labels = c("NC" = "Normal Control", "IS" = "Ischemia",
                               "Day1" = "Day 1", "Day3" = "Day 3", "Day7" = "Day 7")) +
  labs(x = NULL, y = "Median PolyA Length (nt)",
       title = "PolyA Tail Length Distribution per Sample",
       subtitle = sprintf("IQR bars = Q25-Q75 (real per-sample median); Day1 = %.0f nt vs NC = %.0f nt (d = %.0f nt)",
                          polya_stat[timepoint == "Day1", mean(Q50)],
                          polya_stat[timepoint == "NC", mean(Q50)],
                          polya_stat[timepoint == "NC", mean(Q50)] -
                          polya_stat[timepoint == "Day1", mean(Q50)])) +
  theme_uni + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                     legend.position = "right")

traj_pa <- polya_stat[, .(mean_median = mean(Q50), mean_Q25 = mean(Q25),
                           mean_Q75 = mean(Q75)), by = timepoint]
traj_pa[, timepoint := factor(timepoint, levels = c("NC","IS","Day1","Day3","Day7"))]

p_d2 <- ggplot(traj_pa, aes(x = timepoint, y = mean_median, group = 1)) +
  geom_ribbon(aes(ymin = mean_Q25, ymax = mean_Q75), alpha = 0.18, fill = fig3k_line_color) +
  annotate("text", x = 3, y = max(traj_pa$mean_Q75) * 1.02, label = "shaded: IQR (Q25-Q75)",
           size = 2.8, color = "grey50", hjust = 0.5) +
  geom_line(linewidth = 1.3, color = fig3k_line_color) +
  geom_point(size = 4.5, color = fig3k_line_color, shape = 21, fill = "white", stroke = 1.8) +
  geom_text(aes(label = round(mean_median)), vjust = -1.8, size = 3.8,
            fontface = "bold", color = "grey20") +
  scale_x_discrete(limits = c("NC","IS","Day1","Day3","Day7")) +
  labs(x = NULL, y = "Mean of Per-Sample Median PolyA Length (nt)",
       title = "PolyA Tail Length Trajectory After Kidney IRI",
       subtitle = sprintf("Global shortening at Day1 (-%.0f nt vs NC), gradual recovery by Day7",
                          traj_pa[timepoint == "NC", mean_median] -
                          traj_pa[timepoint == "Day1", mean_median])) +
  theme_uni

p_d <- wrap_plots(p_d1, p_d2, ncol = 1, heights = c(0.55, 0.45))
save_std5(p_d, "D_PolyA_Length_Panorama", 10, 7)
cat("Panel D saved\n")

# =============================================================================
# PANEL E: APA Site Shift (was Fig4J)
# =============================================================================
cat("\n=== Panel E: APA Site Shift ===\n")

apa_comps <- c("NC_vs_Day1", "NC_vs_Day3", "NC_vs_Day7", "NC_vs_I",
               "I_vs_Day1", "I_vs_Day3", "I_vs_Day7")
apa_data <- rbindlist(lapply(apa_comps, function(cmp) {
  d <- file.path(POLYA_BASE, "site", cmp)
  if (!dir.exists(d)) return(NULL)
  ff <- list.files(d, pattern = "normalized\\.all\\.txt$", full.names = TRUE)
  if (length(ff) == 0) ff <- list.files(d, pattern = "normalized\\.txt$", full.names = TRUE)
  if (length(ff) == 0) return(NULL)
  dt <- fread(ff[1])
  cols_keep <- intersect(c("log2FoldChange", "pvalue", "padj", "type"), names(dt))
  if (length(cols_keep) == 0) return(NULL)
  sub <- dt[, ..cols_keep]
  sub[, Comparison := cmp]
  sub
}), fill = TRUE)
cat(sprintf("APA data loaded: %d rows across %d comparisons\n",
            nrow(apa_data), uniqueN(apa_data$Comparison)))

apa_summary <- apa_data[type %in% c("3UTR", "ext_3UTR"),
  .(mean_FC = mean(log2FoldChange, na.rm = TRUE),
    n_sites = .N,
    pct_up = sum(log2FoldChange > 0, na.rm = TRUE) / .N * 100,
    pct_down = sum(log2FoldChange < 0, na.rm = TRUE) / .N * 100),
  by = .(Comparison, type)]
apa_summary[, direction := fifelse(mean_FC > 0, "Lengthen", "Shorten")]
apa_summary[, Comparison := factor(Comparison, levels = apa_comps)]

apa_summary_plot <- apa_summary[Comparison %in% c("NC_vs_Day1", "NC_vs_Day3",
                                                   "NC_vs_Day7", "NC_vs_I")]

p_e <- ggplot(apa_summary_plot, aes(x = Comparison, y = mean_FC, fill = type)) +
  geom_col(position = position_dodge(0.7), width = 0.55, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50", lty = "dashed") +
  scale_fill_manual(values = c("3UTR" = "#E64B35", "ext_3UTR" = "#4DBBD5"),
                    labels = c("3UTR" = "3'UTR PAS", "ext_3UTR" = "Extended 3'UTR"),
                    name = NULL) +
  scale_x_discrete(labels = c("NC_vs_Day1" = "D1", "NC_vs_Day3" = "D3",
                              "NC_vs_Day7" = "D7", "NC_vs_I" = "IS")) +
  labs(x = NULL, y = "Mean log2FC (PAS Usage)",
       title = "APA Site Usage Shift (Real DRS PolyA Site Data)",
       subtitle = paste0("3'UTR PAS shift vs Extended 3'UTR; n=",
                         sum(apa_summary_plot$n_sites), " sites total (4 NC-anchored comparisons shown)")) +
  theme_uni

save_std5(p_e, "E_APA_Site_Shift", 9, 4.5)
cat("Panel E saved\n")

# =============================================================================
# PANEL F: PolyA x Expression decoupling (was Fig4I)
# =============================================================================
cat("\n=== Panel F: PolyA Expression Decoupling ===\n")

polya_expr <- polya_merged[!is.na(log2fc) & !is.na(Effect_size)]
cor_pe <- cor.test(polya_expr$Effect_size, polya_expr$log2fc, method = "spearman")
polya_expr[, abs_eff := abs(Effect_size)]

p_f <- ggplot(polya_expr, aes(x = log2fc, y = Effect_size)) +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70", lty = "dashed") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70", lty = "dashed") +
  geom_point(aes(color = abs_eff), alpha = 0.4, size = 1.5) +
  scale_color_gradient2(low = "#4DBBD5", mid = "grey85", high = "#E64B35",
                         midpoint = median(polya_expr$abs_eff), guide = "none") +
  annotate("text", x = min(polya_expr$log2fc) * 0.85, y = max(polya_expr$Effect_size) * 0.9,
    label = sprintf("rho=%.3f  p=%.2e  n=%d",
                    cor_pe$estimate, cor_pe$p.value, nrow(polya_expr)),
    hjust = 0, size = 3, fontface = "bold") +
  labs(x = "Expression log2FC (Day1 vs NC)", y = "PolyA Length Effect Size",
       title = "PolyA and Expression are Effectively Independent",
       subtitle = sprintf("rho=%.3f, negligible correlation despite large n=%d",
                          cor_pe$estimate, nrow(polya_expr))) +
  theme_pub

save_std5(p_f, "F_PolyA_Expression", 8, 5)
cat(sprintf("Panel F saved: rho=%.3f, n=%d\n", cor_pe$estimate, nrow(polya_expr)))

# =============================================================================
# PANEL G: PolyA x CellType Heatmap (was Fig4K) — PolyA-side of the AS/PolyA
# celltype enrichment table; the AS-side lives in Fig4G. This is the one
# authoritative location for the PolyA version (Fig7's duplicate was cut).
# =============================================================================
cat("\n=== Panel G: PolyA CellType Heatmap ===\n")

enrich_dt <- fread("../../fig7_multiomics/data/enrichment_AS_PolyA_Novel_zscore.csv")
enrich_mat_all <- as.matrix(enrich_dt, rownames = "V1")
polya_z <- enrich_mat_all[c("PolyA_Shorten", "PolyA_Lengthen"), , drop = FALSE]
rownames(polya_z) <- c("PolyA Shorten", "PolyA Lengthen")

p_g <- pheatmap(polya_z,
  cluster_rows = FALSE, cluster_cols = TRUE,
  display_numbers = TRUE, number_format = "%.2f", fontsize_number = 9, fontsize = 8,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  main = "PolyA Length Change Genes x Cell Type Enrichment",
  angle_col = 45, silent = TRUE,
  legend_breaks = c(-2, 0, 2), legend_labels = c("-2", "0", "2"))

save_phm5(p_g, "G_PolyA_CellType_Heatmap", 9, 3.5)
cat("Panel G saved\n")

# =============================================================================
# PANEL H: AS x PolyA Crosstalk (was Fig4L) — COMPRESSED per 2026-07-13 plan:
# keep Venn + effect-size scatter (2 sub-panels); the "Top 15 co-targeted
# genes" lollipop is dropped from the figure and exported as a CSV table
# instead (same information, no third redundant sub-panel).
# =============================================================================
cat("\n=== Panel H: AS x PolyA Crosstalk (compressed to 2 sub-panels) ===\n")

as_genes <- unique(as_sig$gene_sym)
polya_sig_filt <- polya_all[FDR < 0.05 & abs(Effect_size) > 0.3 & Name != "-"]
polya_genes <- unique(norm_gene(polya_sig_filt$Name))

overlap_genes <- intersect(as_genes, polya_genes)
both_n <- length(overlap_genes)

bg_as <- unique(norm_gene(as_all$gene_sym))
bg_polya <- unique(norm_gene(polya_all[Name != "-"]$Name))
all_bg <- length(unique(c(bg_as, bg_polya)))
expected_overlap <- length(as_genes) * length(polya_genes) / all_bg
p_enrich <- phyper(both_n - 1, length(polya_genes), all_bg - length(polya_genes),
                   length(as_genes), lower.tail = FALSE)
or_val <- (both_n / max(1, length(as_genes) - both_n)) /
          (length(polya_genes) / max(1, all_bg - length(polya_genes)))

cat(sprintf("AS=%d | PolyA=%d | overlap=%d | expected=%.0f | OR=%.2f | p=%.2e\n",
            length(as_genes), length(polya_genes), both_n, expected_overlap, or_val, p_enrich))

n1 <- length(as_genes); n2 <- length(polya_genes)
scale_factor <- 0.85
r1 <- sqrt(n1) / sqrt(max(n1, n2)) * scale_factor
r2 <- sqrt(n2) / sqrt(max(n1, n2)) * scale_factor
expected_ratio <- both_n / min(n1, n2)
d <- r1 + r2 - (r1 * expected_ratio * 0.6)

circle_df <- data.frame(
  x = c(-d/2, d/2), y = 0, r = c(r1, r2),
  label = c(paste0("AS\n", n1), paste0("PolyA\n", n2)),
  color = c("#4DBBD5", "#9C78B8")
)

p_h_left <- ggplot() +
  ggforce::geom_circle(data = circle_df, aes(x0 = x, y0 = y, r = r, fill = color), alpha = 0.3) +
  ggforce::geom_circle(data = circle_df, aes(x0 = x, y0 = y, r = r, color = color),
                       fill = NA, linewidth = 1.6) +
  geom_text(data = circle_df, aes(x = x, y = y + r + 0.25, label = label),
            size = 4.2, fontface = "bold", lineheight = 0.9, color = "grey20") +
  annotate("text", x = 0, y = 0.12, label = both_n, size = 6, fontface = "bold", color = "#E64B35") +
  annotate("text", x = 0, y = -0.15, label = "shared", size = 3, color = "grey50") +
  annotate("text", x = 0, y = -0.65,
           label = sprintf("OR = %.2f\np = %.2e", or_val, p_enrich),
           size = 3.2, fontface = "bold", color = "grey30", lineheight = 0.9) +
  annotate("text", x = 0, y = -1.0,
           label = sprintf("%.1fx enrichment", both_n / expected_overlap),
           size = 3, color = "#E64B35", fontface = "bold") +
  annotate("text", x = 0, y = -1.5,
           label = paste0("Background: ", format(all_bg, big.mark = ","), " genes"),
           size = 2.8, color = "grey60") +
  scale_fill_identity() + scale_color_identity() +
  coord_fixed(ylim = c(-1.8, 1.8), xlim = c(-1.5, 1.5)) +
  theme_void()

as_gene_dpsi <- as_sig[, .(max_abs_dpsi = max(abs(dpsi))), by = gene_sym]
polya_eff <- polya_sig_filt[, .(max_abs_eff = max(abs(Effect_size))), by = .(Name_up = norm_gene(Name))]
scatter_data <- merge(as_gene_dpsi, polya_eff, by.x = "gene_sym", by.y = "Name_up", all.x = TRUE)
scatter_data[is.na(max_abs_eff), max_abs_eff := 0]
scatter_data[, group := "Neither"]
scatter_data[gene_sym %in% as_genes & !gene_sym %in% polya_genes, group := "AS Only"]
scatter_data[!gene_sym %in% as_genes & gene_sym %in% polya_genes, group := "PolyA Only"]
scatter_data[gene_sym %in% overlap_genes, group := "Both"]
scatter_data[, group := factor(group, levels = c("AS Only", "PolyA Only", "Both", "Neither"))]
scatter_data[, label := ifelse(group == "Both", gene_sym, "")]
top_both <- scatter_data[group == "Both"][order(-max_abs_dpsi)][1:12]

p_h_right <- ggplot(scatter_data, aes(x = max_abs_dpsi, y = max_abs_eff, color = group)) +
  geom_hline(yintercept = 0.3, linewidth = 0.35, lty = "dashed", color = "grey60") +
  annotate("text", x = max(scatter_data$max_abs_dpsi) * 0.95, y = 0.33, label = "|Effect| = 0.3",
           size = 2.6, color = "grey60", hjust = 1) +
  geom_vline(xintercept = 0.05, linewidth = 0.35, lty = "dashed", color = "grey60") +
  annotate("text", x = 0.055, y = max(scatter_data$max_abs_eff) * 0.95, label = "|dPSI| = 0.05",
           size = 2.6, color = "grey60", hjust = 0, angle = 90) +
  geom_point(data = scatter_data[group != "Both"], alpha = 0.2, size = 1.0) +
  geom_point(data = scatter_data[group == "Both"], alpha = 0.88, size = 2.8, stroke = 0.5) +
  geom_text_repel(data = top_both, aes(label = label), size = 3.0, max.overlaps = 12,
                  box.padding = 0.35, min.segment.length = 0.1, seed = 42) +
  scale_color_manual(values = c("AS Only" = "#4DBBD5", "PolyA Only" = "#9C78B8",
                                 "Both" = "#E64B35", "Neither" = "grey90"),
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3.5))) +
  labs(x = "Max |dPSI| per Gene (absolute splicing change)",
       y = "Max |PolyA Effect| per Gene (absolute length change)",
       title = "AS vs PolyA Effect Size per Gene",
       subtitle = paste0(length(as_genes), " AS genes x ", length(polya_genes),
                         " PolyA genes: ", both_n, " co-targeted (",
                         round(both_n/expected_overlap, 1), "x enrichment, p=",
                         formatC(p_enrich, format="e", digits=1), ")")) +
  theme_uni + theme(legend.position = "right", legend.title = element_blank(),
                     axis.title = element_text(size = 8))

p_h <- wrap_plots(p_h_left, p_h_right, ncol = 2, widths = c(0.32, 0.68))
save_std5(p_h, "H_AS_PolyA_Crosstalk", 12, 5.5)

# Top co-targeted genes -> supplementary table, not a 3rd sub-panel
top_co_table <- scatter_data[group == "Both"][order(-max_abs_dpsi)][1:15,
                 .(Gene = gene_sym, Max_abs_dPSI = round(max_abs_dpsi, 3),
                   Max_abs_PolyA_Effect = round(max_abs_eff, 3))]
fwrite(top_co_table, file.path(out_dir, "Fig5H_top15_cotargeted_genes.csv"))
cat("Panel H saved (2 sub-panels); top-15 co-targeted genes exported as CSV, not a figure panel\n")

# =============================================================================
# PANEL I: Star Genes Iceberg — SLIM VERSION (2026-07-20)
# 4 sub-plots + legend banner -> single row (iceberg + star-gene dot matrix).
# Demoted to supplement: Top 2-Layer Combos (about the 402 two-layer genes, not
# star genes) and cell-type enrichment heatmap (7/8 rows are single-layer sets,
# redundant with Fig5B PT-specific matrix). Legend banner dropped: dot-matrix
# color-coded column headers are self-explanatory; "Novel" column auto-drops
# (no star gene has a Novel-layer hit).
# =============================================================================
cat("\n=== Panel I: Star Genes Iceberg (slim) ===\n")

# Validated multi-layer table (wave framing retired 2026-07-22)
ml_table <- fread("../../fig7_multiomics/data/Table_MultiLayer_Genes_validated.csv")
ml_table[, N_Layers := as.integer(N_Layers)]

tier_dt <- ml_table[, .(Tier = fcase(
  N_Layers == 1, "1 Layer",
  N_Layers == 2, "2 Layers",
  N_Layers == 3, "3 Layers",
  N_Layers == 4, "4 Layers",
  N_Layers >= 5, "5+ Layers\n(elite)"))]
layer_dist <- tier_dt[, .N, by = Tier]
layer_dist[, Tier := factor(Tier, levels = c("5+ Layers\n(elite)", "4 Layers", "3 Layers", "2 Layers", "1 Layer"))]
layer_dist[, pct := N / sum(N) * 100]
setorder(layer_dist, Tier)
n_ge3 <- ml_table[N_Layers >= 3, .N]; n_ge5 <- ml_table[N_Layers >= 5, .N]
cat(sprintf("MultiLayer(validated): %d total | >=3: %d (%.1f%%) | >=5: %d\n",
            nrow(ml_table), n_ge3, 100*n_ge3/nrow(ml_table), n_ge5))

pyramid_colors <- c("1 Layer" = "#B0B0B0", "2 Layers" = "#E69F5B", "3 Layers" = "#C9774A",
                    "4 Layers" = "#C94A4A", "5+ Layers\n(elite)" = "#8E2F2F")

p_i1 <- ggplot(layer_dist, aes(x = Tier, y = N, fill = Tier)) +
  geom_col(width = 0.65, alpha = 0.92, color = "white", linewidth = 0.3) +
  geom_text(aes(label = scales::comma(N)), vjust = -0.4, size = 4.2, fontface = "bold") +
  geom_text(aes(label = paste0(sprintf("%.1f", pct), "%")), vjust = 1.8, size = 3.2,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = pyramid_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of Genes",
       title = "Multi-Layer Regulation:\nAn Iceberg Structure",
       subtitle = sprintf("%s genes (%.1f%%) in 3+ layers; %d reach 5+",
                          scales::comma(n_ge3), 100*n_ge3/nrow(ml_table), n_ge5)) +
  theme_uni +
  theme(plot.margin = margin(8, 8, 8, 8))

# Elite dot-matrix: genes with >=5 validated layers
star3 <- ml_table[N_Layers >= 5][order(-N_Layers, Gene)]
layer_map <- list("m6A↑"="m6A↑","m5C↓"="m5C↓","PSU↑"="PSU↑","inosine↑"="inosine↑",
                  "Nm↑"="Nm↑","Nm↓"="Nm↓","AS"="AS_diff","PolyA"="PolyA_diff","Novel"="Novel_Tx")
layer_cat_colors <- c("m6A↑"="#D6404C","m5C↓"="#E8753A","PSU↑"="#3C5488","inosine↑"="#00A087",
                      "Nm↑"="#6A3D9A","Nm↓"="#4B2673",
                      "AS"="#984EA3","PolyA"="#FF7F00","Novel"="#6BAED6")

star_parsed <- rbindlist(lapply(seq_len(nrow(star3)), function(i) {
  gene <- star3$Gene[i]
  layers <- strsplit(star3$Layers[i], "; ")[[1]]
  rbindlist(lapply(names(layer_map), function(lm_name) {
    data.table(Gene = gene, Layer = lm_name, Hit = layer_map[[lm_name]] %in% layers)
  }))
}))
star_parsed[, Layer := factor(Layer, levels = rev(names(layer_map)))]
star_parsed[, Gene := factor(Gene, levels = rev(star3$Gene))]
star_hits <- star_parsed[Hit == TRUE]

p_i2 <- ggplot(star_parsed, aes(x = Layer, y = Gene)) +
  geom_vline(xintercept = seq(0.5, 9.5, 1), color = "grey93", linewidth = 0.2) +
  geom_point(data = star_hits, aes(fill = Layer, color = Layer), size = 4.2, shape = 21, stroke = 0.7) +
  scale_fill_manual(values = layer_cat_colors, guide = "none") +
  scale_color_manual(values = c("m6A↑"="#A93226","m5C↓"="#C05A2A","PSU↑"="#2A3A5E","inosine↑"="#00755C",
                                "Nm↑"="#4A2A6B","Nm↓"="#33194D",
                                "AS"="#7D3C98","PolyA"="#D68910","Novel"="#5499C7"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = sprintf("%d Elite Genes: 5-6 Layer Multi-Regulation (validated)", nrow(star3)),
       subtitle = "Solid dot = layer affected; Dag1/Kcnj15/Tlk1 reach all 6 layers") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7, face = "bold", lineheight = 0.9, color = "grey20"),
        axis.text.y = element_text(size = 6, face = "italic"),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7.5, color = "grey50", hjust = 0.5),
        plot.margin = margin(8, 8, 8, 8))

star2 <- ml_table[N_Layers == 2]
combo2 <- star2[, .N, by = Layers][order(-N)][1:min(6, .N)]
combo2[, Layers_short := gsub("; ", " +\n", Layers)]
combo2[, Layers_short := factor(Layers_short, levels = Layers_short)]

p_i3 <- ggplot(combo2, aes(x = Layers_short, y = N, fill = N)) +
  geom_col(width = 0.6, alpha = 0.88, color = "white", linewidth = 0.3) +
  geom_text(aes(label = N), vjust = -0.3, size = 3.8, fontface = "bold") +
  scale_fill_gradientn(colors = c("#B09C85", "#C94A4A"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Number of Genes",
       title = "Top 2-Layer Combinations",
       subtitle = sprintf("Total: %d genes with 2-layer regulation", nrow(star2))) +
  theme_uni +
  theme(axis.text.x = element_text(size = 6.8, lineheight = 0.9, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))

enrich_all <- fread("../../fig7_multiomics/data/enrichment_AS_PolyA_Novel_zscore.csv")
enrich_mat_i <- as.matrix(enrich_all, rownames = "V1")
# 2026-09-23: "m6A_Down" row REMOVED — that set is the pre-arbitration legacy
# call set (m6A loss wave judged artifact by dual-anchor ModVerdict; validated
# framework carries no m6A↓ row, cf. main Fig 2K). Remaining rows are validated
# direction sets (m5C_Down/PSU_Up/Um_Up) + structural-layer sets.
layer_rows_avail <- intersect(c("AS_SE", "PolyA_Shorten", "Multi_Layer",
                "m5C_Down", "PSU_Up", "Um_Up", "Novel_Tx"), rownames(enrich_mat_i))
layer_labels_map <- c("AS_SE"="Alt. Splicing (SE)","PolyA_Shorten"="PolyA Shortening",
                       "Multi_Layer"="3-Layer Combined","m5C_Down"="m5C Down",
                       "PSU_Up"="PSU Up","Um_Up"="Um Up","Novel_Tx"="Novel Transcript")
n_z <- enrich_mat_i[layer_rows_avail, , drop = FALSE]
rownames(n_z) <- layer_labels_map[layer_rows_avail]

p_i4 <- pheatmap(n_z,
  cluster_rows = TRUE, cluster_cols = TRUE,
  display_numbers = TRUE, number_format = "%.1f", fontsize_number = 9, fontsize = 9,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  breaks = seq(-2.5, 3.5, length.out = 101),
  main = "Multi-Layer Gene Set Enrichment Across Cell Types",
  angle_col = 45, silent = TRUE,
  legend_breaks = c(-2, 0, 2, 4), legend_labels = c("-2", "0", "2", "4"),
  border_color = NA)

# --- Main panel: single row, iceberg (narrow) + star-gene dot matrix ---
p_i_main <- wrap_plots(p_i1, p_i2, ncol = 2, widths = c(0.30, 0.70))
ggsave(file.path(out_dir, "Fig5I_StarGenes_MultiLayer.pdf"), p_i_main,
       width = 13, height = 4.8, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig5I_StarGenes_MultiLayer.png"), p_i_main,
       width = 13, height = 4.8, dpi = 600)

# --- Supplement: demoted sub-plots exported standalone ---
ggsave(file.path(out_dir, "Fig5I_supp_2Layer_Combos.pdf"), p_i3,
       width = 8, height = 4, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig5I_supp_2Layer_Combos.png"), p_i3,
       width = 8, height = 4, dpi = 600)
save_phm5(p_i4, "I_supp_MultiLayer_CellType_Enrichment", 10, 5)

cat(sprintf("Panel I (slim) saved: %d total genes, %d elite genes (5+ layers); combos + enrichment -> supplement\n",
            nrow(ml_table), nrow(star3)))

cat("\n=== Fig5 D-I complete ===\n")
