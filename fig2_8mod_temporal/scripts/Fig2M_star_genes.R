#!/usr/bin/env Rscript
# Fig2 Panel M: Star Gene Multi-Modification Tracks
# ============================================================
# RE-ANCHORED 2026-07-22 (wave framing retired): validated star genes
#   Tmem176b: Nm turnover — Um +4.97pp (100% sites up) vs Gm -1.26pp (paired)
#   Ggt1:     PT brush-border marker — Cm/Um/inosine gains, 5-layer elite
#   Sec31a:   inosine +2.7pp (87% up) + Um +8.6pp, 5-layer elite
# Park7/Fis1/Malat1 dropped: zero membership in any validated direction set.
#
# Without BAM files, use:
#   - Expression trajectory (5-group violin/point)
#   - Modification fraction change (NC→Day1, 4 mods)
#   - Gene-level summary: m6A vs activating mod ratio

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("../../shared/palette/color_palette.R")

mod_display_map <- c(m6a="m6A", m5c="m5C", psu="PSU", inosine="Inosine",
                     Am="Am", Cm="Cm", Gm="Gm", Um="Um")

out_dir <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load expression data for star genes
# ============================================================
cat("Loading expression data...\n")

expr_file <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/03.expression/gene/gene.TPM.matrix.add_info.xls"
expr_all <- fread(expr_file)

# Star genes
star_genes <- list(
  Tmem176b = "MGP_C57BL6NJ_G0031149",
  Ggt1     = "MGP_C57BL6NJ_G0017753",
  Sec31a   = "MGP_C57BL6NJ_G0030216"
)

sample_order <- c("NCs1","NCs2","NCs3",
                  "Is1","Is2","Is3",
                  "Day1s1","Day1s2","Day1s3",
                  "Day3s1","Day3s2","Day3s3",
                  "Day7s1","Day7s2","Day7s3")

expr_star <- expr_all[ID %in% unlist(star_genes)]
# Extract gene name
expr_star[, gene_name := names(star_genes)[match(ID, star_genes)]]

# Melt to long format
expr_long <- melt(expr_star, id.vars = c("ID", "gene_name"),
                  measure.vars = sample_order,
                  variable.name = "Sample", value.name = "TPM")
expr_long[, group := gsub("s[0-9]$", "", Sample)]
expr_long[group == "I", group := "IS"]
expr_long[, group := factor(group, levels = group_order)]
expr_long[, gene_name := factor(gene_name, levels = names(star_genes))]

# ============================================================
# 2. Load modification data (NC vs Day1, per-transcript)
# ============================================================
cat("Loading modification data...\n")

# Canonical transcripts for each gene
canonical <- list(
  Tmem176b = "ENSMUSG00000029810.t1",
  Ggt1     = "ENSMUSG00000006345.t14",
  Sec31a   = "ENSMUSG00000035325.t5"
)

mod_data <- rbindlist(lapply(c("m6a","m5c","psu","inosine"), function(mod) {
  f <- sprintf("/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/%s/trans_diff/NC_vs_Day1/NC_vs_Day1.%s.per_trans.add_name.xls.gz",
               mod, mod)
  if (!file.exists(f)) return(NULL)
  dt <- fread(cmd = paste("zcat", f))
  dt[, mod := mod_display_map[mod]]
  dt[, gene := sub("\\..*", "", Trans)]  # Strip .t suffix
  return(dt)
}))

# Filter for canonical transcripts
mod_star <- mod_data[Trans %in% unlist(canonical)]
mod_star[, gene_name := names(canonical)[match(Trans, canonical)]]
setnames(mod_star, "name", "GeneSymbol")

cat(sprintf("Loaded %d modification records for star genes\n", nrow(mod_star)))

# ============================================================
# 3. Compute per-group mean expression
# ============================================================
expr_summary <- expr_long[, .(mean_TPM = mean(TPM), se_TPM = sd(TPM)/sqrt(.N)),
                          by = .(gene_name, group)]

# ============================================================
# Panel F1: Expression trajectory
# ============================================================
cat("\nDrawing Panel F1: Expression trajectory...\n")

p_expr <- ggplot(expr_long, aes(x = group, y = TPM, fill = group)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_point(aes(color = group), position = position_jitter(width = 0.1),
             size = 1.5, alpha = 0.6) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.8, color = "grey30") +
  stat_summary(fun = mean, geom = "point", size = 2.5, shape = 18, color = "grey30") +
  facet_wrap(~gene_name, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  labs(x = NULL, y = "TPM",
       title = "Star Gene Expression Trajectory",
       subtitle = "Tmem176b: high (14-41 TPM) | Ggt1: moderate (4-18) | Sec31a: low (0.7-2.0)") +
  theme_drs

# ============================================================
# Panel F2: Modification fraction change (NC → Day1)
# ============================================================
cat("Drawing Panel F2: Modification fraction change...\n")

mod_plot <- mod_star[, .(gene_name, mod,
                          NC_pct = Group1_pct_modified * 100,
                          Day1_pct = Group2_pct_modified * 100,
                          log2fc, diffType)]

mod_long <- melt(mod_plot, id.vars = c("gene_name", "mod", "log2fc", "diffType"),
                 measure.vars = c("NC_pct", "Day1_pct"),
                 variable.name = "Group", value.name = "PctModified")
mod_long[, Group := gsub("_pct", "", Group)]
mod_long[, mod := factor(mod, levels = c("m6A","m5C","PSU","Inosine"))]
mod_long[, gene_name := factor(gene_name, levels = names(star_genes))]

p_mod <- ggplot(mod_long, aes(x = mod, y = PctModified, fill = Group)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.85) +
  facet_wrap(~gene_name, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(NC = "#4DBBD5", Day1 = "#E64B35"),
                    labels = c("NC", "Day1")) +
  labs(x = NULL, y = "% Modified",
       title = "Modification Fraction: NC vs Day1",
       subtitle = "Canonical transcript per gene") +
  theme_drs

# M2 → Supplementary figure (moved per redundancy audit 2026-06-20)
pdf(file.path(out_dir, "FigSx_StarGenes_ModFraction.pdf"), width = 10, height = 5)
print(p_mod)
dev.off()
png(file.path(out_dir, "FigSx_StarGenes_ModFraction.png"), width = 10, height = 5, units = "in", res = 300)
print(p_mod)
dev.off()
cat("  Saved (Supp): FigSx_StarGenes_ModFraction.{pdf,png}\n")

# ============================================================
# Panel F3: Per-gene modification shifts (VALIDATED paired sites)
# ============================================================
cat("Drawing Panel F3: validated per-gene paired deltas...\n")

pd_dt <- fread("../data/Fig2F3_per_gene_paired_deltas.tsv")
MOD_ORDER <- c("m6A","m5C","PSU","inosine","Um","Am","Cm","Gm")
MOD_CLASS <- c("m6A"="m6A","m5C"="m5C","PSU"="PSU","inosine"="inosine",
               "Um"="2'-O-Me (Nm)","Am"="2'-O-Me (Nm)","Cm"="2'-O-Me (Nm)","Gm"="2'-O-Me (Nm)")
pd_dt <- pd_dt[gene %in% names(star_genes)]
pd_dt[, mod := factor(mod, levels = MOD_ORDER)]
pd_dt[, class := factor(MOD_CLASS[as.character(mod)], levels = unique(MOD_CLASS))]
pd_dt[, gene := factor(gene, levels = names(star_genes))]
pd_dt[, lab := sprintf("%+.1fpp\n(%d%% up, n=%d)", median_delta_pp, round(pct_up), n_paired)]

p_bidirect <- ggplot(pd_dt, aes(x = mod, y = gene, fill = median_delta_pp)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = lab), size = 2.6, color = "grey15", lineheight = 0.9) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#D6404C", midpoint = 0,
                       name = "Median paired\nfrac delta (pp)") +
  facet_wrap(~class, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = NULL,
       title = "Per-Gene Modification Shifts: Validated Paired Sites (Sham -> Day1)",
       subtitle = "Caller-independent paired sites (3/3 reps each condition); Tmem176b: Um up / Gm down = Nm turnover") +
  theme_drs +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold", size = 9))

# ============================================================
# Panel F4: Gene structure schematic (simplified)
# ============================================================
cat("Drawing Panel F4: Gene structure schematic...\n")

# Load m6A site data per gene (estimated from all_m6a.rds summary)
# Since RDS is 797MB, use trans_diff summary data
# Extract total m6A sites per transcript from the per_trans files
gene_summary <- rbindlist(lapply(names(canonical), function(gn) {
  trans <- canonical[[gn]]

  # Get site count from per_trans file
  mod_info <- mod_star[Trans == trans, .(mod, Group1_counts, Group2_counts,
                                          Group1_total, Group2_total)]
  mod_info[, gene_name := gn]

  # Extract site count from "prefix:number" format (a:/m:/17802:/17596:)
  extract_count <- function(x) {
    x <- sub("^[a-zA-Z0-9]+:", "", x)
    x[x == "." | x == ""] <- "0"
    as.numeric(x)
  }
  nc_sites <- sum(extract_count(mod_info$Group1_counts), na.rm = TRUE)
  d1_sites <- sum(extract_count(mod_info$Group2_counts), na.rm = TRUE)
  nc_coverage <- max(mod_info$Group1_total, na.rm = TRUE)
  d1_coverage <- max(mod_info$Group2_total, na.rm = TRUE)

  data.table(
    gene_name = gn,
    transcript = trans,
    m6A_sites_NC = nc_sites,
    m6A_sites_Day1 = d1_sites,
    coverage_NC = nc_coverage,
    coverage_Day1 = d1_coverage,
    m6A_density_NC = nc_sites / nc_coverage * 1000,
    m6A_density_Day1 = d1_sites / d1_coverage * 1000
  )
}))

# Combine with expression (per group mean)
expr_wide <- dcast(expr_summary, gene_name ~ group, value.var = "mean_TPM")
expr_wide[, mean_expr := rowMeans(.SD, na.rm = TRUE), .SDcols = group_order]
expr_wide <- expr_wide[, .(gene_name, mean_expr)]

gene_summary <- merge(gene_summary, expr_wide, by = "gene_name")

# Multi-panel gene card (lollipop only — expression bar removed per redundancy audit,
# M1 already provides full expression trajectory)
gene_cards <- lapply(names(star_genes), function(gn) {
  gd <- gene_summary[gene_name == gn]
  md <- mod_plot[gene_name == gn]

  ggplot(md, aes(x = mod, y = log2fc, color = mod)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(xend = mod, yend = 0), linewidth = 1.5) +
    geom_point(size = 5) +
    scale_color_manual(values = all_8mod_colors, guide = "none") +
    labs(x = NULL, y = "log2FC (Day1/NC)", title = gn,
         subtitle = sprintf("m6A density: %.1f → %.1f sites/kb",
                           gd$m6A_density_NC, gd$m6A_density_Day1)) +
    theme_drs +
    theme(plot.title = element_text(face = "bold", size = 12))
})

p_cards <- wrap_plots(gene_cards, nrow = 1)

# ============================================================
# Combine: M1(expression) + M3(bidirectional heatmap) + M4(gene lollipops)
# M2(mod fraction bar) moved to Supplementary
# ============================================================
p_f_combined <- p_expr / p_bidirect / p_cards +
  plot_layout(heights = c(1.0, 0.6, 0.8)) +
  plot_annotation(
    title = "Figure 2M: Star Gene Multi-Modification Architecture",
    subtitle = "Tmem176b (Nm turnover) | Ggt1 (PT marker, Cm/Um gain) | Sec31a (inosine/Um gain) — validated paired-site basis",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

pdf(file.path(out_dir, "Fig2M_StarGenes.pdf"), width = 12, height = 10)
print(p_f_combined)
dev.off()
png(file.path(out_dir, "Fig2M_StarGenes.png"), width = 12, height = 10, units = "in", res = 300)
print(p_f_combined)
dev.off()
png(file.path(out_dir, "Fig2M_StarGenes.png"), width = 14, height = 16, units = "in", res = 300)
print(p_f_combined)
dev.off()
cat("  Saved: Fig2M_StarGenes.{pdf,png}\n")

# ============================================================
# Export data
# ============================================================
fwrite(gene_summary, file.path(out_dir, "Table_Fig2M_StarGenes_Summary.csv"), sep = ",")
fwrite(mod_plot, file.path(out_dir, "Table_Fig2M_StarGenes_Modifications.csv"), sep = ",")

cat("\n=== Fig2 Panel M Complete ===\n")
cat("Key findings:\n")
for (gn in names(star_genes)) {
  md <- mod_plot[gene_name == gn]
  m6a_fc <- md[mod == "m6A", log2fc]
  others_up <- md[mod != "m6A" & log2fc > 0, .N]
  cat(sprintf("  %s: m6A log2FC=%.2f, %d/3 other mods ↑ (bidirectional=%s)\n",
              gn, m6a_fc, others_up, ifelse(others_up >= 2, "YES", "partial")))
}
