#!/usr/bin/env Rscript
# Fig1 Panels C, E, F: DRS Platform QC & Data Quality
# ============================================================
# Panel C: N50 violin + Mapping rate bar + Read length density
# Panel E: PCA + 8-mod per-sample site count correlation heatmap
# Panel F: m6A site count bar + Jaccard stability + 8-mod CV%

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})

source("../../shared/palette/color_palette.R")

# ============================================================
# Data paths
# ============================================================
data_align    <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/02.alignment"
data_expr     <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/03.expression"
data_methy    <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
out_dir       <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_order <- c("NCs1","NCs2","NCs3",
                  "Is1","Is2","Is3",
                  "Day1s1","Day1s2","Day1s3",
                  "Day3s1","Day3s2","Day3s3",
                  "Day7s1","Day7s2","Day7s3")

mods_8 <- c("m6a","m5c","psu","inosine","Am","Cm","Gm","Um")
mod_display_map <- c(m6a="m6A", m5c="m5C", psu="PSU", inosine="Inosine",
                     Am="Am", Cm="Cm", Gm="Gm", Um="Um")

# ============================================================
# 1. Load QC data
# ============================================================
cat("Loading QC data...\n")

# Collapse stats (N50)
# [2026-07-20 fix] File uses comma as thousands separator ("3,926"); fread
# auto-detects dec="," and parses it as 3.926. Force dec="." so the field
# stays character, then strip commas and convert bp -> kb.
collapse_dt <- fread(file.path(data_align, "collapse_stats.xls"), dec = ".")
setnames(collapse_dt, old = "Sample Name", new = "Sample")
collapse_dt[, Sample := gsub(" ", "", Sample)]
num_clean <- function(x) as.numeric(gsub("[, ]", "", as.character(x)))
collapse_dt[, `:=`(
  N50_kb     = num_clean(`N50(bp)`) / 1000,
  MeanLen_kb = num_clean(`Mean Length(bp)`) / 1000,
  SeqNum     = num_clean(`Seq Num`)
)]
collapse_dt[, group := gsub("s[0-9]$", "", Sample)]
collapse_dt[group == "I", group := "IS"]

# [2026-07-20] Unified figure naming: NC -> Sham (matches composed figures);
# defined locally so the global palette's NC naming stays untouched.
# [2026-07-20 v2] Colors extracted from composed fig1.ai (NPG family):
#   Sham #E64B35 / IS #F39B7F / Day1 #4DBBD5 / Day3 #00A087 / Day7 #3C5488
sham_levels <- c("Sham", "IS", "Day1", "Day3", "Day7")
sham_colors <- c(Sham = "#E64B35", IS = "#F39B7F", Day1 = "#4DBBD5",
                 Day3 = "#00A087", Day7 = "#3C5488")
# Classic theme matching the composed Fig1 panels (white bg, axis lines, Arial)
theme_fig1c <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(axis.line  = element_line(linewidth = 0.4, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.4, colour = "grey30"),
        axis.ticks.length = unit(2.5, "pt"),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.title.position = "panel")
collapse_dt[, group := factor(ifelse(group == "NC", "Sham", as.character(group)),
                              levels = sham_levels)]

# Alignment stats (mapping rate) — dec="." for the same reason as above
align_dt <- fread(file.path(data_align, "align_stats.xls"), fill = TRUE, dec = ".")
align_dt[, Sample := gsub(" ", "", Sample)]
align_dt[, map_rate := as.numeric(gsub("%", "", `map_rate`))]

# Merge QC
qc_dt <- merge(collapse_dt[, .(Sample, group, N50_kb, MeanLen_kb, SeqNum)],
               align_dt[, .(Sample, `Total reads`, `mapped reads`, map_rate)],
               by = "Sample")

cat(sprintf("Loaded QC for %d samples\n", nrow(qc_dt)))

# ============================================================
# 2. Load 8-mod site counts
# ============================================================
cat("Loading 8-mod site counts...\n")

site_counts <- matrix(NA, nrow = length(mods_8), ncol = length(sample_order))
rownames(site_counts) <- mods_8
colnames(site_counts) <- sample_order

for (m in mods_8) {
  fpath <- file.path(data_methy, m, "methy_detect", paste0("all_sample.", m, ".0.7.xls"))
  if (!file.exists(fpath)) {
    cat(sprintf("  [WARN] %s: file not found: %s\n", m, fpath))
    next
  }
  dt <- fread(fpath, sep = "\t", showProgress = FALSE)
  if (!("type" %in% names(dt))) next
  for (i in seq_len(nrow(dt))) {
    if (dt$type[i] %in% sample_order) {
      site_counts[m, dt$type[i]] <- dt$number[i]
    }
  }
  cat(sprintf("  [OK] %s: %d sites (mean per sample)\n", m, round(mean(dt$number, na.rm = TRUE))))
}

# ============================================================
# Panel C1: N50 by group — dots + median crossbar (n=3/group, no violin)
# ============================================================
cat("\n=== Panel C: N50 + Mapping Rate + Read Length ===\n")

p_n50 <- ggplot(qc_dt, aes(x = group, y = N50_kb, color = group)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4,
               linewidth = 0.5, color = "grey30") +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.9, shape = 16) +
  scale_color_manual(values = sham_colors, guide = "none") +
  labs(x = NULL, y = "N50 (kb)", title = "Read N50 by Group") +
  theme_fig1c

# Panel C2: Mapping rate — dots + mean +/- SD
p_map <- ggplot(qc_dt, aes(x = group, y = map_rate, color = group)) +
  stat_summary(fun.min = function(y) mean(y) - sd(y),
               fun.max = function(y) mean(y) + sd(y),
               geom = "errorbar", width = 0.15, linewidth = 0.4, color = "grey30") +
  stat_summary(fun = mean, geom = "point", size = 3, color = "grey30") +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.9, shape = 16) +
  scale_color_manual(values = sham_colors, guide = "none") +
  scale_y_continuous(limits = c(85, 100)) +
  labs(x = NULL, y = "Mapping Rate (%)", title = "Mapping Rate by Group") +
  theme_fig1c

# Panel C3: real read-length distribution from per-sample collapse isoform
# FASTA (precomputed to data/Fig1C_isoform_lengths.csv; N50 cross-validated
# against collapse_stats — exact match). Replaces the old fake density of
# 3 summary values which rendered as an empty panel.
len_dt <- fread("../data/Fig1C_isoform_lengths.csv")
len_dt[, group := gsub("s[0-9]$", "", Sample)]
len_dt[group == "I", group := "IS"]
len_dt[, group := factor(ifelse(group == "NC", "Sham", as.character(group)),
                         levels = sham_levels)]
len_dt[, length_kb := length / 1000]

p_len <- ggplot(len_dt, aes(x = length_kb, color = group, fill = group)) +
  geom_density(alpha = 0.12, linewidth = 0.5) +
  scale_color_manual(values = sham_colors) +
  scale_fill_manual(values = sham_colors) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "Collapsed isoform length (kb)", y = "Density",
       title = "Read Length Distribution by Group",
       color = NULL, fill = NULL) +
  theme_fig1c + theme(legend.position = "bottom")

# Combine C panels
p_c_combined <- (p_n50 | p_map | p_len) +
  plot_annotation(title = "Figure 1C: DRS RNA Integrity & Data Quality",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))

cairo_pdf(file.path(out_dir, "Fig1C_N50_Map_ReadLen.pdf"), width = 14, height = 5)
print(p_c_combined)
dev.off()
png(file.path(out_dir, "Fig1C_N50_Map_ReadLen.png"), width = 14, height = 5, units = "in", res = 600)
print(p_c_combined)
dev.off()
cat("  Saved: Fig1C_N50_Map_ReadLen.{pdf,png}\n")

# ============================================================
# Panel E1: Gene-level PCA
# ============================================================
cat("\n=== Panel E: PCA + Site Count Correlation ===\n")

pca_gene <- fread(file.path(data_expr, "gene", "All_samples.PCA.pcx.txt"))
pca_dt <- melt(pca_gene, id.vars = "V1", variable.name = "PC", value.name = "coord")
pca_wide <- dcast(pca_dt, V1 ~ PC, value.var = "coord")
setnames(pca_wide, "V1", "Sample")
pca_wide[, group := gsub("s[0-9]$", "", Sample)]
pca_wide[group == "I", group := "IS"]
pca_wide[, group := factor(ifelse(group == "NC", "Sham", as.character(group)),
                           levels = sham_levels)]

pvar <- fread(file.path(data_expr, "gene", "All_samples.PCA.proportion.txt"))
pc1_var <- round(as.numeric(pvar[V1 == "Proportion of Variance", PC1]) * 100, 1)
pc2_var <- round(as.numeric(pvar[V1 == "Proportion of Variance", PC2]) * 100, 1)

p_pca <- ggplot(pca_wide, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  # [2026-07-20] stat_ellipse removed: n=3/group is too few to compute an
  # ellipse (silently renders nothing + warnings); composed-figure ellipses
  # around Sham/IS are drawn in Illustrator if needed
  geom_text_repel(aes(label = Sample), size = 2.8, show.legend = FALSE, max.overlaps = 15) +
  scale_color_manual(values = sham_colors, name = "Group") +
  scale_fill_manual(values = sham_colors, guide = "none") +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       title = "Gene-Level PCA (15 samples)") +
  theme_fig1c

# ============================================================
# Panel E2: 8-mod per-sample site count correlation heatmap
# ============================================================
# Transpose to samples × mods, compute Spearman between samples
site_mat <- t(site_counts)
colnames(site_mat) <- mod_display_map[colnames(site_mat)]

# Remove mods with too much missing data
valid_mods <- colnames(site_mat)[colSums(!is.na(site_mat)) > 0]
site_mat_clean <- site_mat[, valid_mods, drop = FALSE]

# Sample-sample correlation (15x15)
cor_samples <- cor(t(site_mat_clean), method = "spearman", use = "pairwise.complete.obs")

# Group annotation (NC -> Sham unified naming + composed-figure colors)
sample_groups <- gsub("s[0-9]$", "", rownames(cor_samples))
sample_groups[sample_groups == "I"] <- "IS"
sample_groups[sample_groups == "NC"] <- "Sham"
names(sample_groups) <- rownames(cor_samples)

annot_col <- data.frame(Group = sample_groups, row.names = rownames(cor_samples))
annot_colors_list <- list(Group = sham_colors)

# Also compute site count matrix for export
site_long <- as.data.table(site_mat_clean, keep.rownames = "Sample")
site_long[, group := gsub("s[0-9]$", "", Sample)]
site_long[group == "I", group := "IS"]

p_heat <- pheatmap(cor_samples,
                   annotation_col = annot_col,
                   annotation_row = annot_col,
                   annotation_colors = annot_colors_list,
                   display_numbers = TRUE,
                   number_format = "%.2f",
                   number_color = "grey40",
                   fontsize_number = 7,
                   color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
                   border_color = NA,
                   clustering_method = "ward.D2",
                   main = "Sample-Sample Correlation (8-mod site counts, Spearman)",
                   angle_col = 45,
                   fontsize = 8,
                   silent = TRUE)

cairo_pdf(file.path(out_dir, "Fig1E_SampleCor_Heatmap.pdf"), width = 10, height = 9)
grid::grid.newpage()
grid::grid.draw(p_heat$gtable)
dev.off()
png(file.path(out_dir, "Fig1E_SampleCor_Heatmap.png"), width = 10, height = 9, units = "in", res = 600)
grid::grid.newpage()
grid::grid.draw(p_heat$gtable)
dev.off()
cat("  Saved: Fig1E_SampleCor_Heatmap.{pdf,png}\n")

# PCA + heatmap combined
p_e_combined <- p_pca +
  plot_annotation(title = "Figure 1E: Biological Reproducibility",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))

cairo_pdf(file.path(out_dir, "Fig1E_Gene_PCA.pdf"), width = 8, height = 7)
print(p_e_combined)
dev.off()
png(file.path(out_dir, "Fig1E_Gene_PCA.png"), width = 8, height = 7, units = "in", res = 600)
print(p_e_combined)
dev.off()
cat("  Saved: Fig1E_Gene_PCA.{pdf,png}\n")

# ============================================================
# Panel F1: m6A per-sample site count bar
# ============================================================
cat("\n=== Panel F: Site Count Stability ===\n")

m6a_counts <- data.table(
  Sample = sample_order,
  Sites = site_counts["m6a", ],
  group = gsub("s[0-9]$", "", sample_order)
)
m6a_counts[group == "I", group := "IS"]
m6a_counts[, group := factor(ifelse(group == "NC", "Sham", as.character(group)),
                             levels = sham_levels)]

p_m6abar <- ggplot(m6a_counts, aes(x = Sample, y = Sites, fill = group)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = sham_colors, guide = "none") +
  labs(y = "m6A Site Count", title = "m6A Detection per Sample (Fraction > 0.7)") +
  theme_fig1c +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ============================================================
# Panel F2: NC internal vs NC-D1 Jaccard by modification
# ============================================================
cat("Computing per-mod Jaccard stability...\n")

compute_mod_jaccard <- function(row_vec, all_sample_names) {
  nc_idx <- grep("^NC", all_sample_names)
  d1_idx <- grep("^Day1", all_sample_names)

  if (length(nc_idx) < 2) return(c(NA, NA))

  nc_pairs <- combn(nc_idx, 2)
  nc_jac <- mean(sapply(seq_len(ncol(nc_pairs)), function(i) {
    a <- nc_pairs[1, i]; b <- nc_pairs[2, i]
    inter <- min(row_vec[a], row_vec[b])
    union <- max(row_vec[a], row_vec[b])
    if (union == 0) return(NA) else inter / union
  }), na.rm = TRUE)

  if (length(d1_idx) == 0) return(c(nc_jac, NA))
  nc_d1_jac <- mean(sapply(nc_idx, function(nc) {
    mean(sapply(d1_idx, function(d1) {
      inter <- min(row_vec[nc], row_vec[d1])
      union <- max(row_vec[nc], row_vec[d1])
      if (union == 0) return(NA) else inter / union
    }), na.rm = TRUE)
  }), na.rm = TRUE)
  c(nc_jac, nc_d1_jac)
}

jac_results <- t(sapply(seq_len(nrow(site_counts)), function(i) {
  compute_mod_jaccard(site_counts[i, ], colnames(site_counts))
}))
colnames(jac_results) <- c("NC_Internal", "NC_Day1")
rownames(jac_results) <- mods_8

jac_dt <- as.data.table(jac_results, keep.rownames = "Mod")
jac_dt[, Mod_display := mod_display_map[Mod]]
jac_dt <- jac_dt[!is.na(NC_Internal)]
jac_dt[, Mod_display := factor(Mod_display, levels = intersect(all_8mod_order, jac_dt$Mod_display))]

jac_long <- melt(jac_dt, id.vars = "Mod_display", measure.vars = c("NC_Internal", "NC_Day1"),
                 variable.name = "Comparison", value.name = "Jaccard")

p_jaccard <- ggplot(jac_long, aes(x = Mod_display, y = Jaccard, fill = Comparison)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.85) +
  # [2026-07-20] recolored: old #4DBBD5/#E64B35 now collide with Day1/Sham
  # group colors in the composed Fig1 — use non-group hues instead
  scale_fill_manual(values = c(NC_Internal = "#8491B4", NC_Day1 = "#D6404C"),
                    labels = c("Sham Internal", "Sham vs Day1")) +
  labs(x = NULL, y = "Jaccard Similarity",
       title = "Site Detection Stability: Sham Internal vs Sham-Day1",
       subtitle = "Per-modification Jaccard index") +
  theme_fig1c

# ============================================================
# Panel F3: 8-mod per-sample CV%
# ============================================================
cat("Computing per-mod CV%...\n")

cv_mod <- apply(site_counts, 1, function(x) {
  x_clean <- x[!is.na(x) & x > 0]
  if (length(x_clean) < 3) return(NA)
  sd(x_clean) / mean(x_clean) * 100
})

cv_dt <- data.table(
  Mod = mod_display_map[names(cv_mod)],
  CV_pct = as.numeric(cv_mod)
)
cv_dt <- cv_dt[!is.na(CV_pct)]
cv_dt[, Mod := factor(Mod, levels = intersect(all_8mod_order, cv_dt$Mod))]

# Assign mod type colors
cv_dt[, mod_type := ifelse(Mod %in% c("m6A","m5C","PSU","Inosine"),
                            "Base modification", "2'-O-Me")]
cv_dt[, fill_color := ifelse(mod_type == "Base modification", "#E64B35", "#984EA3")]

p_cv <- ggplot(cv_dt, aes(x = Mod, y = CV_pct, fill = mod_type)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", CV_pct)), vjust = -0.4, size = 3.2) +
  # [2026-07-20] recolored to the paper's modification palette family
  # (m6A red #D6404C / Nm teal #2E7D6F) instead of Sham-red/purple
  scale_fill_manual(values = c("Base modification" = "#D6404C", "2'-O-Me" = "#2E7D6F"),
                    name = "Modification Type") +
  labs(x = NULL, y = "Coefficient of Variation (%)",
       title = "Per-Modification Site Count Variability (CV%)",
       subtitle = "Across 15 samples; higher CV = more biological/capture variation") +
  theme_fig1c +
  theme(axis.text.x = element_text(face = "bold"))

# ============================================================
# Combine Panel F sub-panels
# ============================================================
p_f_top <- p_m6abar | p_jaccard
p_f_combined <- p_f_top / p_cv +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(title = "Figure 1F: Site Detection Stability & Reproducibility",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))

cairo_pdf(file.path(out_dir, "Fig1F_Jaccard_Stability.pdf"), width = 14, height = 10)
print(p_f_combined)
dev.off()
png(file.path(out_dir, "Fig1F_Jaccard_Stability.png"), width = 14, height = 10, units = "in", res = 600)
print(p_f_combined)
dev.off()
cat("  Saved: Fig1F_Jaccard_Stability.{pdf,png}\n")

# ============================================================
# Export data tables
# ============================================================
cat("\nExporting data tables...\n")

fwrite(qc_dt, file.path(out_dir, "Table_Fig1_QC_Merged.csv"), sep = ",")
fwrite(site_long, file.path(out_dir, "Table_Fig1_8mod_SiteCounts_Long.csv"), sep = ",")
fwrite(as.data.table(cor_samples, keep.rownames = "Sample"),
       file.path(out_dir, "Table_Fig1_SampleCor.csv"), sep = ",")
fwrite(jac_dt, file.path(out_dir, "Table_Fig1_Jaccard_Stability.csv"), sep = ",")
fwrite(cv_dt[, .(Mod, CV_pct, mod_type)],
       file.path(out_dir, "Table_Fig1_8mod_CV.csv"), sep = ",")

# Save NCs3 deviation stats
cat("\n=== NCs3 Outlier Check ===\n")
ncs3_col <- which(colnames(site_counts) == "NCs3")
if (length(ncs3_col) > 0) {
  for (m in mods_8) {
    nc_vals <- site_counts[m, grep("^NC", colnames(site_counts))]
    if (all(is.na(nc_vals))) next
    nc_mean <- mean(nc_vals[-which(names(nc_vals) == "NCs3")], na.rm = TRUE)
    ncs3_val <- site_counts[m, "NCs3"]
    z <- (ncs3_val - nc_mean) / sd(nc_vals, na.rm = TRUE)
    if (!is.na(z) && abs(z) > 1.5) {
      cat(sprintf("  NCs3 outlier: %s Z=%.2f (NCs3=%d, other NC mean=%.0f)\n",
                  m, z, ncs3_val, nc_mean))
    }
  }
}

cat("\n=== Fig1 Panels C/E/F Complete ===\n")
cat("Output:", out_dir, "\n")
cat("  Fig1C_N50_Map_ReadLen.{pdf,png}\n")
cat("  Fig1E_Gene_PCA.{pdf,png}\n")
cat("  Fig1E_SampleCor_Heatmap.{pdf,png}\n")
cat("  Fig1F_Jaccard_Stability.{pdf,png}\n")
cat("Tables: Table_Fig1_*.csv (5 files)\n")
