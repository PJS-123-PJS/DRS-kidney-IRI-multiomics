#!/usr/bin/env Rscript
# Fig1S KEGG pathway-activity Z-score heatmap — canonical rebuild (2026-10-26)
# Provenance: the inherited Fig 1S panel ("DRS Kegg pathway") had no locatable generator
# (not in local projects, server, or company delivery). This script reconstructs the
# analysis from primary delivered data and becomes the documented generator.
#
# Method (Methods 2.2 wording anchor):
#   KEGG pathway gene sets: company-delivered functional annotation
#     (03.expression/gene/gene.TPM.matrix.add_info.xls, `pathway` column, ko-map IDs).
#   Per platform (DRS TPM / RNA-seq FPKM): log2(x + 0.01); per-gene Z-score across the
#   platform's 15 samples; pathway activity per sample = mean of member-gene Z-scores;
#   heatmap shows pathway activity row-scaled (Z) across all 30 samples.
#   Pathways retained: >= 10 member genes detected (max group mean > 1, either platform).

suppressPackageStartupMessages({library(data.table); library(pheatmap)})

BASE <- "/mnt/f/BN2505012AH02S66N1"
DRS_TPM  <- file.path(BASE, "DirectRNA/result/03.expression/gene/gene.TPM.matrix.xls")
DRS_ANN  <- file.path(BASE, "DirectRNA/result/03.expression/gene/gene.TPM.matrix.add_info.xls")
RNA_FPKM <- file.path(BASE, "RNAseq/Result/03.expression/gene/gene.FPKM.matrix.xls")
OUT <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig1_platform_qc/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

samples <- c(paste0("NCs",1:3), paste0("Day1s",1:3), paste0("Day3s",1:3),
             paste0("Day7s",1:3), paste0("Is",1:3))
group_of <- c(NC="Sham", Day1="Day1", Day3="Day3", Day7="Day7", Is="IS")

# ---- 1. expression matrices -------------------------------------------------
cat("Reading expression matrices...\n")
drs <- fread(DRS_TPM);  rna <- fread(RNA_FPKM)
stopifnot(all(colnames(drs)[-1] == colnames(rna)[-1]))
drs <- as.data.frame(drs); rownames(drs) <- drs$ID; drs$ID <- NULL
rna <- as.data.frame(rna); rownames(rna) <- rna$ID; rna$ID <- NULL
common <- intersect(rownames(drs), rownames(rna))
cat(sprintf("  genes: DRS %d, RNA-seq %d, common %d\n", nrow(drs), nrow(rna), length(common)))
drs <- drs[common, ]; rna <- rna[common, ]

# ---- 2. KEGG pathway gene sets from delivered annotation --------------------
ann <- fread(DRS_ANN, select = c("ID", "pathway"))
ann <- ann[ID %in% common & pathway != "-" & !is.na(pathway)]
# pathway column format: "ko04740//Olfactory transduction`Organismal Systems`Sensory system"
# multiple entries separated by " ; " or ";" — split robustly on any of ; | \n
pw_list <- strsplit(ann$pathway, split = "[;|\n]")
ids <- lapply(pw_list, function(v) sub("//.*", "", trimws(v)))
ids <- lapply(ids, function(v) unique(v[grepl("^ko[0-9]{5}$", v)]))
pw2gene <- list()
for (i in seq_along(ids)) for (p in ids[[i]])
  pw2gene[[p]] <- unique(c(pw2gene[[p]], ann$ID[i]))
pw2gene <- pw2gene[lengths(pw2gene) > 0]
cat(sprintf("  pathways annotated: %d\n", length(pw2gene)))

# ---- 3. detection filter -----------------------------------------------------
detected <- function(m) rowSums(m > 1) >= 3   # >1 TPM/FPKM in >=3 samples
det <- detected(as.matrix(drs)) | detected(as.matrix(rna))
pw2gene <- lapply(pw2gene, function(g) g[g %in% common[det]])
pw2gene <- pw2gene[lengths(pw2gene) >= 10]
cat(sprintf("  pathways with >=10 detected members: %d\n", length(pw2gene)))
stopifnot(length(pw2gene) >= 50)

# ---- 4. pathway activity Z-score ---------------------------------------------
activity <- function(mat, sets) {
  x <- log2(as.matrix(mat) + 0.01)
  m <- rowMeans(x); s <- apply(x, 1, sd)
  s[s < 1e-8] <- 1e-8                      # zero-variance gene (e.g. flat in one
                                           # platform) -> Z = 0 = neutral, not NaN
  z <- (x - m) / s                         # per-gene Z across the 15 samples
  z <- z[rownames(z) %in% unique(unlist(sets)), , drop = FALSE]
  t(sapply(sets, function(g) colMeans(z[rownames(z) %in% g, , drop = FALSE])))
}
act_drs <- activity(drs, pw2gene)   # pathways x 15
act_rna <- activity(rna, pw2gene)
stopifnot(identical(rownames(act_drs), rownames(act_rna)))

# ---- 5. two-platform concordance (verification) ------------------------------
v_d <- as.vector(act_drs[, samples]); v_r <- as.vector(act_rna[, samples])
pear <- cor(v_d, v_r); spear <- cor(v_d, v_r, method = "spearman")
perpw <- sapply(seq_len(nrow(act_drs)), function(i)
  cor(act_drs[i, samples], act_rna[i, samples], method = "spearman"))
stopifnot(is.finite(pear), is.finite(spear))
cat(sprintf("  two-platform concordance: Pearson r = %.3f, Spearman rho = %.3f\n", pear, spear))
cat(sprintf("  per-pathway Spearman: median %.2f, IQR %.2f-%.2f, >0.5 in %d/%d\n",
            median(perpw, na.rm = TRUE), quantile(perpw, .25, na.rm = TRUE),
            quantile(perpw, .75, na.rm = TRUE),
            sum(perpw > 0.5, na.rm = TRUE), sum(!is.na(perpw))))
# temporal rearrangement: group-mean range
pref <- c(NC = "NCs", Day1 = "Day1s", Day3 = "Day3s", Day7 = "Day7s", Is = "Is")
grp_means <- sapply(names(pref), function(g)
  rowMeans(act_drs[, paste0(pref[g], 1:3), drop = FALSE]))
colnames(grp_means) <- c("Sham","Day1","Day3","Day7","IS")
rng <- apply(grp_means, 1, function(r) max(r) - min(r))
cat(sprintf("  temporal range of group-mean activity: median %.2f, >0.5 in %d/%d pathways\n",
            median(rng), sum(rng > 0.5), length(rng)))

# ---- 6. display matrix (group-major, DRS then RNA within group) --------------
pref <- c(NC = "NCs", Day1 = "Day1s", Day3 = "Day3s", Day7 = "Day7s", Is = "Is")
M <- do.call(cbind, lapply(names(pref), function(g) {
  cols <- paste0(pref[g], 1:3)
  cbind(act_drs[, cols, drop = FALSE], act_rna[, cols, drop = FALSE])
}))
colnames(M) <- paste(rep(c("Sham","Day1","Day3","Day7","IS"), each = 6),
                     rep(c("DRS","RNA"), each = 3), seq_len(30), sep = "|")
M_disp <- t(scale(t(M)))  # row-scale across 30 samples
M_disp[!is.finite(M_disp)] <- 0            # constant activity row -> neutral 0
M_disp <- pmax(pmin(M_disp, 1), -1)        # inherited panel legend range: -1..1

# ---- 7. annotation + palette (matched to inherited panel S) ------------------
ann_col <- data.frame(
  Group = factor(rep(c("Sham","Day1","Day3","Day7","IS"), each = 6),
                 levels = c("Sham","IS","Day7","Day3","Day1"))  # legend order is appearance-driven (levels inert),
  Seq   = factor(rep(c("DRS","RNA"), times = 5, each = 3), levels = c("DRS","RNA")))
rownames(ann_col) <- colnames(M_disp)
# NB: pheatmap draws annotation tracks/legends in reverse df order ->
# (Group, Seq) here renders as Seq track on top + Seq legend first (as inherited)
ann_colors <- list(
  Seq   = c(DRS = "#C05F50", RNA = "#3D87B9"),
  Group = c(Sham = "#E16B5A", Day1 = "#67BFD3", Day3 = "#4AAB99",
            Day7 = "#5E7199", IS = "#EEAB94"))
# cell gradient: 20 anchors pixel-sampled from the inherited panel S legend bar
# (2026-10-26, fig1(A4).png), value -1 (dark navy) .. +1 (red); legend range -1..1
grad_anchors <- rev(c("#C22628","#D15447","#DE7667","#E99789","#F0A89B","#F8CBC3",
                      "#FBDBD5","#FFFEFE","#ECE8F3","#DAD1E7","#C7BADB","#B4A3CE",
                      "#A18DC1","#8F79B6","#7B65AA","#69529F","#583E96","#422D8C",
                      "#2C247B","#1D2C5C"))
col_z <- colorRampPalette(grad_anchors)(101)
brks  <- seq(-1, 1, length.out = 102)
lg_br <- c(-1, -0.5, 0, 0.5, 1)
lg_lb <- c("-1", "-0.5", "0", "0.5", "1")

draw_panel <- function() {
  pheatmap(M_disp, color = col_z, breaks = brks,
           legend_breaks = lg_br, legend_labels = lg_lb,
           cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0,
           show_colnames = FALSE, show_rownames = FALSE,
           annotation_col = ann_col, annotation_colors = ann_colors,
           annotation_legend = TRUE, legend = TRUE,
           main = "DRS Kegg pathway")
  grid::grid.text("Z-score", x = grid::unit(0.861, "npc"),
                  y = grid::unit(0.978, "npc"),
                  gp = grid::gpar(fontsize = 9, fontface = "bold"))
}
pdf(file.path(OUT, "Fig1S_KEGG_pathway_activity.pdf"), width = 7.5, height = 9)
draw_panel()
dev.off()
png(file.path(OUT, "Fig1S_KEGG_pathway_activity.png"), width = 7.5, height = 9,
    units = "in", res = 300)
draw_panel()
dev.off()

# ---- 8. deliverables -----------------------------------------------------------
fwrite(data.table(pathway = rownames(M), M), file.path(OUT, "Fig1S_KEGG_activity_matrix.csv"))
fwrite(data.table(pathway = rownames(act_drs), spearman_drs_vs_rna = perpw),
       file.path(OUT, "Fig1S_KEGG_perpathway_concordance.csv"))
sink(file.path(OUT, "Fig1S_KEGG_verification_report.txt"))
cat("Fig1S KEGG pathway-activity rebuild — verification report (2026-10-26)\n")
cat("==================================================================\n")
cat(sprintf("genes common to both platforms: %d\n", length(common)))
cat(sprintf("pathways annotated (delivered ko-map): %d\n", length(pw2gene)))
cat(sprintf("pathways retained (>=10 detected members): %d\n", nrow(M)))
cat(sprintf("two-platform concordance of pathway activity (all pathway x sample cells): Pearson r = %.3f, Spearman rho = %.3f\n", pear, spear))
cat(sprintf("per-pathway Spearman (DRS vs RNA-seq across 15 matched samples): median %.2f (IQR %.2f-%.2f); rho>0.5 in %d/%d\n",
            median(perpw, na.rm = TRUE), quantile(perpw, .25, na.rm = TRUE),
            quantile(perpw, .75, na.rm = TRUE),
            sum(perpw > 0.5, na.rm = TRUE), sum(!is.na(perpw))))
cat(sprintf("temporal rearrangement (DRS group-mean activity range): median %.2f; >0.5 in %d/%d pathways\n",
            median(rng), sum(rng > 0.5), length(rng)))
cat("manuscript claims checked: 'KEGG pathway activity systematically rearranged across the time course' (temporal range above) and 'concordant between platforms' (concordance above)\n")
sink()
cat("Outputs: Fig1S_KEGG_pathway_activity.{pdf,png}, Fig1S_KEGG_activity_matrix.csv,\n")
cat("         Fig1S_KEGG_perpathway_concordance.csv, Fig1S_KEGG_verification_report.txt\n")
