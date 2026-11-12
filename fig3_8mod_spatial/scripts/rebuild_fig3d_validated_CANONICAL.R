# Rebuild Fig3D spatial projection using REAL per-timepoint Visium tissue
# (replaces single-Sham-reused-3x approach flagged as misleading)
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")  # fix: default C locale mangles unicode (x/em-dash) into "..." in output
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork); library(scales)
})

setwd("/home/pjs/GSE269622")
out_dir <- "real_5mod_validated"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load 3 REAL Visium tissue objects (Sham/Hour12/Day14, all GSE269622 CytAssist, same batch) ----
cond_map <- list(Sham = "mouse_sham", Hour12 = "mouse_hour12", Day14 = "mouse_day14")
vis_list <- list()
for (cond in names(cond_map)) {
  cat(sprintf("\n=== Loading %s ===\n", cond))
  v <- Load10X_Spatial(file.path(cond_map[[cond]], "outs"), slice = cond)
  cat(sprintf("  Raw: %d spotsFeature_Spatial > 200)
  v <- NormalizeData(v, assay = "Spatial", verbose = FALSE)
  cat(sprintf("  After QC: %d spots\n", ncol(v)))
  vis_list[[cond]] <- v
}

# ---- 2. Cell-type marker gene sets (same 12 types as original spatial_projection_5mod.R) ----
ct_markers <- list(
  PT            = c("Slc34a1","Lrp2","Slc5a2","Slc22a6","Ggt1","Slc13a1","Slc22a8","Cubn","Miox"),
  TAL           = c("Slc12a1","Umod","Kcnj1","Cldn16","Slc8a1"),
  DCT           = c("Slc12a3","Pvalb","Trpm6"),
  CNT_CD        = c("Calb1","Aqp2","Hsd11b2"),
  ICB           = c("Slc4a1","Slc26a4","Atp6v1b1"),
  Podo          = c("Nphs1","Nphs2","Podxl","Wt1","Synpo"),
  Endo          = c("Pecam1","Cdh5","Eng","Kdr","Flt1"),
  Macrophage    = c("Cd68","Csf1r","Adgre1","Itgam"),
  T_NK          = c("Cd3e","Cd3d","Cd8a","Nkg7"),
  B_cell        = c("Cd79a","Cd79b","Ms4a1"),
  Neutrophil    = c("S100a8","S100a9","Retnlg"),
  Proliferating = c("Mki67","Top2a","Pcna")
)

coverage <- rbindlist(lap- intersect(genes, rownames(vis_list[[1]]))
  data.table(CellType = ct, N_markers = length(genes), N_found = length(found),
             Found = paste(found, collapse = ","))
}))
cat("\n=== Marker gene coverage in REAL Visium (whole-transcriptome, should be much better than Xenium 300-panel) ===\n")
print(coverage)
fwrite(coverage, file.path(out_dir, "marker_coverage_check.csv"))

ct_scores <- list()
for (cond in names(vis_list)) {
  v <- vis_list[[cond]]
  for (ct in names(ct_markers)) {
    genes <- intersect(ct_markers[[ct]], rownames(v))
    if (length(genes) < 2) { cat(sprintf("  [%s] %s: too few markers (%d), skipping\n", cond, ct, length(genes))); next }
    score_name <- paste0("CT_", ct, "_tmp")
    v <- AddModuleScore(v, features = list(genes), name = score_name, ctrl = min(20, length(genes)))
    colnames(v@meta.data)[colnames(v@meta.data) == paste0(score_name, "1")] <- paste0("CT_", ct)
  }
  vis_list[[cond]] <- v
}

# ---- 3. Load REAL scRNA-deriveal/5mod_diff/mod_celltype_timepoint_zscore_validated.csv")
MODS <- c("m6A_UP","m6A_DOWN_resid","m5C_DOWN","PSU_UP","Inosine_UP","Nm_UP","Nm_DOWN")

# Real Visium timepoint -> closest scRNA z-score bucket (documented mapping, NOT claimed identical)
tp_to_zbucket <- list(Sham = "Sham", Hour12 = "Early", Day14 = "Late")
tp_display <- list(Sham = "Sham (real tissue)", Hour12 = "12h (real tissue) × scRNA-Early weights",
                    Day14 = "14d (real tissue) × scRNA-Late weights")

theme_spatial <- theme(
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "grey40"),
  legend.position = "bottom", legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.6, "cm"),
  legend.text = element_text(size = 6), legend.title = element_text(size = 7)
)
spatial_palette <- c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")
esid = expression(bold("m"^"6"*"A"*"↓ (residual)")), m5C_DOWN = expression(bold("m"^"5"*"C"*"↓")), PSU_UP = expression(bold("PSU ↑")), Inosine_UP = expression(bold("Inosine ↑")), Nm_UP = expression(bold("Nm ↑ (2′-O-Me)")), Nm_DOWN = expression(bold("Nm ↓ (2′-O-Me)")))

plot_list <- list()
summary_rows <- list()

for (mod in MODS) {
  for (cond in names(vis_list)) {
    v <- vis_list[[cond]]
    zbucket <- tp_to_zbucket[[cond]]
    proj_name <- paste0(mod, "_", cond, "_Proj")
    v@meta.data[[proj_name]] <- 0
    ct_used <- 0
    for (ct in names(ct_markers)) {
      ct_col <- paste0("CT_", ct)
      if (!ct_col %in% colnames(v@meta.data)) next
      z_row <- mod_z[CellType == ct & Timepoint == zbucket]
      if (nrow(z_row) != 1 || !(mod %in% colnames(z_row))) next
      z_val <- as.numeric(z_row[[mod]])
      if (is.na(z_val)) next
      ct_vals <- v@meta.data[[ct_col]]
      ct_vals_scaled <- (ct_vals - mean(ct_vals, na.rm = TRUE)) / sd(ct_val + z_val * ct_vals_scaled
      ct_used <- ct_used + 1
    }
    vis_list[[cond]] <- v
    summary_rows[[paste(mod, cond)]] <- data.table(Mod = mod, Timepoint = cond, N_CellTypes_Used = ct_used,
      Min = min(v@meta.data[[proj_name]], na.rm=TRUE), Max = max(v@meta.data[[proj_name]], na.rm=TRUE))

    p <- SpatialFeaturePlot(v, features = proj_name, pt.size.factor = 1.8, alpha = c(0.05, 1), image.alpha = 0.7) +
      scale_fill_gradientn(colors = spatial_palette, name = "Score") +
      labs(title = mod_labels[[mod]], subtitle = tp_display[[cond]]) + theme_spatial
    key <- paste(mod, cond, sep = "_")
    plot_list[[key]] <- p
    ggsave(file.path(out_dir, sprintf("Visium_%s_%s_VALIDATED.pdf", mod, cond)), p, width = 4.0, height = 3.8, device = "pdf")
    ggsave(file.path(out_dir, sprintf("Visium_%s_%s_VALIDATED.png", mod, cond)), p, width = 4.0, height = 3.8, dpi = 400)
    cat(sprintf("  Done: %s @ %s (%d cell types)\n", mod, cond, ct_used))
  }
}

fwrite(r-- 4. Combined figure ----
tp_labels <- c(Sham = "Sham", Hour12 = "12h (real)", Day14 = "14 day (real)")
row_list <- list(); panel_idx <- 0
for (mod in MODS) {
  mod_plots <- list()
  for (cond in names(vis_list)) {
    panel_idx <- panel_idx + 1
    key <- paste(mod, cond, sep = "_")
    mod_plots[[cond]] <- plot_list[[key]] + labs(tag = LETTERS[panel_idx]) + theme(plot.tag = element_text(size = 10, face = "bold"))
  }
  row_list[[mod]] <- wrap_plots(mod_plots, ncol = 3)
}
combined <- wrap_plots(row_list, ncol = 1) +
  plot_annotation(
    title = "Spatial Dynamics of Validated RNA-Modification Direction Sets — REAL Independent Tissue per Timepoint",
    subtitle = paste0("Rows: validated direction sets | Columns: Sham / 12h / 14-day — three GENUINELY DIFFERENT tissue sections (GSE269622 Visium CytAssist)\n",
                       "Per-spot score = real cell-type composition (this tissue, this timepoint) × scRNA-derived (celltype×timepoint×mod) weight.'Late' bucket (closest available match, not identical timepoints)."),
    caption = "Visium: GSE269622 (Sham/12h/14d, real independent sections). scRNA weights: GSE274819 IRI time course. Validated direction-set weights (dual-anchor arbitration basis); supersedes legacy-set version.",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 7, hjust = 0.5, color = "grey30", margin = margin(t=4,b=8)),
                  plot.caption = element_text(size = 6, hjust = 0.5, color = "grey50", margin = margin(t=8))))
ggsave(file.path(out_dir, "Fig3D_Projection_5mod_3timepoints_VALIDATED.pdf"), combined, width = 13, height = 18, device = "pdf")
ggsave(file.path(out_dir, "Fig3D_Projection_5mod_3timepoints_VALIDATED.png"), combined, width = 13, height = 18, dpi = 400)

saveRDS(vis_list, file.path(out_dir, "visium_3timepoint_real_processed.rds"))
cat(sprintf("\n=== Done! %d panels, %d mods x %rey30", margin = margin(t=4,b=8)),
                  plot.caption = element_text(size = 6, hjust = 0.5, color = "grey50", margin = margin(t=8))))
ggsave(file.path(out_dir, "Fig3D_Projection_5mod_3timepoints_VALIDATED.pdf"), combined, width = 13, height = 18, device = "pdf")
ggsave(file.path(out_dir, "Fig3D_Projection_5mod_3timepoints_VALIDATED.png"), combined, width = 13, height = 18, dpi = 400)

saveRDS(vis_list, file.path(out_dir, "visium_3timepoint_real_processed.rds"))
cat(sprintf("\n=== Done! %d panels, %d mods x %d real timepoints ===\n", panel_idx, 