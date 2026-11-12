# =============================================================================
# Fig 3F extension: SIX-timepoint Visium projection (Sham/4h/12h/2d/14d/42d)
# 2026-09-28 — direct module-score caliber (self-contained per section):
#   7 validated direction sets scored on every spot (AddModuleScore),
#   PT compartment = CT_PT score >= per-section upper quartile
#   (CT markers identical to canonical rebuild_fig3d_validated.R).
# Anchor check: Sham/12h/14d PT m6A_UP trend vs canonical scRNA-weight
#   +0.445 / +1.156 / -0.060 (reported, not forced).
# Runs on server (Seurat present). Output: GSE269622/sixtp_projection/
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})

setwd("/home/pjs/GSE269622")
out_dir <- "sixtp_projection"; dir.create(out_dir, showWarnings = FALSE)

cond_map <- list(Sham = "mouse_sham", Hour4 = "mouse_hour4", Hour12 = "mouse_hour12",
                 Day2 = "mouse_day2", Day14 = "mouse_day14", Week6 = "mouse_week6")
tp_levels <- c("Sham","Hour4","Hour12","Day2","Day14","Week6")

ct_markers <- list(
  PT = c("Slc34a1","Lrp2","Slc5a2","Slc22a6","Ggt1","Slc13a1","Slc22a8","Cubn","Miox"),
  TAL = c("Slc12a1","Umod","Kcnj1","Cldn16","Slc8a1"),
  DCT = c("Slc12a3","Pvalb","Trpm6"),
  CNT_CD = c("Calb1","Aqp2","Hsd11b2"),
  ICB = c("Slc4a1","Slc26a4","Atp6v1b1"),
  Podo = c("Nphs1","Nphs2","Podxl","Wt1","Synpo"),
  Endo = c("Pecam1","Cdh5","Eng","Kdr","Flt1"),
  Macrophage = c("Cd68","Csf1r","Adgre1","Itgam"),
  T_NK = c("Cd3e","Cd3d","Cd8a","Nkg7"),
  B_cell = c("Cd79a","Cd79b","Ms4a1"),
  Neutrophil = c("S100a8","S100a9","Retnlg"),
  Proliferating = c("Mki67","Top2a","Pcna"))

gs_dir <- "/home/pjs/scRNA_analysis/gene_sets_5mod"
mod_files <- c(m6A_UP = "M6A_UP.txt", m6A_DOWN_resid = "M6A_DOWN.txt",
               m5C_DOWN = "M5C_DOWN.txt", PSU_UP = "PSU_UP.txt",
               Inosine_UP = "INOSINE_UP.txt", Nm_UP = "NM_UP.txt", Nm_DOWN = "NM_DOWN.txt")
mod_sets <- lapply(mod_files, function(f)
  trimws(readLines(file.path(gs_dir, f))))

res_list <- list()
for (cond in tp_levels) {
  cat(sprintf("\n=== %s (%s) ===\n", cond, cond_map[[cond]]))
  v <- Load10X_Spatial(file.path(cond_map[[cond]], "outs"), slice = cond)
  v <- subset(v, subset = nFeature_Spatial > 200)          # canonical QC
  v <- NormalizeData(v, assay = "Spatial", verbose = FALSE)
  # CT scores (canonical markers)
  for (ct in names(ct_markers)) {
    genes <- intersect(ct_markers[[ct]], rownames(v))
    if (length(genes) < 2) next
    v <- AddModuleScore(v, features = list(genes), name = paste0("CT_", ct, "_x"),
                        ctrl = min(20, length(genes)))
    colnames(v@meta.data)[colnames(v@meta.data) == paste0("CT_", ct, "_x1")] <- paste0("CT_", ct)
  }
  # mod scores
  for (m in names(mod_sets)) {
    genes <- intersect(mod_sets[[m]], rownames(v))
    if (length(genes) < 5) { cat(sprintf("  %s: only %d genes\n", m, length(genes))); next }
    v <- AddModuleScore(v, features = list(genes), name = paste0(m, "_x"),
                        ctrl = min(50, length(genes)))
    colnames(v@meta.data)[colnames(v@meta.data) == paste0(m, "_x1")] <- m
  }
  md <- v@meta.data
  pt_thr <- quantile(md$CT_PT, 0.75, na.rm = TRUE)
  is_pt <- md$CT_PT >= pt_thr
  row <- data.table(timepoint = cond, n_spots = ncol(v), n_PT = sum(is_pt))
  for (m in names(mod_sets)) {
    if (!m %in% colnames(md)) next
    z <- as.numeric(scale(md[[m]]))                        # per-section pooled Z
    row[, (paste0(m, "_PT_Z")) := mean(z[is_pt])]
    row[, (paste0(m, "_all_Z")) := mean(z)]
    # also dominant-CT call per spot for context
  }
  # full CT x mod matrix (for heatmap extension of Fig 3F)
  for (ct in names(ct_markers)) {
    ccol <- paste0("CT_", ct)
    if (!ccol %in% colnames(md)) next
    ct_thr <- quantile(md[[ccol]], 0.75, na.rm = TRUE)
    is_ct <- md[[ccol]] >= ct_thr
    for (m in names(mod_sets)) {
      if (!m %in% colnames(md)) next
      z <- as.numeric(scale(md[[m]]))
      res_list[[paste(cond, ct, m)]] <- data.table(
        timepoint = cond, celltype = ct, mod = m,
        n_spots_ct = sum(is_ct), meanZ = mean(z[is_ct]))
    }
  }
  res_list[[paste(cond, "PTsummary")]] <- row
  print(as.data.frame(row))
  rm(v); gc(verbose = FALSE)
}

ctmod <- rbindlist(res_list[setdiff(names(res_list), paste(tp_levels, "PTsummary"))], fill = TRUE)
ptsum <- rbindlist(res_list[paste(tp_levels, "PTsummary")], fill = TRUE)
fwrite(ctmod, file.path(out_dir, "visium_6tp_ctmod_z.csv"))
fwrite(ptsum, file.path(out_dir, "visium_6tp_PT_summary.csv"))

cat("\n=== PT-compartment trajectory (direct module-score caliber) ===\n")
pz <- dcast(ptsum, timepoint ~ paste0("PT"), value.var = "m6A_UP_PT_Z")
print(ptsum[, .(timepoint, n_spots, n_PT, m6A_UP_PT_Z, m6A_DOWN_resid_PT_Z,
                m5C_DOWN_PT_Z, Nm_UP_PT_Z, Inosine_UP_PT_Z, Nm_DOWN_PT_Z, PSU_UP_PT_Z)])
cat("\nAnchor check (canonical scRNA-weight): Sham +0.445 / Hour12 +1.156 / Day14 -0.060\n")
cat("=== DONE ===\n")
