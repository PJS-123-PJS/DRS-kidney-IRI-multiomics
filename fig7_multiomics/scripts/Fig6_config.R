# Figure 6: PT-S3 Multi-Omics Integration — DRS Single-Platform Advantage
# ============================================================================
# Narrative v8.0 (2026-07-09): Centered on PT-S3 — the most IRI-vulnerable segment.
# DRS captures expression + modifications + splicing + PolyA in ONE run.
# PT-S3 is where multi-omics converge: acute m6A/m5C suppression + delayed AS/PolyA.
#
# Architecture: 5 panels (A-E), PT-S3 focused
# ============================================================================

source("../../shared/palette/color_palette.R")

# ============================================================
# Panel Specifications v8.0 — PT-S3 Centered
# ============================================================
# Panel A: PT-S3 Multi-Omics Circos — 4 omics layers with PT gene annotation
# Panel B: Temporal Decoupling in PT — acute (m6A+Expr) vs delayed (AS+PolyA+ATAC)
# Panel C: PT-S3 Hub Gene Network — multi-layer regulatory hubs
# Panel D: SAM-m6A Metabolic Link — methionine cycle → PT m6A dynamics
# Panel E: DRS Single-Platform Advantage — layer coverage per gene, PT vs non-PT

# ============================================================
# PT-S3 Gene Set
# ============================================================
# PT segment markers (mouse kidney)
PT_S1 <- c("Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc22a8")           # early PT
PT_S2 <- c("Slc13a1","Slc22a7","Slc13a3","Slc6a19")                  # mid PT
PT_S3 <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1")  # late PT (IRI-vulnerable)
PT_ALL <- unique(c(PT_S1, PT_S2, PT_S3))

# ============================================================
# Data Paths
# ============================================================
DATA_MULTI   <- "../data/Table_MultiLayer_Genes_validated.csv"
DATA_ATAC    <- "../data/ATAC_DRS_merged.csv"
DATA_LAYERS  <- "../data/Table_Layer_Sizes.csv"
DATA_METAB   <- "../output/Fig7J_Sig_Metabolites_Full.csv"

# ============================================================
# Key Findings
# ============================================================
# - ATAC-DRS Expression r=0.008 (temporal decoupling confirmed)
# - DRS-RNAseq gene expression r=0.645 (cross-platform validation)
# - 20 hub genes in >=3 omics layers
# - Methionine -1.83, SAH -1.38 (SAM cycle depletion → m6A decrease)
# - PT-S3 is the convergence zone: m6A targets + AS/PolyA changes + metabolic stress

# ============================================================
# Output
# ============================================================
OUTDIR <- "../output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("Fig7 v8.0 config loaded — PT-S3 centered, 5-panel architecture")
