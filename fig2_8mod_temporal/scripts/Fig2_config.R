# Figure 2: 8-Modification Temporal Dynamics
# ============================================================
# ⚠️ LIVE INFRASTRUCTURE (2026-07-24注): 本文件被5个正典panel脚本source()
#    (Fig2C/G/I/J_v6/K)——头部注释是v3.0历史遗留，但palette/helper为现行依赖。
#    清理时勿归档，归档将连锁破坏5个panel。
# Narrative v3.0 (2026-06-27): L1 Chemical Stoichiometry + L2 Cross-Modification Network
# DRS USP: 1 DRS run = 8 traditional experiments (batch-free)
# "双向修饰波" framework retired — replaced by distribution + correlation

source("../../shared/palette/color_palette.R")

# ============================================================
# Panel Specifications v3.0 — 10 panels (A-L, skipping C and F)
# ============================================================
# L1 — 化学计量分布动力学
#   Panel A: 8-mod × 5-group Z-score panorama heatmap
#   Panel B: |Δ| distribution ridge plot (8 mods, absolute effect only)
#   Panel D: UMAP CellType atlas reference (scRNA-seq)

# L2 — 跨修饰协变网络
#   Panel E: 8×8 partial correlation matrix × 5 groups (TPM+nsites controlled)
#   Panel G: Enzyme timecourse + 28-enzyme Day1 FC (merged old E+F)

# 支撑 + 过渡
#   Panel H: Writer Heatmap — enzyme × cell type (old G→H)
#   Panel I: 9-Quadrant — modification × expression (old I→I, keep)
#   Panel J: Cell-type enrichment (old K→J)
#   Panel K: Composition Drift Exclusion (old L→K)
#   Panel L: Co-targeting genes → L3 bridge (old M+F genes → L)

# ============================================================
# Deleted / Moved
# ============================================================
# Old B (Suppressive DoubleBar)  → deleted
# Old C (Activating DoubleBar)   → deleted
# Old J (InterMod 4-in-1)        → correlation→E, Jaccard→Suppl
# Old F (28-Enzyme Day1 FC)      → merged into G
# Coverage stratification         → Supplementary Fig Sx
# Network rewiring diagram        → Supplementary Fig Sx

# ============================================================
# Panel Letter Mapping (Old → New)
# ============================================================
# A→A  B(ridge)→B  (old B/C deleted)  D→D
# E(Enzyme)→G  (old F→G)  G(Writer)→H  H(9Quad)→I
# K(CellType)→J  L(Drift)→K  M(Star)→L
# J(4in1) deleted, E(corr)→E(new)

# ============================================================
# Color Palette (shared across all panels)
# ============================================================
# 5-Mod Architecture (v7.0, 2026-07-08): Nm = Am+Cm+Gm+Um merged
MOD_COLORS <- c(
  "m6A"     = "#D6404C",  # red
  "m5C"     = "#E8753A",  # orange
  "PSU"     = "#3D7EA6",  # blue
  "Inosine"  = "#7B5294",  # purple
  "Nm"       = "#2E7D6F"   # teal-green
)
# Backward-compat aliases for old scripts
MOD_COLORS_ALIAS <- c(
  "m6a"="m6A","m5c"="m5C","psu"="PSU","inosine"="Inosine",
  "Gm"="Nm","Um"="Nm","Am"="Nm","Cm"="Nm"
)

MOD_ORDER <- c("m6A","m5C","PSU","Inosine","Nm")

# Derived palettes — all panels source from here for consistency
GRADIENT_LOW  <- "#3B4F8A"  # Gm deep blue → decrease/depletion end
GRADIENT_HIGH <- "#CC5A4C"  # m6A terracotta → increase/enrichment end
GRADIENT_MID  <- "white"

# Cell class colors (D, J, K panels)
CELL_CLASS_COLORS <- c(
  Epithelial    = "#D4A08A",  # muted warm
  Immune        = "#8CB5C8",  # soft blue-grey
  Stromal       = "#7FB5A0",  # muted sage
  Proliferating = "#B09CC0"   # soft lavender
)

# Enzyme function colors (K panel)
ENZYME_FUNC_COLORS <- c(
  Writer = "#6B7FAE",   # muted blue
  Eraser = "#C27A6A",   # muted terracotta
  Reader = "#A09580"    # warm grey
)

# Enzyme category colors (G panel) — harmonized with MOD_COLORS
ENZYME_CAT_COLORS <- c(
  "m6A Writers"     = "#B05A4C",  # darker terracotta (m6A family)
  "m6A Readers"     = "#D49585",  # light terracotta
  "m6A Erasers"     = "#E8B0A5",  # lighter terracotta
  "m5C Writers"     = "#8B9DAF",  # muted blue-grey
  "Nm Writers"      = "#6B5C8A",  # muted purple (Nm family, =Am+Cm+Gm+Um)
  "PSU Synthases"   = "#5C9E8A"   # muted green (PSU family)
)

# 9-Quadrant colors (I panel) — 2026-09-23 unified to project canonical
QUAD_COLORS <- c(
  Concordant = "#D6404C",  # project red (wet-lab Day1 canonical)
  Discordant = "#8A8A8A"   # neutral grey (IS canonical)
)

# ============================================================
# Data Paths
# ============================================================
BASE_METHY <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
BASE_EXPR  <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/03.expression"

# ============================================================
# Output
# ============================================================
OUTDIR <- "../outputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
