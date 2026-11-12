# Figure 1: DRS Platform Overview & Data Quality
# ============================================================
# Narrative: "DRS data is reliable" — Establish reader trust
# DRS USP: PCR-free, direct RNA seq, native strand-specificity

source("../../shared/palette/color_palette.R")

# ============================================================
# Panel Specifications
# ============================================================
# Panel A: DRS原理图 + 能力对比矩阵 (需Illustrator绘制)
# Panel B: 实验设计总览 + 数据产出汇总 (Fig1A雏形升级)
# Panel C: N50 violin + Mapping rate bar + Read length density
# Panel D: 5'-3' metagene coverage + strand specificity
# Panel E: PCA (gene + transcript) + per-sample site correlation heatmap
# Panel F: m6A per-sample site count + Jaccard (NC vs NC-Day1) + 8-mod CV%

# ============================================================
# Data Sources
# ============================================================
# QC Table:   fig1/data/TableS_Fig1_QC_Summary.csv (from main paper)
# Collapse stats: E:\...\02.alignment\collapse_stats.xls
# Align stats:    E:\...\02.alignment\align_stats.xls
# Expression:     E:\...\03.expression\gene\gene.TPM.matrix.xls
# 8-mod site counts: E:\...\05.methy\{mod}\methy_detect\all_sample.{mod}.0.7.xls

# ============================================================
# Existing Assets (from main paper)
# ============================================================
# - Fig1F_QC_overview.R (N50, mapping rate)
# - fig1_PCA_final.R (Gene + Transcript PCA)
# - analysis_n50_rna_integrity.R (N50 analysis)
# - analysis_coverage_uniformity.R (5'-3' coverage)
# - analysis_dorado_caller_doc.R (Dorado caller docs)
# - TableS_Fig1_QC_Summary.csv (15-sample QC table)

# ============================================================
# New Work Needed
# ============================================================
# 1. Panel A: Illustrator schematic (new)
# 2. Panel B: Upgrade Fig1A with data output summary table
# 3. Panel E: Generate 8-mod 15x15 site count correlation heatmap
# 4. Panel F: Extend Jaccard from m6A-only to all 8 modifications
# 5. FigS1C: 8-mod per-sample site count bar chart
# 6. FigS1D: DRS vs RNA-seq Bland-Altman + edge distribution

message("Fig1 config loaded. Panels: A(Illustrator) B(Design) C(N50+Map) D(Coverage) E(PCA+Heatmap) F(Jaccard)")
