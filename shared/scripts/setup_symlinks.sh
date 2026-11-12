#!/bin/bash
# DRS Tech Showcase — Shared Data Symlinks
# ============================================================
# Creates symlinks from main paper data/cache to tech showcase shared/data/

SHARED_DIR="/mnt/d/desktop/DRS/DRS_TechShowcase/shared/data"
MAIN_DIR="/mnt/d/desktop/DRS/正刊-DRS-1/new"
E_DRS="/mnt/e/BN2505012AH02S66N1/DirectRNA"

echo "=== Creating shared data symlinks ==="
echo "Target: $SHARED_DIR"
mkdir -p "$SHARED_DIR"

# --- QC Data ---
echo "[1/8] QC tables..."
ln -sf "$MAIN_DIR/fig1/TableS_Fig1_QC_Summary.csv" "$SHARED_DIR/QC_Summary.csv" 2>/dev/null
ln -sf "$MAIN_DIR/supplement/panels/TableS_RNA_Integrity_QC.csv" "$SHARED_DIR/RNA_Integrity_QC.csv" 2>/dev/null

# --- Expression Data ---
echo "[2/8] Expression matrices..."
ln -sf "$E_DRS/result/03.expression/gene/gene.TPM.matrix.xls" "$SHARED_DIR/DRS_gene_TPM.xls" 2>/dev/null
ln -sf "$E_DRS/result/03.expression/transcript/isoforms.TPM.matrix.xls" "$SHARED_DIR/DRS_transcript_TPM.xls" 2>/dev/null

# --- Modification Cache (Fig3/Fig4 from main paper) ---
echo "[3/8] Spatial caches..."
ln -sf "$MAIN_DIR/fig4/cross_group_spatial_cache.rds" "$SHARED_DIR/cross_group_spatial_cache.rds" 2>/dev/null
ln -sf "$MAIN_DIR/fig4/spatial_2ome_cache.rds" "$SHARED_DIR/spatial_2ome_cache.rds" 2>/dev/null
ln -sf "$MAIN_DIR/fig4/frac_2ome_cache.rds" "$SHARED_DIR/frac_2ome_cache.rds" 2>/dev/null
ln -sf "$MAIN_DIR/fig4/chr_group_cache.rds" "$SHARED_DIR/chr_group_cache.rds" 2>/dev/null

# --- CDS-m6A Analysis ---
echo "[4/8] CDS-m6A data..."
ln -sf "$MAIN_DIR/fig3/CDS_m6A_analysis/cds_delta_annotated.csv" "$SHARED_DIR/cds_delta_annotated.csv" 2>/dev/null
ln -sf "$MAIN_DIR/fig3/CDS_m6A_analysis/complete_loss_annotated.csv" "$SHARED_DIR/complete_loss_annotated.csv" 2>/dev/null
ln -sf "$MAIN_DIR/fig3/CDS_m6A_analysis/m6a_sites_by_region_NC_vs_Day1.csv" "$SHARED_DIR/m6a_sites_by_region.csv" 2>/dev/null

# --- ATAC Data ---
echo "[5/8] ATAC data..."
ln -sf "$MAIN_DIR/fig6/data/differential_accessibility.csv" "$SHARED_DIR/ATAC_diff_accessibility.csv" 2>/dev/null
ln -sf "$MAIN_DIR/fig6/data/genome_wide_5kb_bins_sample5000.csv" "$SHARED_DIR/ATAC_5kb_bins.csv" 2>/dev/null
ln -sf "$MAIN_DIR/fig6/data/promoter_accessibility_by_gene.csv" "$SHARED_DIR/ATAC_promoter_by_gene.csv" 2>/dev/null
ln -sf "$MAIN_DIR/fig6/data/DRS_mod_ATAC_colocalization_sample5000.csv" "$SHARED_DIR/DRS_ATAC_colocalization.csv" 2>/dev/null
ln -sf "$MAIN_DIR/atac_results/" "$SHARED_DIR/atac_results" 2>/dev/null

# --- AS Data ---
echo "[6/8] AS data..."
ln -sf "$MAIN_DIR/fig2/data/all_sample_AS.xls" "$SHARED_DIR/all_sample_AS.xls" 2>/dev/null
ln -sf "$MAIN_DIR/fig2/data/novel_trans_code_stat.xls" "$SHARED_DIR/novel_trans_code_stat.xls" 2>/dev/null
ln -sf "$MAIN_DIR/fig2/data/number_stats.xls" "$SHARED_DIR/number_stats.xls" 2>/dev/null

# --- PolyA Data ---
echo "[7/8] PolyA data..."
ln -sf "$MAIN_DIR/fig2/data/polyA_stat.xls" "$SHARED_DIR/polyA_stat.xls" 2>/dev/null
ln -sf "$MAIN_DIR/fig5/data/new_transcripts.all_annotation.xls" "$SHARED_DIR/new_transcripts_annotation.xls" 2>/dev/null

# --- Novel Transcripts ---
echo "[8/8] Novel transcript data..."
ln -sf "$E_DRS/result/07.Structure_Analysis/02.novel_trans/new_transcripts.gtf" "$SHARED_DIR/new_transcripts.gtf" 2>/dev/null
ln -sf "$E_DRS/result/07.Structure_Analysis/02.novel_trans/new_transcripts.fa" "$SHARED_DIR/new_transcripts.fa" 2>/dev/null

echo ""
echo "=== Symlink creation complete ==="
echo "Note: Some symlinks may show errors if E: drive is not mounted."
echo "Run: ls -la $SHARED_DIR/ to verify."
ls -la "$SHARED_DIR/" 2>/dev/null | head -30
