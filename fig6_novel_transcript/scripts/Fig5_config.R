# Figure 5: Novel Transcript Discovery
# ============================================================
# Narrative: "DRS discovers new transcripts beyond reference annotation"
# DRS USP: Full-length reads enable novel isoform detection from one experiment
# Audit 2026-06-19: All panels verified with real data. rlnorm() removed.
#   Merged: A+C→B (classification panel). Deleted: Fig5L/M (orphans).
#   Renumbered: B-K→A-G (7 panels, 2026-06-23: I deleted (dup→F), H→Fig6I)

source("../../shared/palette/color_palette.R")

# ============================================================
# Panel Specifications (A-G, continuous)
# ============================================================
# Panel A: Chromosome distribution by class code (from structural CSV)
# Panel B: Class code pie + length + exon count combined (from structural CSV)
# Panel C: m6A modification density: Novel vs Known (from server CSV)
# Panel D: GO enrichment parsed from annotation XLS
# Panel E: Class code × cell type enrichment (server scRNA enrichment CSV)
# Panel F: Temporal expression NC→Day1 (from novel_transcript_expression.csv)
# Panel G: Expression comparison Novel vs Known (from expression CSVs)
# (Panel I deleted: duplicate of F; Panel H moved to Fig6I: multi-layer overlap)

# ============================================================
# Data Sources (verified 2026-06-19)
# ============================================================
# Core:    novel_transcript_structural.csv (19,487 transcripts, 7 class codes)
# GO:      new_transcripts.all_annotation.xls (14,381 transcripts, parsed GO)
# Density: Fig5E_real_m6A_density.csv → renamed to Fig5C_real_m6A_density.csv
# Enrich:  enrichment_class_code_zscore.csv (server: scRNA class code × cell type)
# Express: novel_transcript_expression.csv + known_mRNA_expression.csv
# GeneSets: geneset_class_{i,j,k,m,n,u,x}_symbols.txt (per-class-code host genes)
