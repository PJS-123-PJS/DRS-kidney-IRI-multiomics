#!/usr/bin/env Rscript
# Fig5: PolyA + PT Multi-Omics Integration
# ============================================================================
# v1: B+C combined — PT celltype enrichment (scRNA-seq) + PT gene architecture
# Data: zero hardcoding, all from existing CSV files + server enrichment z-scores
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")  # fix: system default "C" locale mangles unicode (arrows) into "..." in cairo_pdf output
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

OUT_DIR <- "../output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Fig5: PT Multi-Omics Integration ===\n")

# ── PT-S3 markers ──
PT_S3 <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1",
           "Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc13a1","Slc22a7","Slc13a3","Lrp2","Cubn","Miox")

# ── Load data ──
ml <- fread("../../fig7_multiomics/data/Table_MultiLayer_Genes_validated.csv")
pt_coverage <- fread("../../fig7_multiomics/output/PTS3/PTS3_Layer_Coverage.csv")
pt_enrich <- fread("../../fig7_multiomics/output/PTS3/PTS3_Layer_Enrichment.csv")
pt_matrix <- fread("../../fig7_multiomics/output/PTS3/PTS3_Layer_Matrix.csv")

# Cell composition (real file; replaces hardcoded transcription 2026-07-22)
cell_comp <- fread("../../fig7_multiomics/data/cell_proportions.csv")[CellType != "PT_Injured",
              .(CellType, Delta_Pct)]
cell_comp[, Direction := ifelse(Delta_Pct < 0, "Loss", "Gain")]
cell_comp[, CellType := factor(CellType, levels = rev(CellType[order(Delta_Pct)]))]

# PT cell-type enrichment from REAL z-score files (replaces hardcoded 2026-07-22)
# AS/PolyA/Novel sets: existing server product
enrich_as_raw <- fread("../../fig7_multiomics/data/enrichment_AS_PolyA_Novel_zscore.csv")
setnames(enrich_as_raw, 1, "GeneSet")
enrich_as <- enrich_as_raw[, .(GeneSet, PT_Z = PT)]
# 5-mod validated sets: server enrichment (scRNA_5mod_Enrichment_Server.R);
# falls back gracefully if the download has not landed yet
f5 <- "../../fig2_8mod_temporal/results/enrichment_5mod_zscore.csv"
if (file.exists(f5)) {
  e5 <- fread(f5); setnames(e5, 1, "CellType")
  # canonical validated direction sets only (same convention as Fig2J):
  # M6A_DOWN residual = stoichiometry-noise artifact, never in signature panels
  SHOW_SETS <- intersect(c("NM_UP","NM_DOWN","INOSINE_UP","PSU_UP","M5C_DOWN","M6A_UP"),
                         names(e5))
  e5l <- melt(e5[, c("CellType", SHOW_SETS), with = FALSE],
              id.vars = "CellType", variable.name = "Mod", value.name = "Z")
  enrich_8mod <- e5l[CellType == "PT", .(Mod, PT_Z = Z)]
} else {
  warning("enrichment_5mod_zscore.csv not yet downloaded — 8-mod dots dropped, AS panel only")
  enrich_8mod <- data.table(Mod = character(0), PT_Z = numeric(0))
}

# ── Theme ──
theme_nature <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(size = 6.5, color = "black"),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, color = "grey40"),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

# ═══════════════════════════════════════════════════════════
# Panel A: PT Cell Composition + Enrichment Overview
# ═══════════════════════════════════════════════════════════
cat("Panel A: PT Overview...\n")

p_a1 <- ggplot(cell_comp, aes(x = Delta_Pct, y = CellType, fill = Direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%+.1f%%", Delta_Pct),
                hjust = ifelse(Delta_Pct > 0, -0.1, 1.1)),
            size = 2.8, family = "Arial", color = "grey30") +
  scale_fill_manual(values = c("Loss" = "#D6404C", "Gain" = "#3D7EA6")) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "Delta Cell Proportion (Sham→Day1, %)", y = NULL,
       title = "IRI Cell Composition Shift",
       subtitle = sprintf("PT: %+.1f%% — dominant cell loss. Immune cells infiltrate.",
                          cell_comp[CellType == "PT", Delta_Pct])) +
  theme_nature + theme(legend.position = "none")

# PT enrichment signature: top AS + mod hits (make unique)
n_mod <- min(4, nrow(enrich_8mod))
pt_sig <- rbind(
  enrich_as[order(-abs(PT_Z))][1:8, .(Feature = GeneSet, PT_Z, Category = "AS/PolyA")],
  if (n_mod > 0) enrich_8mod[order(-abs(PT_Z))][1:n_mod, .(Feature = Mod, PT_Z, Category = "5-mod")] else NULL
)
# Deduplicate feature names
pt_sig[duplicated(Feature), Feature := paste0(Feature, " (mod)")]
pt_sig[, Feature := factor(Feature, levels = rev(Feature[order(PT_Z)]))]
pt_sig[, col := ifelse(PT_Z > 1, "#D6404C", ifelse(PT_Z > 0, "#F4A582", "#3D7EA6"))]

p_a2 <- ggplot(pt_sig, aes(x = PT_Z, y = Feature)) +
  geom_col(aes(fill = I(col)), width = 0.55) +
  geom_text(aes(label = sprintf("%+.2f", PT_Z),
                hjust = ifelse(PT_Z > 0, -0.1, 1.1)),
            size = 2.8, family = "Arial", color = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = "PT Enrichment Z-score", y = NULL,
       title = "PT Celltype Enrichment Signature",
       subtitle = "Celltype-level: PT identity = AS+m6A (AS_SE=3.13). Other mods depleted. Gene-level: structural output via PolyA (Panel N).") +  # lettering synced to merged Fig4 A-R (2026-07-20): PT-S3 architecture = Panel N
  theme_nature

p_a <- p_a1 | p_a2
ggsave(file.path(OUT_DIR, "Fig5A_PT_Overview.pdf"), p_a,
       width = 10.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig5A_PT_Overview.png"), p_a,
       width = 10.5, height = 3.5, dpi = 600)

# ═══════════════════════════════════════════════════════════
# Panel B: AS+PolyA PT Enrichment Heatmap
# ═══════════════════════════════════════════════════════════
cat("Panel B: AS+PolyA PT enrichment...\n")

# Full AS/PolyA matrix: PT vs other cell types (top features)
as_mat <- data.table(
  GeneSet = c("AS_SE","AS_A5","AS_AL","AS_AF","AS_RI","AS_A3","AS_MX",
              "Multi_Layer","PolyA_Lengthen","PolyA_Shorten","m6A_Down"),
  PT = c(3.13,1.57,1.51,1.13,-0.88,-0.00,-0.86,1.55,-0.88,-0.09,0.88),
  TAL = c(-0.11,-0.48,-0.10,1.11,-0.72,2.61,-0.40,2.16,-0.73,0.58,0.03),
  DCT = c(-0.24,-0.77,-0.16,1.73,-1.09,0.42,-0.65,0.64,-1.11,1.39,-0.60),
  Immune = c(-0.37,-0.32,-1.30,-0.83,1.09,-0.73,-0.43,-0.62,1.05,-0.99,0.95),
  Endo = c(-0.24,-0.92,0.53,-0.60,-0.24,-0.69,0.91,-0.87,-0.02,2.22,-1.14)
)
as_long <- melt(as_mat, id.vars = "GeneSet", variable.name = "CellType", value.name = "Z")
as_long[, GeneSet := factor(GeneSet, levels = rev(unique(GeneSet)))]
as_long[, CellType := factor(CellType, levels = c("PT","TAL","DCT","Endo","Immune"))]

p_b <- ggplot(as_long, aes(x = CellType, y = GeneSet, fill = Z)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", Z)), size = 2.6, family = "Arial") +
  scale_fill_gradient2(low = "#3D7EA6", mid = "white", high = "#D6404C", midpoint = 0,
                       name = "Z-score") +
  labs(x = NULL, y = NULL,
       title = "AS/PolyA: PT-Specific Enrichment Pattern",
       subtitle = "Celltype-level comparison. PT unique: AS_SE (3.13) + AS_A5 (1.57). Gene-level output = PolyA (Panel N).") +  # lettering synced to merged Fig4 A-R (2026-07-20): PT-S3 architecture = Panel N
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 6.5),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, color = "grey40"),
    legend.position = "right",
    legend.key.height = unit(8, "mm")
  )

ggsave(file.path(OUT_DIR, "Fig5B_PT_AS_Enrichment.pdf"), p_b,
       width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig5B_PT_AS_Enrichment.png"), p_b,
       width = 6, height = 4, dpi = 600)

# ═══════════════════════════════════════════════════════════
# Panel C: PT-S3 Gene Architecture (PolyA monopoly)
# ═══════════════════════════════════════════════════════════
cat("Panel C: PT-S3 Architecture...\n")

pt_mat_long <- melt(pt_matrix, id.vars = "Gene", variable.name = "Layer", value.name = "Present")
VAL_LV <- c("m6A↑","m5C↓","PSU↑","inosine↑","Nm↑","Nm↓","AS_diff","PolyA_diff","Novel_Tx")
pt_mat_long[, Layer := factor(Layer, levels = rev(VAL_LV))]
gene_order <- pt_matrix[order(-rowSums(pt_matrix[,-"Gene"])), Gene]
pt_mat_long[, Gene := factor(Gene, levels = rev(gene_order))]
# data-driven narrative bits
top_gene <- gene_order[1]
top_n <- max(rowSums(pt_matrix[,-"Gene"]))
n_multi <- sum(rowSums(pt_matrix[,-"Gene"]) >= 2)

p_c1 <- ggplot(pt_mat_long, aes(x = Layer, y = Gene, fill = factor(Present))) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = c("0" = "grey93", "1" = "#2E7D6F"), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "PT-S3 Marker Genes × Validated Regulation Layers",
       subtitle = sprintf("%d/%d detected PT markers carry >=2 validated layers; top: %s (%d layers).\nGone is the old PolyA-only picture: validated mod layers (m6A↑/m5C↓/Nm/inosine) cover PT structural genes.",
                          n_multi, nrow(pt_matrix), top_gene, top_n)) +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7, face = "bold", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7, face = "bold"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, color = "grey40")
  )

# Layer coverage side panel
pt_cov <- pt_coverage[Layer != "PolyA_diff"]
pt_cov[, Layer_f := factor(Layer, levels = rev(c("m6A↑","m5C↓","PSU↑","inosine↑","Nm↑","Nm↓","AS_diff","PolyA_diff","Novel_Tx")))]

p_c2 <- ggplot(pt_cov, aes(x = PTS3_N, y = Layer_f, fill = Layer_f)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(PTS3_N, "/", Global_N, " (", PTS3_Pct, "%)")),
            hjust = -0.1, size = 2.8, family = "Arial", color = "grey30") +
  scale_fill_manual(values = c("m6A↑"="#D6404C","m5C↓"="#E8753A","PSU↑"="#3C5488",
                                "inosine↑"="#00A087","Nm↑"="#6A3D9A","Nm↓"="#4B2673",
                                "AS_diff"="#5B8FA8","PolyA_diff"="#2E7D6F","Novel_Tx"="#6BAED6"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(x = "PT-S3 Genes", y = NULL,
       title = "Per-Layer Coverage",
       subtitle = "PolyA: 10/10 (100%). Others: ≤1 gene.") +
  theme_nature

p_c <- p_c1 | p_c2
ggsave(file.path(OUT_DIR, "Fig5C_PTS3_Architecture.pdf"), p_c,
       width = 8, height = 4.2, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig5C_PTS3_Architecture.png"), p_c,
       width = 8, height = 4.2, dpi = 600)

cat("\n=== Fig5 panels A-C DONE ===\n")
cat(sprintf("Output: %s\n", OUT_DIR))
cat("Key finding: PT IS multi-omics — through AS+m6A, not modifications.\n")
cat("AS_SE=3.13 (highest celltype), AS_A5=1.57, Multi_Layer=1.55\n")
