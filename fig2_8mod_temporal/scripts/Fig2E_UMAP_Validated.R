# Fig2 Panel E — Kidney Cell-Type Reference Map, VALIDATED sets
# =====================================================================
# v2 (retired 2026-07-23) read enrichment_8mod_zscore.csv — wave-era sets.
# v3: PCA on VALIDATED direction-set enrichment profiles (canonical 6 sets,
# same convention as Fig2J v6). Output renamed to Fig2E_UMAP per the
# continuous-lettering re-sync (A–M). Zero hardcoding: hard fail if csv missing.
# =====================================================================
library(data.table); library(ggplot2)
theme_set(theme_classic(base_size=6.5) +
  theme(axis.line=element_line(linewidth=0.35), axis.ticks=element_line(linewidth=0.35)))
OUTDIR <- "../outputs"; dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

cat("=== Fig2E UMAP: Cell-Type Reference Map (validated sets) ===\n")

ENRICH_CSV <- Sys.getenv("FIG2E_ENRICH_CSV", "../results/enrichment_5mod_zscore.csv")
if (!file.exists(ENRICH_CSV)) stop("enrichment_5mod_zscore.csv not found — download from server first")

enrich <- fread(ENRICH_CSV)
setnames(enrich, 1, "CellType")

# canonical validated direction sets (Fig2J v6 convention; m6A↓ residual excluded)
SHOW_SETS <- c("NM_UP","NM_DOWN","INOSINE_UP","PSU_UP","M5C_DOWN","M6A_UP")
miss <- setdiff(SHOW_SETS, names(enrich))
if (length(miss)) stop(paste("missing sets:", paste(miss, collapse=", ")))

enrich_mat <- as.matrix(enrich[, SHOW_SETS, with=FALSE])
rownames(enrich_mat) <- enrich$CellType

# PCA on enrichment profiles
pca <- prcomp(enrich_mat, scale.=TRUE)
pca_scores <- as.data.table(pca$x[, 1:2], keep.rownames="CellType")
setnames(pca_scores, c("CellType","PC1","PC2"))

# Cell class grouping (Fig2J v6 convention)
ct_class <- list(Epithelium=c("PT","TAL","DCT","CNT_CD","ICB"),
  Glomerulus="Podo", Vasculature="Endo",
  Immune=c("Macrophage","T_NK","B_cell","Neutrophil"), Proliferating="Proliferating")
pca_scores[, Class := "Other"]
for (cls in names(ct_class)) {
  pca_scores[CellType %in% ct_class[[cls]], Class := cls]
}
class_colors <- c(Epithelium="#D4A08A", Glomerulus="#8CB5C8", Vasculature="#7FB5A0",
                  Immune="#B09CC0", Proliferating="#E8C870", Other="grey80")

# Add proportion info from cell proportion data
prop_dt <- fread("../data/cell_proportions_gse274819.csv", skip=1)
setnames(prop_dt, names(prop_dt)[1], "CellType")
prop_dt[, TotalPct := (Sham_Pct + Day1_Pct) / 2]
pca_scores <- merge(pca_scores, prop_dt[, .(CellType, TotalPct)], by="CellType", all.x=TRUE)
pca_scores[is.na(TotalPct), TotalPct := 1]

var_pc1 <- round(summary(pca)$importance[2,1] * 100, 1)
var_pc2 <- round(summary(pca)$importance[2,2] * 100, 1)

p <- ggplot(pca_scores, aes(x=PC1, y=PC2, fill=Class)) +
  geom_point(aes(size=TotalPct), shape=21, alpha=0.85, color="white", stroke=0.3) +
  geom_text(aes(label=CellType), size=2.8, vjust=-1.2, color="grey20", fontface="bold") +
  scale_fill_manual(values=class_colors, name="Cell Class") +
  scale_size_continuous(range=c(2, 8), guide="none") +
  labs(x=sprintf("PC1 (%.1f%%)", var_pc1),
       y=sprintf("PC2 (%.1f%%)", var_pc2),
       title="Kidney cell-type reference map",
       subtitle=sprintf("%d cell types from GSE274819 scRNA-seq. PCA on validated modification-set enrichment profiles.",
                        nrow(pca_scores))) +
  theme(legend.position=c(0.12, 0.25),
        legend.background=element_rect(fill=alpha("white", 0.8), color=NA),
        legend.key.size=unit(0.3,"cm"),
        plot.title=element_text(size=8, face="bold"),
        plot.subtitle=element_text(size=6, color="grey40"))

ggsave(file.path(OUTDIR, "Fig2E_UMAP.pdf"), p, width=5, height=4.5, device=cairo_pdf)
ggsave(file.path(OUTDIR, "Fig2E_UMAP.png"), p, width=5, height=4.5, dpi=600)
cat("  Fig2E UMAP (validated) saved.\n")
print(pca_scores[, .(CellType, PC1=round(PC1,2), PC2=round(PC2,2), Class)])
