# Fig2J v6 — 5-mod Cell-Type Enrichment Heatmap, VALIDATED sets
# ============================================================================
# v5 read enrichment_8mod_zscore.csv (noise-era sets: GM_Down/Um_Up/M6A_Down...).
# v6: server scRNA_5mod_Enrichment_Server.R output — 16 validated direction
# sets, mean(colMeans)+scale per set, 1,000 dual permutations (global-high +
# cell-type specificity), BH-FDR stars. HARD FAIL if csv missing.
# Env override: FIG2J_ENRICH_CSV
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("Fig2_config.R")
OUTDIR <- "../outputs"; dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

ENRICH_CSV <- Sys.getenv("FIG2J_ENRICH_CSV", "../results/enrichment_5mod_zscore.csv")
FDR_CSV <- sub("zscore", "permFDR_specific", ENRICH_CSV)
if (!file.exists(ENRICH_CSV)) stop("enrichment_5mod_zscore.csv not found — download from server first")

ez <- fread(ENRICH_CSV); setnames(ez, 1, "CellType")

# Validated headline direction sets (6) — per-class display
SHOW_SETS <- c("NM_UP","NM_DOWN","INOSINE_UP","PSU_UP","M5C_DOWN","M6A_UP")
SET_LABELS <- c(NM_UP="Nm↑",NM_DOWN="Nm↓",INOSINE_UP="inosine↑",PSU_UP="PSU↑",
                M5C_DOWN="m5C↓",M6A_UP="m6A↑")
miss <- setdiff(SHOW_SETS, names(ez))
if (length(miss)) stop(paste("missing sets:", paste(miss, collapse=", ")))

el <- melt(ez[, c("CellType", SHOW_SETS), with=FALSE],
           id.vars="CellType", variable.name="Set", value.name="Zscore")
el[, Mod := factor(SET_LABELS[Set], levels=SET_LABELS[SHOW_SETS])]

if (file.exists(FDR_CSV)) {
  ef <- fread(FDR_CSV); setnames(ef, 1, "CellType")
  fl <- melt(ef[, c("CellType", SHOW_SETS), with=FALSE],
             id.vars="CellType", variable.name="Set", value.name="FDR")
  el <- merge(el, fl, by=c("CellType","Set"), all.x=TRUE)
} else {
  el[, FDR := NA_real_]; warning("permFDR_specific csv missing — no stars")
}
el[, star := fifelse(is.na(FDR), "", fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", "")))]

# Cell-type ordering by nephron/immune class
ct_class <- list(Epithelium=c("PT","TAL","DCT","CNT_CD","ICB"),
  Glomerulus="Podo",Vasculature="Endo",
  Immune=c("Macrophage","T_NK","B_cell","Neutrophil"),Proliferating="Proliferating")
ct_order <- intersect(unlist(ct_class), unique(el$CellType))
extra <- setdiff(unique(el$CellType), ct_order)
ct_order <- c(ct_order, extra)
el[, CellType:=factor(CellType,levels=ct_order)]

theme_nature <- theme_classic(base_size=7,base_family="Arial")+
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
    axis.text.x=element_text(angle=45,hjust=1,size=6.5),
    axis.text.y=element_text(size=6.5),
    plot.title=element_text(size=9,face="bold",hjust=0.5),
    plot.subtitle=element_text(size=7,color="grey30",hjust=0.5),
    plot.caption=element_text(size=6,color="grey50"),
    panel.grid=element_blank())

# Class annotation bar
class_of <- function(ct) names(ct_class)[vapply(ct_class, function(x) ct %in% x, logical(1))]
ct_anno <- data.table(CellType=factor(ct_order,levels=ct_order))
ct_anno[, Class := fifelse(CellType %in% extra, "Other", class_of(as.character(CellType))), by=CellType]
class_colors <- c(Epithelium="#D4A08A",Glomerulus="#8CB5C8",Vasculature="#7FB5A0",
                  Immune="#B09CC0",Proliferating="#E8C870",Other="grey80")

p_class <- ggplot(ct_anno,aes(x=CellType,y=1,fill=Class))+
  geom_tile(color="white",linewidth=0.3)+
  scale_fill_manual(values=class_colors,name=NULL,guide=guide_legend(nrow=1))+
  theme_void()+theme(legend.position="top",legend.key.size=unit(6,"pt"),
    legend.text=element_text(size=6,family="Arial"))

el[, lab := sprintf("%+.1f%s", Zscore, star)]
p_heat <- ggplot(el,aes(x=CellType,y=Mod,fill=Zscore))+
  geom_tile(color="white",linewidth=1)+
  geom_text(aes(label=lab, color=abs(Zscore)>1.5),size=2.6,family="Arial",fontface="bold")+
  scale_color_manual(values=c("TRUE"="white","FALSE"="grey20"),guide="none")+
  scale_fill_gradient2(low="#3B4F8A",mid="white",high="#D6404C",midpoint=0,name="Z-score")+
  labs(x=NULL,y=NULL,
    title="Validated Modification Sets × Cell-Type Enrichment (GSE274819)",
    subtitle="Mean expression of validated direction sets per cell type, z-scaled per set; * FDR<0.05 ** <0.01 (1,000 dual permutations)",
    caption="Sets: dual-anchor paired-arbitration validated directions. m6A↓ residual set excluded (stoichiometry-noise artifact).")+
  theme_nature+theme(legend.position="right",legend.key.height=unit(10,"pt"))

p_out <- p_class/p_heat+plot_layout(heights=c(0.06,1))
ggsave(file.path(OUTDIR,"Fig2J_CellType_Enrichment_v6_validated.pdf"),p_out,width=8,height=4,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig2J_CellType_Enrichment_v6_validated.png"),p_out,width=8,height=4,dpi=600)
cat("Fig2J v6 saved\n")
print(dcast(el, CellType ~ Mod, value.var="Zscore"))
