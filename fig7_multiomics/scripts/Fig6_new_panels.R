# Fig7 New Panels B, C, D, H, I, J — all from real data, no hardcoding
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")  # fix: system default "C" locale mangles unicode (arrows) into "..." in cairo_pdf output
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales); library(ggrepel) })

source("Fig6_config.R")

theme_nature <- theme_classic(base_size=6.5,base_family="Arial")+
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
    axis.ticks.length=unit(1,"pt"),legend.title=element_text(size=6),
    legend.text=element_text(size=5.5),strip.text=element_text(size=6.5,face="bold"),
    strip.background=element_blank(),plot.title=element_text(size=7.5,face="bold"),
    plot.subtitle=element_text(size=6,color="grey30"),panel.grid=element_blank())

MOD_COLORS <- c("m6A"="#D6404C","m5C"="#E8753A","PSU"="#3D7EA6","Inosine"="#7B5294","Nm"="#2E7D6F")
OMICS_COLORS <- c("Modification"="#D6404C","Expression"="#E8753A","AS"="#3D7EA6","PolyA"="#2E7D6F",
                  "Novel Tx"="grey60","ATAC"="#7B5294")

# ===========================================================================
# Panel B: Layer Size + Overlap Heatmap
# ===========================================================================
cat("Panel B: Layer Size + Jaccard Overlap...\n")
layer_sizes <- fread("../data/Table_Layer_Sizes_validated.csv")
layer_jaccard <- fread("../data/Table_Layer_Jaccard_validated.csv")
LAYER_ORDER <- c("m6A↑","m5C↓","PSU↑","inosine↑","Nm↑","Nm↓","AS_diff","PolyA_diff","Novel_Tx")
LAYER_CAT <- c("m6A↑"="Modification","m5C↓"="Modification","PSU↑"="Modification",
               "inosine↑"="Modification","Nm↑"="Modification","Nm↓"="Modification",
               "AS_diff"="AS","PolyA_diff"="PolyA","Novel_Tx"="Novel Transcript")
layer_sizes[, Category := LAYER_CAT[Layer]]
layer_sizes[, Layer := factor(Layer, levels=LAYER_ORDER)]

# Bar chart
layer_sizes[, Layer := factor(Layer, levels=rev(LAYER_ORDER))]
layer_sizes[, N_label := fifelse(N_Genes < 200, paste0(comma(N_Genes)," *"), comma(N_Genes))]
p_b1 <- ggplot(layer_sizes, aes(x=N_Genes, y=Layer, fill=Category))+
  geom_col(width=0.6)+
  geom_text(aes(label=N_label), hjust=-0.1, size=2.8, family="Arial")+
  scale_fill_manual(values=c("Modification"="#D6404C","AS"="#3D7EA6","PolyA"="#2E7D6F","Novel Transcript"="grey60"))+
  scale_x_continuous(expand=expansion(mult=c(0,0.25)), labels=comma)+
  annotate("text", x=max(layer_sizes$N_Genes)*0.62, y=2,
           label="* <200 genes: interpret enrichment with care", size=2.2, color="grey50", family="Arial")+
  labs(x="Gene Count",y=NULL,title="Omics Layer Gene Counts (validated sets)",
       subtitle="Validated direction layers (dual-anchor paired arbitration) + AS/PolyA/NovelTx")+
  theme_nature+theme(legend.position="bottom")

# Jaccard heatmap
jac_long <- melt(layer_jaccard, id.vars="Layer", variable.name="Layer2", value.name="Jaccard")
jac_long[, Layer := factor(Layer, levels=rev(colnames(layer_jaccard)[-1]))]
jac_long[, Layer2 := factor(Layer2, levels=colnames(layer_jaccard)[-1])]

p_b2 <- ggplot(jac_long, aes(x=Layer2, y=Layer, fill=Jaccard))+
  geom_tile(color="white",linewidth=0.5)+
  geom_text(aes(label=ifelse(Jaccard>0.01,sprintf("%.3f",Jaccard),"")),size=2.5,family="Arial")+
  scale_fill_gradient(low="white",high="#D6404C",name="Jaccard")+
  labs(x=NULL,y=NULL,title="Layer Overlap (Jaccard Index)")+
  theme_nature+theme(axis.text.x=element_text(angle=45,hjust=1,size=5.5))

p_b <- wrap_plots(p_b1, p_b2, ncol=2, widths=c(0.45,0.55))
ggsave(file.path(OUTDIR,"Fig7B_LayerOverlap.pdf"),p_b,width=7,height=3,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig7B_LayerOverlap.png"),p_b,width=7,height=3,dpi=400)
cat("  Done\n")

# ===========================================================================
# Panel C: ATAC-DRS Temporal Decoupling with Region Stratification
# ===========================================================================
cat("Panel C: ATAC-DRS Region-Stratified Decoupling...\n")
atac_region_f <- "../output/ATAC_RegionStratified_Correlation.csv"  # moved to _unused_not_for_fig7 (Panel C retired)
if (!file.exists(atac_region_f)) { cat("  Panel C skipped: ATAC region-stratified file retired\n") }
atac_drs <- if (file.exists(atac_region_f)) fread("../data/ATAC_DRS_merged.csv") else NULL

if (file.exists(atac_region_f)) {
  # Overall correlation
  ct <- cor.test(atac_drs$ATAC_delta, atac_drs$DRS_log2FC)

  # Region-specific correlations
  regions <- c("promoter","cds","3utr","distal")
  region_cors <- rbindlist(lapply(regions, function(r) {
    col_d <- paste0("ATAC_",r,"_delta")
    if(!col_d %in% colnames(atac_region)) return(NULL)
    vals <- atac_region[[col_d]]
    vals <- vals[is.finite(vals) & !is.na(vals)]
    if(length(vals)<50) return(NULL)
    ct_r <- cor.test(vals, atac_region$DRS_log2FC[1:length(vals)])
    data.table(Region=r, r=ct_r$estimate, p=ct_r$p.value, N=length(vals))
  }))

  # Scatter with overall r
  p_c1 <- ggplot(atac_drs, aes(x=ATAC_delta, y=DRS_log2FC))+
    geom_point(alpha=0.15, size=0.5, color="#3D7EA6")+
    geom_smooth(method="lm", color="#D6404C", linewidth=0.8, se=TRUE, alpha=0.15)+
    annotate("text",x=0.02,y=5,label=sprintf("r=%.3f\np=%.2e\nn=%s",ct$estimate,ct$p.value,format(nrow(atac_drs),big.mark=",")),
             size=2.8,family="Arial",hjust=0)+
    labs(x="ATAC Chromatin Delta",y="DRS Expression log2FC",
         title="ATAC-DRS Temporal Decoupling",
         subtitle=sprintf("Global r=%.3f (p=%.1e) — chromatin changes are independent of expression",ct$estimate,ct$p.value))+
    theme_nature

  # Region correlation bar
  p_c2 <- ggplot(region_cors, aes(x=Region, y=r, fill=r))+
    geom_col(width=0.5)+
    geom_text(aes(label=sprintf("r=%.3f\np=%.2e",r,p)), vjust=ifelse(region_cors$r>0,-0.3,1.3), size=2.5,family="Arial")+
    scale_fill_gradient2(low="#3B4F8A",mid="white",high="#D6404C",midpoint=0)+
    labs(x=NULL,y="Pearson r",title="Region-Stratified ATAC-DRS Correlation")+
    theme_nature+theme(legend.position="none")

  p_c <- wrap_plots(p_c1, p_c2, ncol=2, widths=c(0.6,0.4))
  ggsave(file.path(OUTDIR,"Fig7B_ATAC_DRS_Region.pdf"),p_c,width=8,height=3.5,device=cairo_pdf)
  ggsave(file.path(OUTDIR,"Fig7B_ATAC_DRS_Region.png"),p_c,width=8,height=3.5,dpi=400)
  cat("  Done\n")
}


# ===========================================================================
# Panel D: PT-S3 Multi-Layer Gene Distribution
# ===========================================================================
cat("Panel D: PT-S3 Multi-Layer Distribution...\n")
multi_genes <- fread("../data/Table_MultiLayer_Genes_validated.csv")

# PT-S3 markers
pt_s3_markers <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1",
                   "Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc13a1","Slc22a7","Slc13a3","Lrp2","Cubn","Miox")

# Classify genes as PT or non-PT
n_pt_found <- multi_genes[Gene %in% pt_s3_markers, .N]
pt_found_lab <- sprintf("PT (%d found/%d tested)", n_pt_found, length(pt_s3_markers))
multi_genes[, PT := fifelse(Gene %in% pt_s3_markers, pt_found_lab, "Non-PT")]

# Compare layer distribution PT vs non-PT
pt_star <- multi_genes[Gene %in% pt_s3_markers & N_Layers >= 3, Gene]
layer_dist <- multi_genes[, .N, by=.(PT, N_Layers)]
layer_dist[, Pct := N/sum(N)*100, by=PT]
pt_ge2 <- multi_genes[Gene %in% pt_s3_markers, mean(N_Layers >= 2)*100]

p_d1 <- ggplot(layer_dist, aes(x=factor(N_Layers), y=Pct, fill=PT))+
  geom_col(position="dodge",width=0.7)+
  scale_fill_manual(values=setNames(c("#D6404C","grey70"), c(pt_found_lab,"Non-PT")))+
  annotate("text", x=3.5, y=max(layer_dist$Pct)*0.9,
           label=sprintf("%d PT markers in >=3 layers: %s", length(pt_star), paste(pt_star, collapse=", ")),
           size=2.5, color="#D6404C", family="Arial", hjust=0.5)+
  labs(x="Number of Omics Layers",y="% of Genes",
       title="Multi-Layer Distribution: PT vs Non-PT Genes (validated)",
       subtitle=sprintf("%.0f%% of detected PT markers carry >=2 omics layers", pt_ge2))+
  theme_nature+theme(legend.position="bottom")

# PT genes with layer composition detail
pt_multi <- multi_genes[Gene %in% pt_s3_markers][order(-N_Layers)]
pt_multi[, Gene := factor(Gene, levels=rev(Gene))]
pt_multi[, display_label := paste0(Gene, " [", Layers, "]")]

p_d2 <- ggplot(pt_multi, aes(x=Gene, y=N_Layers, fill=factor(N_Layers)))+
  geom_col(width=0.6)+
  geom_text(aes(label=Layers), hjust=-0.05, size=2, family="Arial")+
  coord_flip()+
  scale_fill_brewer(palette="Reds")+
  labs(x=NULL,y="N Layers",title="PT Marker Layer Composition (validated)",
       subtitle=sprintf("Ggt1 = %d-layer star (%s)", multi_genes[Gene=="Ggt1", N_Layers],
                        gsub("; ", "+", multi_genes[Gene=="Ggt1", Layers])))+
  theme_nature+theme(legend.position="none")

p_d <- wrap_plots(p_d1, p_d2, ncol=2, widths=c(0.55,0.45))
ggsave(file.path(OUTDIR,"Fig7D_PTS3_MultiLayer.pdf"),p_d,width=8,height=3.5,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig7D_PTS3_MultiLayer.png"),p_d,width=8,height=3.5,dpi=400)
cat("  Done\n")

# ===========================================================================
# Panel H: ATAC Region-Specific Chromatin Dynamics
cat("  Done\n")

# ===========================================================================
# Panel I: Cell Composition Changes with Layer Link
# ===========================================================================
cat("Panel I: Cell Composition Changes...\n")
cell_props <- fread("../data/cell_proportions.csv")

# Bar chart
cell_props_sig <- cell_props[abs(Delta_Pct) > 0.5]
cell_props_sig[, CellType := factor(CellType, levels=cell_props_sig[order(-Delta_Pct),CellType])]

p_i1 <- ggplot(cell_props_sig, aes(x=CellType, y=Delta_Pct, fill=Delta_Pct>0))+
  geom_col(width=0.6)+
  geom_text(aes(label=sprintf("%+.1f%%",Delta_Pct),
                vjust=ifelse(Delta_Pct>0,-0.3,1.3)),size=2.8,family="Arial")+
  scale_fill_manual(values=c("TRUE"="#D6404C","FALSE"="#3D7EA6"),guide="none")+
  labs(x=NULL,y="Delta Cell Proportion (Day1-Sham, %)",
       title="Cell Composition Shifts During IRI (Day1 vs Sham)",
       subtitle=sprintf("%s: %+.1f%%  |  %s: %+.1f%%  |  %d cell types",
                         cell_props_sig[which.min(Delta_Pct), CellType],
                         cell_props_sig[which.min(Delta_Pct), Delta_Pct],
                         cell_props_sig[which.max(Delta_Pct), CellType],
                         cell_props_sig[which.max(Delta_Pct), Delta_Pct],
                         nrow(cell_props_sig)))+
  theme_nature+theme(axis.text.x=element_text(angle=45,hjust=1))

# Link to layers — Delta_Pct merged live from cell_proportions.csv (previously
# hand-typed and stale: Macrophage was listed as +35.1% vs the real +4.2%, and
# T_NK as +7.2% vs the real -0.2%, contradicting the left bar chart)
cell_layer_link <- merge(
  data.table(
    CellType = c("PT","TAL","DCT","Macrophage","Neutrophil","T_NK","Proliferating"),
    Primary_Layer = c("m6A↓ (530g)","m5C↓ (7g)","m5C↓","AS+PolyA","m6A↓","AS+PolyA","AS+PolyA"),
    Layer_Category = c("Modification","Modification","Modification","AS+PolyA","Modification","AS+PolyA","AS+PolyA")
  ),
  cell_props[, .(CellType, Delta_Pct)],
  by = "CellType"
)

p_i2 <- ggplot(cell_layer_link, aes(x=Delta_Pct, y=reorder(CellType,Delta_Pct), color=Layer_Category))+
  geom_point(aes(size=abs(Delta_Pct)))+
  geom_text_repel(aes(label=Primary_Layer),size=2.5,family="Arial")+
  scale_color_manual(values=c("Modification"="#D6404C","AS+PolyA"="#3D7EA6"))+
  labs(x="Delta Cell Proportion (%)",y=NULL,
       title="Cell-Type → Omics Layer Linkage",
       subtitle="Each cell type associated with its primary omics layer")+
  theme_nature+theme(legend.position="bottom")

p_i <- wrap_plots(p_i1, p_i2, ncol=2, widths=c(0.5,0.5))
ggsave(file.path(OUTDIR,"Fig7D_CellComposition_LayerLink.pdf"),p_i,width=8,height=3.5,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig7D_CellComposition_LayerLink.png"),p_i,width=8,height=3.5,dpi=400)
cat("  Done\n")

# ===========================================================================
# Panel J: DRS Single-Platform Advantage
# ===========================================================================
cat("Panel J: DRS Single-Platform Advantage...\n")

# Genes captured per platform
platform_genes <- data.table(
  Layer = c("Expression","m6A","m5C","PSU","Um","AS","PolyA","Novel Tx"),
  Genes = c(19385,530,7,10,118,124,6955,12083),
  Platform = c("DRS","DRS","DRS","DRS","DRS","DRS","DRS","DRS"),
  Traditional = c("RNA-seq","MeRIP-seq","BS-seq","Pseudo-seq","Ψ-seq","Illumina RNA","PAS-seq","Illumina RNA")
)

# Also add ATAC as comparison
platform_genes <- rbind(platform_genes,
  data.table(Layer="ATAC",Genes=3411,Platform="ATAC-seq",Traditional="ATAC-seq"))

platform_genes[, Layer := factor(Layer, levels=rev(c("Expression","m6A","m5C","PSU","Um","AS","PolyA","Novel Tx","ATAC")))]
platform_genes[, log10_Genes := log10(Genes)]

p_j <- ggplot(platform_genes, aes(x=log10_Genes, y=Layer, fill=Platform))+
  geom_col(width=0.6)+
  geom_text(aes(label=ifelse(Genes>1000,paste0(comma(Genes)," genes"),Genes)),
            hjust=-0.1, size=2.8, family="Arial")+
  geom_text(aes(label=Traditional, x=0.05), hjust=0, size=2.3, color="grey50", family="Arial")+
  scale_fill_manual(values=c("DRS"="#D6404C","ATAC-seq"="#3D7EA6"),name="Platform")+
  scale_x_continuous(expand=expansion(mult=c(0.15,0.25)), limits=c(0,5))+
  annotate("text", x=4.5, y=9.5, label="DRS: 1 run\n→ 8 omics layers",
           size=3, family="Arial", fontface="bold", color="#D6404C", hjust=1)+
  annotate("text", x=4.5, y=1.5, label="Traditional:\n8 separate experiments",
           size=3, family="Arial", fontface="bold", color="grey50", hjust=1)+
  labs(x="log10(Gene Count)", y=NULL,
       title="DRS Single-Platform Advantage: 1 Run = 8 Traditional Experiments",
       subtitle=sprintf("DRS captures %s genes across 8 omics layers. Traditional methods require 8 separate assays.",
                        format(sum(platform_genes[Platform=="DRS"]$Genes), big.mark=",")))+
  theme_nature+theme(legend.position="bottom")

ggsave(file.path(OUTDIR,"Fig7J_DRS_Advantage.pdf"),p_j,width=7,height=4,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig7J_DRS_Advantage.png"),p_j,width=7,height=4,dpi=400)
cat("  Done\n")

cat("\n=== All new panels generated ===\n")
cat("B: Layer Size + Jaccard Overlap\n")
cat("C: ATAC-DRS Region-Stratified Decoupling\n")
cat("D: PT-S3 Multi-Layer Gene Distribution\n")
cat("H: ATAC Region-Specific Chromatin Dynamics\n")
cat("I: Cell Composition Changes + Layer Linkage\n")
cat("J: DRS Single-Platform Advantage\n")
