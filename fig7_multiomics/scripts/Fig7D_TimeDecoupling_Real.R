# Fig7 Panel A (headline): Temporal Decoupling of Omics Layers — REAL DATA
# ============================================================================
# Replaces deleted Fig7E (which had hardcoded, unverifiable boxplot + trend stats).
# All 4 layers recomputed from real sources with documented thresholds:
#   DRS Modification : unique mod-transcript pairs diffType Down/Up, 8 mods (E-drive 05.methy)
#   DRS Expression   : genes |log2FC|>1 & meanTPM>1 (E-drive gene.TPM.matrix)
#   Metabolome       : metabolites |log2FC|>1 & p<0.05 (company metabolome xlsx)
#   ATAC Chromatin   : genes with promoter(TSS+/-2kb) accessibility |log2FC|>1 & p<0.05,
#                      recomputed from 15 real ATAC BigWigs (no company diff table existed)
# Story: functional layers (modification/expression/metabolome) peak Day1 (ACUTE);
#        chromatin (ATAC) peaks Day3 (DELAYED) => temporal decoupling.
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
source("../../shared/palette/color_palette.R")

OUT <- "../output"
layers3 <- fread("../data/TimeDecoupling_layers.csv")           # DRS_Modification, DRS_Expression, Metabolome
atac    <- fread("../data/ATAC_promoter_timecourse.csv")        # Timepoint, N_Changed_Genes

# ---- Assemble long table, normalize each layer to its own max (%) ----
TPS <- c("IS","Day1","Day3","Day7")
dt <- rbindlist(list(
  data.table(Time=layers3$Timepoint, Layer="DRS Modification", N=layers3$DRS_Modification),
  data.table(Time=layers3$Timepoint, Layer="DRS Expression",   N=layers3$DRS_Expression),
  data.table(Time=layers3$Timepoint, Layer="Metabolome",       N=layers3$Metabolome),
  data.table(Time=atac$Timepoint,    Layer="ATAC Chromatin",   N=atac$N_Changed_Genes)
))
dt[, Time := factor(Time, levels=TPS)]
dt[, TimeN := as.integer(Time)]   # numeric x so annotate() rects/labels align on same scale
dt[, Pct := N / max(N) * 100, by=Layer]
dt[, Layer := factor(Layer, levels=c("DRS Modification","DRS Expression","Metabolome","ATAC Chromatin"))]

LAYER_COL <- c("DRS Modification"="#E8753A","DRS Expression"="#D6404C",
               "Metabolome"="#2E7D6F","ATAC Chromatin"="#3D7EA6")
LAYER_LTY <- c("DRS Modification"="solid","DRS Expression"="solid",
               "Metabolome"="solid","ATAC Chromatin"="22")
LAYER_SHP <- c("DRS Modification"=16,"DRS Expression"=15,   # distinct point shapes so the
               "Metabolome"=17,"ATAC Chromatin"=18)          # three Day1-overlapping lines separate

# peak marker per layer (uses TimeN for numeric x)
peak_dt <- dt[dt[, .I[which.max(Pct)], by=Layer]$V1]

theme_n <- theme_classic(base_size=8, base_family="Arial") +
  theme(axis.line=element_line(linewidth=0.35), axis.ticks=element_line(linewidth=0.35),
        plot.title=element_text(size=10,face="bold"), plot.subtitle=element_text(size=7,color="grey35"),
        legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7),
        legend.key.width=unit(9,"mm"), panel.grid.major.y=element_line(linewidth=0.15,colour="grey92"))

convergence_note <- data.table(TimeN=2.0, Pct=100)

p <- ggplot(dt, aes(x=TimeN, y=Pct, color=Layer, group=Layer, linetype=Layer, shape=Layer)) +
  # acute (Day1) + delayed (Day3) zones
  annotate("rect", xmin=1.5, xmax=2.5, ymin=-3, ymax=131, fill="#D6404C", alpha=0.06) +
  annotate("rect", xmin=2.5, xmax=3.5, ymin=-3, ymax=131, fill="#3D7EA6", alpha=0.06) +
  annotate("text", x=2.0, y=128, label="Acute\nfunctional", color="#D6404C",
           size=2.9, fontface="bold", lineheight=0.9, vjust=1) +
  annotate("text", x=2.72, y=128, label="Delayed\nstructural", color="#3D7EA6",
           size=2.9, fontface="bold", lineheight=0.9, vjust=1, hjust=0) +
  # white halo behind every point so overlapping markers stay separable at a glance
  geom_point(size=4.6, color="white", alpha=0.9) +
  geom_line(linewidth=0.9, alpha=0.9) +
  geom_point(size=2.6, alpha=0.95, stroke=0.6) +
  geom_point(data=peak_dt, size=4.6, shape=21, color="grey20", fill=NA, stroke=1.0) +
  # name the Day1 convergence explicitly instead of leaving 3 overlapping points unexplained
  annotate("text", x=2.0, y=113, label="3 layers converge",
           size=2.5, color="grey30", fontface="italic", vjust=0) +
  annotate("segment", x=2.0, xend=2.0, y=111.5, yend=101.5, linewidth=0.3,
           color="grey50", arrow=arrow(length=unit(1, "mm"), type="closed")) +
  scale_color_manual(values=LAYER_COL) +
  scale_linetype_manual(values=LAYER_LTY) +
  scale_shape_manual(values=LAYER_SHP) +
  scale_x_continuous(breaks=1:4, labels=TPS, expand=expansion(mult=c(0.07,0.07))) +
  scale_y_continuous(limits=c(-3,131), breaks=seq(0,100,25),
                     labels=function(x) paste0(x,"%")) +
  labs(x=NULL, y="Changed features (% of layer max)",
       title="Temporal Decoupling of Omics Layers",
       subtitle=paste0("Epitranscriptome + metabolome peak Day1 (acute); chromatin (ATAC)\n",
                        "peaks Day3 (delayed). Real per-timepoint counts, max-normalized per layer.")) +
  guides(color=guide_legend(nrow=2), linetype=guide_legend(nrow=2), shape=guide_legend(nrow=2)) +
  theme_n +
  theme(plot.margin=margin(5.5, 8, 5.5, 5.5))

ggsave(file.path(OUT,"Fig7A_TimeDecoupling.pdf"), p, width=170/25.4, height=130/25.4, device=cairo_pdf)
ggsave(file.path(OUT,"Fig7A_TimeDecoupling.png"), p, width=170/25.4, height=130/25.4, dpi=600)

# provenance table
prov <- dcast(dt, Time ~ Layer, value.var="N")
fwrite(prov, file.path(OUT,"Fig7A_TimeDecoupling_counts.csv"))
cat("Fig7A saved. Real per-timepoint counts:\n"); print(prov)
cat("\nPeak timepoint per layer:\n"); print(peak_dt[, .(Layer, PeakTime=Time, N)])
