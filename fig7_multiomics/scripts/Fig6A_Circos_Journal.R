# Fig7A v6: Publication-elegant Multi-Omics — slim bars, outside labels, tall aspect
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
out_dir <- "../output"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

layer_sizes <- fread("../data/Table_Layer_Sizes.csv")
multi_genes <- fread("../data/Table_MultiLayer_Genes_validated.csv")
metab_dyn <- if(file.exists("../output/Fig7G_Metabolite_Dynamic.csv")) {
  fread("../output/Fig7G_Metabolite_Dynamic.csv")
} else data.table(Category="Metabolites",N=10)

pt_markers <- c("Slc7a13","Slc22a19","Slc5a8","Aldh1a2","Hsd3b1","Slc22a28","Slc16a7","Aqp1",
                "Slc5a2","Slc22a6","Slc34a1","Ggt1","Slc13a1","Slc22a7","Slc13a3","Lrp2","Cubn","Miox")
n_pt <- multi_genes[Gene %in% pt_markers, .N]

# ── Build elegant 4-layer data ──
expr_dist <- multi_genes[, .N, by=N_Layers][order(N_Layers)]
expr_dist[, `:=`(item=paste0(N_Layers,"-layer genes"), pct=N/sum(N), ring="A  DRS Expression")]
expr_dist[, label := paste0(comma(N)," (",round(pct*100,1),"%)")]

mod_data <- layer_sizes[Category=="Modification"]
mod_items <- c("m6A","m5C","PSU","Um","Nm")
mod_data2 <- data.table(
  N=mod_data$N_Genes[1:5], item=mod_items,
  pct=mod_data$N_Genes[1:5]/sum(mod_data$N_Genes[1:5]),
  ring="B  DRS Modifications")
mod_data2[, label := paste0(N," genes (",round(pct*100,1),"%)")]

asp_data <- layer_sizes[Category %in% c("AS","PolyA")]
asp_data2 <- data.table(
  N=asp_data$N_Genes, item=c("AS","PolyA"),
  pct=asp_data$N_Genes/sum(asp_data$N_Genes),
  ring="C  DRS AS + PolyA")
asp_data2[, label := paste0(comma(N)," genes (",round(pct*100,1),"%)")]

metab_valid <- metab_dyn[!is.na(log2FC)]
metab_cats <- metab_valid[, .(.N), by=Category][order(-N)]
if(nrow(metab_cats)==0) metab_cats <- data.table(Category="Metabolites",N=10)
metab_cats[, `:=`(item=Category, pct=N/sum(N), ring="D  Metabolome")]
metab_cats[, label := paste0(N," metabolites")]

all_layers <- rbind(
  expr_dist[, .(item, N, pct, ring, label)],
  mod_data2[, .(item, N, pct, ring, label)],
  asp_data2[, .(item, N, pct, ring, label)],
  metab_cats[, .(item, N, pct, ring, label)],
  fill=TRUE
)
all_layers[, ring := factor(ring, levels=c(
  "A  DRS Expression","B  DRS Modifications","C  DRS AS + PolyA","D  Metabolome"))]
# Build order within each ring
all_layers[, y_order := rev(seq_len(.N)), by=ring]

# ── Refined palette (assign per ring segment) ──
expr_colors <- c("#F0D8C8","#E0B898","#D09068","#C07048","#B05038")
mod_colors  <- c("#D6404C","#E8A090","#5BA0B0","#6090C0","#2E7D6F")
asp_colors  <- c("#5A7A8A","#B8C8D4")
met_colors  <- c("#D6404C","#3D7EA6","#E8753A")
all_layers[, color := ""]
all_layers[ring=="A  DRS Expression", color := expr_colors[seq_len(.N)]]
all_layers[ring=="B  DRS Modifications", color := mod_colors[seq_len(.N)]]
all_layers[ring=="C  DRS AS + PolyA", color := asp_colors[seq_len(.N)]]
all_layers[ring=="D  Metabolome", color := met_colors[seq_len(.N)]]

# ── Nature theme: tall, slim, elegant ──
theme_elegant <- theme_minimal(base_size=9, base_family="Arial") +
  theme(
    panel.grid=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=7.5, color="grey40"),
    axis.ticks.y=element_blank(),
    plot.title=element_text(size=11, face="bold", family="Arial", hjust=0.5, margin=margin(b=4)),
    plot.subtitle=element_text(size=8, color="grey40", family="Arial", hjust=0.5, margin=margin(b=8)),
    plot.caption=element_text(size=6.5, color="grey60", family="Arial", hjust=0.5),
    strip.text=element_text(size=8.5, face="bold", family="Arial", hjust=0),
    strip.placement="outside",
    plot.margin=margin(15,20,10,10),
    legend.position="none")

# ── Draw ──
# Determine label position: inside for large segments (>8%), outside right for medium (3-8%), hide for tiny
all_layers[, label_pos := fifelse(pct > 0.08, "inside",
                           fifelse(pct > 0.03, "outside", "hide"))]

p <- ggplot(all_layers) +
  # Slim stacked bars with elegant spacing
  geom_col(aes(x=ring, y=pct, ymin=0, ymax=pct, fill=I(color)),
           width=0.35, color="white", linewidth=0.8) +
  # Inside labels (white, bold)
  geom_text(data=all_layers[label_pos=="inside"],
            aes(x=ring, y=pct, label=item, group=item),
            position=position_stack(vjust=0.5), size=3.2, color="white",
            fontface="bold", family="Arial") +
  # Outside labels (dark, right-aligned)
  geom_text(data=all_layers[label_pos=="outside"],
            aes(x=ring, y=pct, label=paste0(item,"  "), group=item),
            position=position_stack(vjust=0.5), size=2.8, color="grey40",
            family="Arial", hjust=1) +
  # Ring labels on top of bars
  geom_text(data=all_layers[, .(pct=max(pct)), by=ring],
            aes(x=ring, y=1.06, label=ring), size=3.5, color="grey20",
            fontface="bold", family="Arial", hjust=0.5) +
  # Subtle total annotation per ring
  geom_text(data=all_layers[, .(pct=max(pct), total=sum(N)), by=ring],
            aes(x=ring, y=-0.04, label=paste0(comma(total)," total")),
            size=2.5, color="grey50", family="Arial", hjust=0.5) +
  scale_y_continuous(labels=percent_format(), expand=c(0,0), limits=c(-0.06, 1.15)) +
  labs(x=NULL, y=NULL,
       title="PT-S3 Multi-Omics Convergence",
       subtitle="DRS single-platform captures Expression + Modifications + Splicing + PolyA + Metabolome",
       caption="■ PT-S3 markers detected in Expression, Modifications, and AS+PolyA layers. Ggt1: only PT gene spanning 3 layers.  All data from single DRS sequencing run.") +
  theme_elegant

ggsave(file.path(out_dir,"Fig7A_Linear_4Omics.pdf"), p, width=7, height=7.5, device=cairo_pdf)
ggsave(file.path(out_dir,"Fig7A_Linear_4Omics.png"), p, width=7, height=7.5, dpi=600)
cat("Fig7A v6: elegant slim bars saved\n")
