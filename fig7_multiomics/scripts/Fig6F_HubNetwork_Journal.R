# Fig7F: Multi-Omics Hub Gene Network — Journal Quality
# =============================================================================
# Optimized: clean network layout, filtered nodes, larger labels,
#            Nature-standard heatmap, consistent typography, 600dpi
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggraph)
  library(scales)
})
source("../../shared/palette/color_palette.R")
out_dir <- "../output"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# Load & prepare
# =============================================================================
cat("Loading data...\n")
ml <- fread("../data/Table_MultiLayer_Genes_validated.csv")
layer_sizes <- fread("../data/Table_Layer_Sizes_validated.csv")
jaccard <- fread("../data/Table_Layer_Jaccard_validated.csv")

# Clean layer names
layer_sizes[, Layer_clean := gsub("↓"," down", gsub("↑"," up", Layer))]
layer_sizes[, Layer_clean := gsub("_diff","", gsub("_Tx","", Layer_clean))]

# Hub gene threshold
hub_genes <- ml[N_Layers >= 3]
setorder(hub_genes, -N_Layers)
cat(sprintf("Hub genes (>=3 layers): %d\n", nrow(hub_genes)))

# =============================================================================
# Network Construction (optimized)
# =============================================================================
cat("Building optimized network...\n")

# Build gene-layer bipartite, project to gene-gene
layer_list <- list()
for (ln in layer_sizes$Layer) {
  genes <- ml[grepl(ln, Layers, fixed=TRUE), Gene]
  if (length(genes) > 0) layer_list[[ln]] <- genes
}

set.seed(42)
gene_to_layer <- rbindlist(lapply(names(layer_list), function(ln) {
  genes <- layer_list[[ln]]
  if (length(genes) > 300) genes <- sample(genes, 300)
  data.table(Gene=genes, Layer=ln)
}))

gene_pairs <- gene_to_layer[, {
  if (.N < 2) return(data.table(Gene1=character(), Gene2=character(), Shared_Layer=character()))
  pairs <- t(combn(Gene, 2))
  data.table(Gene1=pairs[,1], Gene2=pairs[,2], Shared_Layer=Layer[1])
}, by=Layer]

edge_dt <- gene_pairs[, .(N_Shared=.N, Layers=paste(unique(Shared_Layer), collapse="; ")),
                       by=.(Gene1, Gene2)]
edge_dt <- edge_dt[N_Shared >= 1]
setorder(edge_dt, -N_Shared)

# Keep top 400 edges for visual clarity
edge_dt <- edge_dt[1:min(400, .N)]

all_nodes <- unique(c(edge_dt$Gene1, edge_dt$Gene2))
node_meta <- ml[Gene %in% all_nodes, .(Gene, N_Layers, Layers)]

# Build igraph
g <- graph_from_data_frame(edge_dt[, .(Gene1, Gene2, weight=N_Shared)],
                            vertices=node_meta, directed=FALSE)

# Node attributes
V(g)$is_hub <- V(g)$N_Layers >= 3
V(g)$hub_rank <- V(g)$N_Layers
V(g)$label <- V(g)$name
V(g)$label[!V(g)$is_hub] <- ""

# Assign primary layer color
layer_colors_map <- c("m6A up"="#D6404C","m5C down"="#E8753A","PSU up"="#3C5488",
                      "inosine up"="#00A087","Nm up"="#6A3D9A","Nm down"="#4B2673",
                      "AS diff"="#5B8FA8","PolyA diff"="#2E7D6F",
                      "Novel Tx"="#6BAED6")
get_primary <- function(lstr) strsplit(lstr, "; ")[[1]][1]
V(g)$primary <- sapply(V(g)$Layers, get_primary)
# Clean unicode from primary names for matching
primary_clean <- gsub("↓"," down", gsub("↑"," up", V(g)$primary))
primary_clean <- gsub("_diff","", gsub("_Tx","", primary_clean))
V(g)$color <- ifelse(primary_clean %in% names(layer_colors_map),
                     layer_colors_map[primary_clean], "grey50")

cat(sprintf("Network: v=%d, e=%d, hubs=%d\n", vcount(g), ecount(g), sum(V(g)$is_hub)))

# =============================================================================
# F1: Hub Gene Network
# =============================================================================
cat("Plotting F1: Network...\n")

set.seed(42)
p_f1 <- ggraph(g, layout="nicely") +
  geom_edge_link(aes(alpha=weight, width=weight), color="grey65") +
  geom_node_point(aes(size=hub_rank + 1, fill=color, stroke=as.numeric(is_hub)),
                  shape=21, alpha=0.9) +
  geom_node_text(aes(label=label), size=3.2, repel=TRUE,
                 max.overlaps=25, fontface="bold", box.padding=0.8,
                 segment.color="grey50", segment.size=0.3) +
  scale_edge_width_continuous(range=c(0.3, 2), name="Shared\nLayers") +
  scale_edge_alpha_continuous(range=c(0.1, 0.5), guide="none") +
  scale_size_continuous(range=c(2.5, 14), guide="none") +
  scale_fill_identity() +
  scale_color_manual(values=c("0"="grey70","1"="grey10"), guide="none") +
  labs(title="Multi-Omics Hub Gene Network",
       subtitle=sprintf("%d genes in >=2 omics layers | %d hub genes (>=3 layers) | Node size = N layers affected",
                        nrow(ml[N_Layers >= 2]), nrow(hub_genes))) +
  theme_void(base_size=11) +
  theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=9, color="grey40", hjust=0.5),
        legend.position=c(0.92, 0.8),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8))

# =============================================================================
# F2: Top Hub Genes — horizontal bar with layer annotations
# =============================================================================
cat("Plotting F2: Hub gene bars...\n")

top15 <- hub_genes[1:min(15, .N)]
top15[, Rank := .I]
top15[, Gene := factor(Gene, levels=rev(Gene))]
# Parse layers for color strips (validated 9-layer scheme)
top15[, `m6A↑` := grepl("m6A↑", Layers, fixed=TRUE)]
top15[, `m5C↓` := grepl("m5C↓", Layers, fixed=TRUE)]
top15[, `PSU↑` := grepl("PSU↑", Layers, fixed=TRUE)]
top15[, `inosine↑` := grepl("inosine↑", Layers, fixed=TRUE)]
top15[, `Nm↑` := grepl("Nm↑", Layers, fixed=TRUE)]
top15[, `Nm↓` := grepl("Nm↓", Layers, fixed=TRUE)]
top15[, `AS` := grepl("AS_diff", Layers, fixed=TRUE)]
top15[, `PolyA` := grepl("PolyA", Layers, fixed=TRUE)]
top15[, `Novel` := grepl("Novel_Tx", Layers, fixed=TRUE)]

# Create label with layer indicators
top15[, display_label := paste0(Gene, " (", N_Layers, " layers)")]

p_f2 <- ggplot(top15, aes(x=N_Layers, y=Gene)) +
  geom_col(fill="#E64B35", alpha=0.85, width=0.7) +
  geom_text(aes(label=paste0(N_Layers," layers"), x=0.15),
            hjust=0, size=3.5, color="white", fontface="bold") +
  # Layer indicator dots (3x3 grid = 9 validated layers)
  geom_point(aes(x=N_Layers + 0.35, y=as.numeric(Gene) + 0.25, color=`m6A↑`), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.35, y=as.numeric(Gene) + 0.10, color=`m5C↓`), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.35, y=as.numeric(Gene) - 0.05, color=`PSU↑`), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.62, y=as.numeric(Gene) + 0.25, color=`inosine↑`), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.62, y=as.numeric(Gene) + 0.10, color=`Nm↑`), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.62, y=as.numeric(Gene) - 0.05, color=`Nm↓`), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.89, y=as.numeric(Gene) + 0.25, color=AS), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.89, y=as.numeric(Gene) + 0.10, color=PolyA), size=2.6, shape=16) +
  geom_point(aes(x=N_Layers + 0.89, y=as.numeric(Gene) - 0.05, color=Novel), size=2.6, shape=16) +
  scale_color_manual(values=c("TRUE"="#E64B35","FALSE"="grey85"), guide="none") +
  scale_x_continuous(limits=c(0, max(top15$N_Layers)+1.6), expand=c(0,0),
                     breaks=1:max(top15$N_Layers)) +
  labs(x="Number of Omics Layers Affected", y=NULL,
       title="Top 15 Multi-Omics Hub Genes",
       subtitle="Genes regulated across >=3 layers. Dots: specific omics layers affected (right panel)") +
  theme_classic(base_size=6.5,base_family="Arial")+
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
    axis.text=element_text(size=6),strip.text=element_text(size=6.5,face="bold"),
    strip.background=element_blank(),plot.title=element_text(size=8,face="bold"),
    plot.subtitle=element_text(size=6.5,color="grey30"),panel.grid=element_blank(),
    legend.title=element_text(size=6),legend.text=element_text(size=5.5)) +
  theme(axis.text.y=element_text(size=10, face="bold"),
        panel.grid.major.x=element_line(color="grey90"),
        panel.grid.minor=element_blank())

# =============================================================================
# F3: Layer Overlap Heatmap — improved color scale
# =============================================================================
cat("Plotting F3: Jaccard heatmap...\n")
jaccard_mat <- as.matrix(jaccard[, -"Layer", with=FALSE])
rownames(jaccard_mat) <- colnames(jaccard_mat)

keep <- c("m6A↑","m5C↓","PSU↑","inosine↑","Nm↑","Nm↓","AS_diff","PolyA_diff","Novel_Tx")
keep <- intersect(keep, rownames(jaccard_mat))
jac_mat <- jaccard_mat[keep, keep]

# Clean rownames
rn_clean <- gsub("↓"," down", gsub("↑"," up", rownames(jac_mat)))
rn_clean <- gsub("_diff","", rn_clean)
rownames(jac_mat) <- rn_clean
colnames(jac_mat) <- rn_clean

jac_dt <- as.data.table(jac_mat, keep.rownames="Layer1")
jac_long <- melt(jac_dt, id.vars="Layer1", variable.name="Layer2", value.name="Jaccard")
# Ensure Jaccard is numeric
jac_long[, Jaccard := as.numeric(Jaccard)]

mean_j <- mean(jac_mat[upper.tri(jac_mat)])

# Use discrete color breaks for very small values
jac_long[Jaccard == 0, Jaccard_bin := "0"]
jac_long[Jaccard > 0 & Jaccard < 0.005, Jaccard_bin := "0.001-0.005"]
jac_long[Jaccard >= 0.005 & Jaccard < 0.01, Jaccard_bin := "0.005-0.01"]
jac_long[Jaccard >= 0.01 & Jaccard < 0.02, Jaccard_bin := "0.01-0.02"]
jac_long[Jaccard >= 0.02, Jaccard_bin := ">0.02"]

p_f3 <- ggplot(jac_long, aes(x=Layer1, y=Layer2)) +
  geom_tile(aes(fill=Jaccard), color="white", linewidth=1.5) +
  geom_text(aes(label=ifelse(Jaccard==1, "1",
                              ifelse(Jaccard==0, "0", sprintf("%.4f", Jaccard))),
                fontface=ifelse(Jaccard > 0.01, "bold", "plain")),
            size=4) +
  scale_fill_gradientn(
    colors=c("white","#FFF5F0","#FCBBA1","#FB6A4A","#CB181D"),
    values=c(0, 0.005, 0.015, 0.03, 1),
    limits=c(0, max(0.05, max(jac_long$Jaccard[jac_long$Jaccard < 1], na.rm=TRUE))),
    oob=squish,
    name="Jaccard\nIndex"
  ) +
  labs(x=NULL, y=NULL,
       title="Layer Overlap: Jaccard Similarity Index",
       subtitle=sprintf("Mean pairwise Jaccard = %.4f — minimal overlap confirms orthogonal multi-layer regulation",
                        mean_j)) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1, size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        panel.grid=element_blank(),
        legend.position="right",
        legend.text=element_text(size=8),
        plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"))

# =============================================================================
# Assemble Panel F
# =============================================================================
cat("Assembling final figure...\n")

# Top row: Network + Hub Bars
# Bottom: Heatmap
design <- "
AABB
AABB
CCCC
"
p_final <- p_f1 + p_f2 + p_f3 +
  plot_layout(design=design) +
  plot_annotation(
    title="Figure 6F: Multi-Omics Hub Gene Network and Layer Orthogonality",
    subtitle=paste0("DRS resolves 9 validated molecular layers: m6A↑, m5C↓, PSU↑, inosine↑, Nm↑/Nm↓ | AS | PolyA | Novel transcripts. ",
                    sprintf("%d genes affected in >=2 layers, %d hub genes in >=3 layers. Jaccard %.4f confirms orthogonal regulation.",
                            nrow(ml[N_Layers>=2]), nrow(hub_genes), mean_j)),
    theme=theme(plot.title=element_text(size=15, face="bold", hjust=0.5),
                plot.subtitle=element_text(size=9, color="grey40", hjust=0.5))
  )

# Export
pdf(file.path(out_dir, "Fig7F_HubGene_Network_Journal.pdf"), width=16, height=14)
print(p_final)
dev.off()

png(file.path(out_dir, "Fig7F_HubGene_Network_Journal.png"), width=16, height=14, units="in", res=600)
print(p_final)
dev.off()

# Also overwrite 300dpi version
png(file.path(out_dir, "Fig7F_HubGene_Network.png"), width=16, height=14, units="in", res=300)
print(p_final)
dev.off()

# Save hub gene table
fwrite(hub_genes, file.path(out_dir, "Fig7F_HubGenes_Table.csv"))

cat(sprintf("\nFig7F Journal-quality complete.\n"))
cat(sprintf("Key: %d hub genes, mean Jaccard=%.4f (orthogonal layers)\n",
            nrow(hub_genes), mean_j))
