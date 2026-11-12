# Fig5J: Multi-Layer Overlap with REAL gene set data
# Computes actual overlap between Novel_Tx host genes and DRS layers
library(data.table); library(ggplot2); library(patchwork)
source("../../shared/palette/color_palette.R")
source("fig5_muted_palette.R")
out_dir <- "../output"

# ── Load gene lists ──────────────────────────────────────────
novel_genes <- unique(readLines("../data/geneset_Novel_Tx_converted.txt"))
novel_genes <- novel_genes[novel_genes != "" & novel_genes != "-"]

as_genes <- unique(readLines("../../fig4_as_landscape/data/geneset_AS_diff.txt"))
as_genes <- as_genes[as_genes != "" & as_genes != "-"]

polya_genes <- unique(readLines("../../fig4_as_landscape/data/geneset_PolyA_diff.txt"))
polya_genes <- polya_genes[polya_genes != "" & polya_genes != "-"]

VDIR <- "../../fig3_8mod_spatial/data/upload_5mod_validated"
rd_set <- function(nm) { g <- unique(readLines(file.path(VDIR, nm))); g[!g %in% c("", "-")] }
m6a_genes  <- rd_set("M6A_UP.txt")        # validated: DRACH+ HC50 gain
m5c_genes  <- rd_set("M5C_DOWN.txt")      # validated: gate+ HC50 loss
psu_genes  <- rd_set("PSU_UP.txt")        # validated: HC50 (modest band gain)
ino_genes  <- rd_set("INOSINE_UP.txt")    # validated: HC50 editing gain
nm_genes   <- rd_set("NM_UP.txt")         # validated: Um/Am/Cm gain union
nmd_genes  <- rd_set("NM_DOWN.txt")       # validated: Gm-led loss union

cat(sprintf("Novel_Tx (converted): %d genes\n", length(novel_genes)))
cat(sprintf("m6A↑(val): %d | m5C↓(val): %d | PSU↑(val): %d | inosine↑(val): %d | Nm↑(val): %d | Nm↓(val): %d\n",
            length(m6a_genes), length(m5c_genes), length(psu_genes),
            length(ino_genes), length(nm_genes), length(nmd_genes)))
ml_genes <- fread("../../fig7_multiomics/data/Table_MultiLayer_Genes_validated.csv")[N_Layers >= 3]$Gene
cat(sprintf("AS_diff: %d | PolyA_diff: %d | Multi-Layer(>=3, validated): %d\n",
            length(as_genes), length(polya_genes), length(ml_genes)))

# ── Compute overlaps ─────────────────────────────────────────
compute_overlap <- function(set1, set2, label) {
  n1 <- length(set1)
  n2 <- length(set2)
  n_overlap <- length(intersect(set1, set2))
  data.table(
    Layer = label,
    Layer_Size = n2,
    Overlap_Novel = n_overlap,
    Pct_Novel = round(n_overlap / n1 * 100, 2),
    Pct_Layer = round(n_overlap / n2 * 100, 2)
  )
}

results <- rbindlist(list(
  compute_overlap(novel_genes, m6a_genes,  "m6A Up"),
  compute_overlap(novel_genes, m5c_genes,  "m5C Down"),
  compute_overlap(novel_genes, ino_genes,  "inosine Up"),
  compute_overlap(novel_genes, nm_genes,   "Nm Up"),
  compute_overlap(novel_genes, nmd_genes,  "Nm Down"),
  compute_overlap(novel_genes, psu_genes,  "PSU Up"),
  compute_overlap(novel_genes, as_genes,   "AS Diff"),
  compute_overlap(novel_genes, polya_genes, "PolyA Diff"),
  compute_overlap(novel_genes, ml_genes,   "Multi-Layer\n(>=3 validated)")
))

cat("\n── Real Overlap Results ──\n")
print(results)

# ── Panel J: Multi-Layer Overlap Bar ─────────────────────────
# Validated direction layers + AS + PolyA (wave framing retired 2026-07-22)
LV <- c("m6A Up","m5C Down","inosine Up","Nm Up","Nm Down","PSU Up","AS Diff","PolyA Diff")
plot_layers <- results[Layer %in% LV]
plot_layers[, Layer := factor(Layer, levels=LV)]

# Color scheme: per-mod family colors, AS/PolyA grey
layer_colors <- c(
  "m6A Up"     = "#D6404C",
  "m5C Down"   = "#E8753A",
  "inosine Up" = "#00A087",
  "Nm Up"      = "#6A3D9A",
  "Nm Down"    = "#4B2673",
  "PSU Up"     = "#3C5488",
  "AS Diff"    = "#5E8C8C",
  "PolyA Diff" = "#9C78B8"
)

p_j <- ggplot(plot_layers, aes(x=Layer, y=Overlap_Novel, fill=Layer)) +
  geom_col(width=0.5, alpha=0.9) +
  geom_text(aes(label=paste0(Overlap_Novel,"\n(",Pct_Layer,"% of layer\n",Pct_Novel,"% of Novel)")),
            vjust=-0.2, size=3.5, lineheight=0.85) +
  scale_fill_manual(values=layer_colors, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.35))) +
  labs(x=NULL, y="Novel Transcript Host Genes Overlapping",
       title="Novel Transcript Host Genes × DRS Modification Layers",
       subtitle=paste0("Based on ", length(novel_genes), " novel transcript host genes (converted to gene symbols)")) +
  theme_drs + theme(
    axis.text.x = element_text(size=9, face="bold"),
    plot.subtitle = element_text(size=8, color="grey50")
  )

# Fig5H output removed (moved to Fig6I — 2026-06-23)
# ── Save overlap data ────────────────────────────────────────
fwrite(results, file.path(out_dir, "Fig6J_overlap_data.csv"))
cat("\nOverlap data saved. Fig5H panel moved to Fig6I (NovelTx_LayerOverlap).\n")
