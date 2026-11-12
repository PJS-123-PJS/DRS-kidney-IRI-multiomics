# Fig2 Panel K v2: Cell Composition Drift Exclusion (5-mod optimized)
# ============================================================================
# v2 fixes: enzymes colored by 5-mod class | Tet2 removed | Nm defined | Nature theme
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

source("Fig2_config.R")
data_dir <- "../data"
OUTDIR <- "../outputs"; dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# ---- 1. Cell proportion data ----
cell_props <- fread(file.path(data_dir, "cell_proportions_gse274819.csv"))
cell_props <- cell_props[Delta_Pct != 0]

cell_class <- list(
  Epithelial=c("PT","TAL","DCT","CNT_CD","ICB"),
  Immune=c("Macrophage","T_NK","B_cell","Neutrophil"),
  Stromal=c("Endo","Podo"), Proliferating="Proliferating")
class_map <- setNames(rep(names(cell_class), lengths(cell_class)), unlist(cell_class))
cell_props[, Class:=factor(class_map[CellType], levels=names(cell_class))]
cell_props[, CellType:=factor(CellType, levels=cell_props[order(-Delta_Pct), CellType])]
class_colors <- c(Epithelial="#D4A08A", Immune="#8CB5C8", Stromal="#7FB5A0", Proliferating="#B09CC0")

# ---- 2. Enzyme data with 5-mod classification ----
enzyme_dt <- fread(file.path(data_dir, "Table_Enzyme_Day1_FC.csv"))
MOD_COLORS <- c("m6A"="#D6404C","m5C"="#E8753A","PSU"="#3D7EA6","Nm"="#2E7D6F","Other"="grey60")

enzyme_mod <- list(
  m6A=c("Mettl3","Mettl14","Wtap","Rbm15","Zc3h13","Virma","Cbll1","Rbm15b",
        "Fto","Alkbh5","Hnrnpa2b1","Ythdf1","Ythdf2","Ythdf3","Ythdc1","Ythdc2","Hnrnpc","Eif3a"),
  m5C=c("Nsun2","Nsun3","Nsun5","Trmt1","Mettl16","TruB1","TruB2"),
  PSU=c("Pus1","Pus3","Pus7","Pus10","Dkc1","Gar1","Nhp2","Nop10"),
  Nm=c("Fbl","Nop56","Nop58","Ftsj3"),
  Other=c("Tet2"))
mod_map <- setNames(rep(names(enzyme_mod), lengths(enzyme_mod)), unlist(enzyme_mod))
enzyme_dt[, ModClass:=factor(mod_map[enzyme], levels=c("m6A","m5C","PSU","Nm","Other"))]

# Load primary cell type
ct_file <- file.path(data_dir, "enzyme_primary_celltype_server.csv")
if(file.exists(ct_file)) {
  raw <- readLines(ct_file); raw <- raw[!grepl("^cat |^~", raw)]
  enzyme_ct <- fread(paste(raw, collapse="\n"))
} else {
  enzyme_ct <- data.table(enzyme=character(), primary_celltype=character())
}

enzyme_merged <- merge(enzyme_dt, enzyme_ct[,.(enzyme,primary_celltype)], by="enzyme", all.x=TRUE)
enzyme_merged <- merge(enzyme_merged, cell_props[,.(CellType,Delta_Pct)],
                       by.x="primary_celltype", by.y="CellType", all.x=TRUE)

cat(sprintf("Enzymes matched: %d / %d\n",
            sum(!is.na(enzyme_merged$primary_celltype)), nrow(enzyme_merged)))

# ---- Panel K1: Cell proportion lollipop ----
p_bars <- ggplot(cell_props, aes(x=CellType, y=Delta_Pct, color=Class)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey70", linewidth=0.3) +
  geom_segment(aes(xend=CellType, yend=0), linewidth=1.2) +
  geom_point(size=4.5) +
  geom_text(aes(label=sprintf("%+.1f%%",Delta_Pct),
                hjust=ifelse(Delta_Pct>0,-0.2,1.2)), size=3.8, color="grey20", fontface="bold") +
  scale_color_manual(values=class_colors, guide="none") +
  scale_y_continuous(limits=c(-30,18), breaks=seq(-30,15,15)) +
  labs(x=NULL, y=expression(Delta*" Cell Proportion (%)"),
       title="Cell Composition Changes (Sham->Day1)") +
  coord_flip() + theme_drs +
  theme(panel.grid.major.x=element_line(linewidth=0.2, color="grey90"),
        panel.grid.major.y=element_blank(),
        axis.text.y=element_text(size=10), axis.text.x=element_text(size=9),
        plot.title=element_text(size=11, face="bold"), plot.margin=margin(5,15,5,5))

# ---- Panel K2: Enzyme FC vs Cell Composition (5-mod colored) ----
key_label <- c("Fbl","Dkc1","Nop56","Nop58","Pus1","Ftsj3","Nsun2")
enzyme_merged[, Label:=ifelse(enzyme %in% key_label, enzyme, "")]
enzyme_merged[, HasLabel:=enzyme %in% key_label]

p_scatter <- ggplot(enzyme_merged, aes(x=Delta_Pct, y=log2FC)) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=0.15, ymax=Inf, fill="#D6404C", alpha=0.05) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=-0.15, fill="#3D7EA6", alpha=0.05) +
  annotate("text", x=-20, y=1.15, label="Activation despite cell loss",
           size=3.5, color="#D6404C", hjust=0, fontface="italic") +
  annotate("text", x=3, y=-1.15, label="Suppression despite cell gain",
           size=3.5, color="#3D7EA6", hjust=0, fontface="italic") +
  geom_vline(xintercept=0, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_hline(yintercept=0, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_point(data=enzyme_merged[(!HasLabel)], aes(fill=ModClass),
             size=3.5, shape=21, color="white", stroke=0.5, alpha=0.7) +
  geom_point(data=enzyme_merged[(HasLabel)], aes(fill=ModClass),
             size=6, shape=21, color="grey20", stroke=1.0, alpha=1.0) +
  geom_text_repel(data=enzyme_merged[(HasLabel)], aes(label=Label),
                  size=4.2, box.padding=0.6, max.overlaps=15,
                  min.segment.length=0.2, force=5, fontface="bold", color="grey20") +
  scale_fill_manual(values=MOD_COLORS, name="Modification\nClass") +
  scale_x_continuous(limits=c(-30,18)) +
  labs(x="Primary Cell Type Delta Proportion (%)",
       y="Enzyme log2FC (Day1 / NC)",
       title="Enzyme FC vs Cell Composition Change",
       subtitle="Nm=Am+Cm+Gm+Um merged. Tet2 (DNA) excluded from labeling.") +
  theme_drs +
  theme(legend.position="bottom", legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        panel.grid.major=element_line(linewidth=0.2, color="grey92"),
        plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=7.5, color="grey40"))

# ---- Combine ----
combined <- wrap_plots(p_bars, p_scatter, ncol=1, heights=c(2,3)) +
  plot_annotation(
    title="Composition Drift Exclusion: Enzyme Changes Oppose Cell Shifts (5-Mod View)",
    subtitle=paste0("Red zone: enzyme up despite cell loss | Blue zone: enzyme down despite cell gain | ",
                    "Enzymes colored by target modification class. Directional opposition = not cell drift."),
    theme=theme(plot.title=element_text(size=12, face="bold"),
                plot.subtitle=element_text(size=8, color="grey40")))

pdf(file.path(OUTDIR, "Fig2K_Composition_Drift_v2.pdf"), width=10, height=9)
print(combined); dev.off()
png(file.path(OUTDIR, "Fig2K_Composition_Drift_v2.png"), width=10, height=9, units="in", res=600)
print(combined); dev.off()

cat("Fig2K v2 saved (5-mod enzyme classification)\n")
