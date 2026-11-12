# Fig7 Panel E: Temporal Decoupling, Decoded by Cell Type — schematic/synthesis panel
# ============================================================================
# Not a new analysis. Synthesizes two already-real findings that sit in separate
# panels and are never drawn as one causal chain:
#   - Fig7A: DRS Modification/Expression/Metabolome peak Day1; ATAC peaks Day3
#            (real per-timepoint counts, TimeDecoupling_layers.csv / ATAC_promoter_timecourse.csv)
#   - Fig7D + Results §8.3: "acute layer (DRS mod+expr) mainly from PT + immune
#            cells; delayed layer (ATAC) mainly from endothelial + fibroblast cells"
# Metabolome is bulk-tissue (not cell-type resolved), so it is shown on the
# timeline but deliberately NOT wired into the cell-type arrows below — the
# manuscript's cell-type decoding (§8.3) only covers DRS mod+expr and ATAC.
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(ggplot2) })
source("../../shared/palette/color_palette.R")

OUT <- "../output"

LAYER_COL <- c("DRS Modification"="#E8753A","DRS Expression"="#D6404C",
               "Metabolome"="#2E7D6F","ATAC Chromatin"="#3D7EA6")
ACUTE_RED  <- "#D6404C"
DELAY_BLUE <- "#3D7EA6"

theme_n <- theme_classic(base_size=8, base_family="Arial") +
  theme(axis.line=element_blank(), axis.ticks=element_blank(),
        axis.text=element_blank(), axis.title=element_blank(),
        plot.title=element_text(size=10.5, face="bold"),
        plot.subtitle=element_text(size=7.3, color="grey35"),
        plot.caption=element_text(size=5.6, color="grey50", hjust=0,
                                   margin=margin(t=5), lineheight=1.15),
        plot.margin=margin(6,8,6,6),
        panel.grid=element_blank())

p <- ggplot() +
  coord_cartesian(xlim=c(0,180), ylim=c(0,122), expand=FALSE) +

  # ---- shared timeline (same acute/delayed convention as Fig7A) ----
  annotate("rect", xmin=56, xmax=88,  ymin=70, ymax=101, fill=ACUTE_RED,  alpha=0.07) +
  annotate("rect", xmin=96, xmax=128, ymin=70, ymax=101, fill=DELAY_BLUE, alpha=0.07) +
  annotate("segment", x=20, xend=160, y=74, yend=74, linewidth=0.4, color="grey40") +
  annotate("text", x=c(32,72,112,152), y=71.5, label=c("IS","Day1","Day3","Day7"),
           size=2.5, color="grey30", family="Arial") +
  annotate("segment", x=c(32,72,112,152), xend=c(32,72,112,152), y=73, yend=75,
           linewidth=0.35, color="grey40") +
  annotate("text", x=72, y=99.5, label="Acute", color=ACUTE_RED, size=2.9,
           fontface="bold") +
  annotate("text", x=112, y=99.5, label="Delayed", color=DELAY_BLUE, size=2.9,
           fontface="bold") +

  # timeline markers: 3 acute layers spread inside the red zone, ATAC alone in blue zone
  annotate("point", x=c(63,72,81,112), y=rep(85,4), size=3.2,
           color=LAYER_COL[c("DRS Modification","DRS Expression","Metabolome","ATAC Chromatin")]) +
  annotate("segment", x=c(63,72,81,112), xend=c(63,72,81,112), y=85, yend=76,
           linewidth=0.4, color=LAYER_COL[c("DRS Modification","DRS Expression","Metabolome","ATAC Chromatin")],
           alpha=0.55) +

  # ---- left-side legend: what each color/marker means, with real counts ----
  annotate("point", x=10, y=57.5, size=2.6, color=LAYER_COL["DRS Modification"]) +
  annotate("text", x=13.5, y=57.5, label="DRS Modification — 1,180 mod-transcripts (Day1 peak)",
           size=2.3, hjust=0, family="Arial", color="grey20") +
  annotate("point", x=10, y=53.3, size=2.6, color=LAYER_COL["DRS Expression"]) +
  annotate("text", x=13.5, y=53.3, label="DRS Expression — 1,377 DEGs (Day1 peak)",
           size=2.3, hjust=0, family="Arial", color="grey20") +
  annotate("point", x=10, y=49.1, size=2.6, color=LAYER_COL["Metabolome"]) +
  annotate("text", x=13.5, y=49.1, label="Metabolome — 8,809 metabolites (Day1 peak, bulk tissue)",
           size=2.3, hjust=0, family="Arial", color="grey20") +
  annotate("point", x=10, y=44.9, size=2.6, color=LAYER_COL["ATAC Chromatin"]) +
  annotate("text", x=13.5, y=44.9, label="ATAC Chromatin — 756 genes (Day3 peak; r=0.214 w/ DRS, Fig7B)",
           size=2.3, hjust=0, family="Arial", color="grey20") +

  # ---- cell-type decode: acute cluster -> PT+Immune ; ATAC -> Endo+Fibroblast ----
  annotate("segment", x=c(63,72), xend=68, y=85, yend=38, linewidth=0.35,
           color=ACUTE_RED, alpha=0.5,
           arrow=arrow(length=unit(1.4,"mm"), type="closed")) +
  annotate("segment", x=112, xend=112, y=85, yend=38, linewidth=0.35,
           color=DELAY_BLUE, alpha=0.5,
           arrow=arrow(length=unit(1.4,"mm"), type="closed")) +
  # metabolome deliberately gets a short dotted stub instead of an arrow into the
  # cell-type boxes below — bulk-tissue data was not part of the §8.3 cell-type mapping
  annotate("segment", x=81, xend=81, y=85, yend=76, linewidth=0.3, linetype="dotted",
           color=LAYER_COL["Metabolome"], alpha=0.6) +
  annotate("text", x=81, y=74.3, label="no cell-type link",
           size=1.8, color=LAYER_COL["Metabolome"], fontface="italic", family="Arial") +

  annotate("rect", xmin=38, xmax=94, ymin=18, ymax=36, fill=ACUTE_RED, alpha=0.08,
           color=ACUTE_RED, linewidth=0.3) +
  annotate("text", x=66, y=30.5, label="PT + Immune cells", size=3.1, fontface="bold",
           color=ACUTE_RED, family="Arial") +
  annotate("text", x=66, y=25.5, label="drive the acute functional layer",
           size=2.3, color="grey30", fontface="italic", family="Arial") +
  annotate("text", x=66, y=21.5, label="(Fig7D cell-type mapping)",
           size=2.0, color="grey50", family="Arial") +

  annotate("rect", xmin=98, xmax=154, ymin=18, ymax=36, fill=DELAY_BLUE, alpha=0.08,
           color=DELAY_BLUE, linewidth=0.3) +
  annotate("text", x=124, y=30.5, label="Endothelial + Fibroblast", size=3.1,
           fontface="bold", color=DELAY_BLUE, family="Arial") +
  annotate("text", x=124, y=25.5, label="drive the delayed chromatin layer",
           size=2.3, color="grey30", fontface="italic", family="Arial") +
  annotate("text", x=124, y=21.5, label="(Fig7D cell-type mapping)",
           size=2.0, color="grey50", family="Arial") +

  labs(title="Temporal Decoupling, Decoded by Cell Type",
       subtitle=paste0("Four measured layers converge on one timing pattern (Fig7A);\n",
                        "cell-type mapping (Fig7D) assigns the timing of each layer to specific cell types"),
       caption="Peak counts: Fig7A. ATAC-DRS coupling uses raw BigWig coverage delta, not DESeq2 lfcShrink (Results §8.2).") +
  theme_n

ggsave(file.path(OUT,"Fig7E_TemporalCellType_Model.pdf"), p, width=180/25.4, height=122/25.4, device=cairo_pdf)
ggsave(file.path(OUT,"Fig7E_TemporalCellType_Model.png"), p, width=180/25.4, height=122/25.4, dpi=600)
cat("Fig7E saved: 180 x 122 mm\n")
