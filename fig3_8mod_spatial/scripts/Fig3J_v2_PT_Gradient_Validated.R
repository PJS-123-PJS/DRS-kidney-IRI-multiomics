# Fig3 Panel J v2 — PT-subtype gradient, VALIDATED direction sets (GSE139107)
# ============================================================================
# v1 used CDS_m6A_Down headline — that set is the stoichiometry-noise artifact
# (DOWN_all 2.2% DRACH). v2: validated sets from the dual-anchor arbitration:
#   M6A_CDS_UP (DRACH+ HC50, real gain) | M5C_CDS_DOWN (real loss)
#   M6A_CDS_DOWN (79-gene residual, reference) | NM_UP (Nm gain union)
# Data: server validation_GSE139107_5mod_zscore.csv (same scoring as v1:
# mean(colMeans) per cell type -> scale per set). Zero hardcoding.
# Mock: FIG3J_CSV env var.
# ============================================================================

Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
GSE_CSV <- Sys.getenv("FIG3J_CSV",
  "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/validation_GSE139107_5mod_zscore.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_nature <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line=element_line(linewidth=0.3), axis.ticks=element_line(linewidth=0.3),
        axis.ticks.length=unit(1,"pt"), legend.title=element_text(size=6),
        legend.text=element_text(size=5.5), legend.key.size=unit(8,"pt"),
        plot.title=element_text(size=7.5, face="bold"),
        plot.subtitle=element_text(size=6,color="grey30"),
        plot.caption=element_text(size=5.5,color="grey50",hjust=0),
        panel.grid=element_blank())

gse <- fread(GSE_CSV)
setnames(gse, 1, "CellType")

PT_ORDER  <- c("PTS1","PTS2","PTS3","NewPT1","NewPT2")
PT_LABELS <- c(PTS1="PT-S1", PTS2="PT-S2", PTS3="PT-S3",
               NewPT1="NewPT1", NewPT2="NewPT2")
PT_ANAT <- c(PTS1="cortex", PTS2="cortex", PTS3="outer medulla",
             NewPT1="regenerating", NewPT2="regenerating")

pt <- gse[CellType %in% PT_ORDER]
stopifnot(nrow(pt) == length(PT_ORDER))
pt[, Subtype := factor(CellType, levels = PT_ORDER, labels = PT_LABELS)]

CAT_LEVELS <- c("M6A_CDS_UP", "M5C_CDS_DOWN", "NM_UP", "M6A_CDS_DOWN")
stopifnot(all(CAT_LEVELS %in% names(pt)))
CAT_LABELS <- c(M6A_CDS_UP   = "m6A↑ CDS targets (validated)",
                M5C_CDS_DOWN = "m5C↓ CDS targets (validated)",
                NM_UP        = "Nm↑ targets (validated)",
                M6A_CDS_DOWN = "m6A↓ CDS residual (n=79)")
CAT_COLORS <- c(M6A_CDS_UP   = "#D6404C",
                M5C_CDS_DOWN = "#E8753A",
                NM_UP        = "#2E7D6F",
                M6A_CDS_DOWN = "grey60")
CAT_LW     <- c(M6A_CDS_UP = 0.7, M5C_CDS_DOWN = 0.6, NM_UP = 0.5,
                M6A_CDS_DOWN = 0.35)

long <- melt(pt, id.vars = "Subtype",
             measure.vars = CAT_LEVELS,
             variable.name = "Category", value.name = "Z")
long[, Category := factor(Category, levels = CAT_LEVELS, labels = CAT_LABELS)]
CAT_COLORS_L <- setNames(CAT_COLORS[CAT_LEVELS], CAT_LABELS[CAT_LEVELS])
CAT_LW_L     <- setNames(CAT_LW[CAT_LEVELS],    CAT_LABELS[CAT_LEVELS])

long[, xpos := as.numeric(Subtype)]
DIVIDER_X <- 3.5

headline <- CAT_LABELS[["M6A_CDS_UP"]]
peak_row <- long[Category == headline][which.max(Z)]
peak_z   <- peak_row$Z
peak_lab <- as.character(peak_row$Subtype)
newpt_z  <- long[Category == headline & Subtype %in% c("NewPT1","NewPT2"), Z]

y_max <- max(long$Z); y_min <- min(long$Z)

p <- ggplot(long, aes(x = xpos, y = Z, color = Category, group = Category)) +
  annotate("rect", xmin = DIVIDER_X, xmax = length(PT_ORDER) + 0.45,
           ymin = -Inf, ymax = Inf, fill = "grey92", alpha = 0.6) +
  annotate("text", x = (DIVIDER_X + length(PT_ORDER) + 0.45)/2, y = y_max + 0.42,
           label = "recovery phase", size = 2.2, color = "grey40",
           family = "Arial", fontface = "italic") +
  annotate("text", x = 1.85, y = y_max + 0.42,
           label = "proximal tubule axis S1→S3", size = 2.2, color = "grey40",
           family = "Arial", fontface = "italic") +
  geom_vline(xintercept = DIVIDER_X, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70",
             linewidth = 0.25) +
  geom_line(aes(linewidth = Category)) +
  geom_point(size = 1.8, stroke = 0) +
  geom_point(data = long[Category == headline],
             size = 2.6, shape = 21, fill = "white", stroke = 0.8) +
  annotate("segment", x = peak_row$xpos, xend = peak_row$xpos,
           y = peak_z + 0.22, yend = peak_z + 0.62,
           color = "#D6404C", linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  annotate("text", x = peak_row$xpos, y = peak_z + 0.72,
           label = sprintf("peak at %s: Z = %+.2f\n(most IRI-sensitive segment)",
                           peak_lab, peak_z),
           size = 2.3, color = "#D6404C", family = "Arial",
           fontface = "bold", lineheight = 0.9, vjust = 0) +
  scale_color_manual(values = CAT_COLORS_L, name = NULL) +
  scale_linewidth_manual(values = CAT_LW_L, guide = "none") +
  scale_x_continuous(
    breaks = seq_along(PT_ORDER),
    labels = paste0(PT_LABELS[PT_ORDER], "\n", PT_ANAT[PT_ORDER]),
    expand = expansion(mult = c(0.06, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.42))) +
  labs(x = NULL, y = "Enrichment Z-score",
       title = "PT-Subtype Gradient, Validated Direction Sets (GSE139107)",
       subtitle = sprintf(paste0(
         "validated m6A↑ CDS targets: PT-S1→S2→S3 = %+.2f → %+.2f → %+.2f,",
         "\nthen %s in regenerating NewPT1/2; grey = 79-gene m6A↓ residual."),
         pt[CellType=="PTS1", M6A_CDS_UP], pt[CellType=="PTS2", M6A_CDS_UP],
         pt[CellType=="PTS3", M6A_CDS_UP],
         paste(sprintf("%+.2f", newpt_z), collapse = " / ")),
       caption = paste0(
         "Z-scores: enrichment of validated DRS direction sets in GSE139107 (Kirita et al. 2020)",
         "\nPT subtypes. Sets: DRACH+/HC50 m6A gain; gate+/HC50 m5C loss; Nm gain union. PT-S3 = most IRI-sensitive segment.")) +
  theme_nature +
  theme(legend.position = "bottom", legend.justification = "left",
        legend.margin = margin(t = -1),
        axis.text.x = element_text(size = 5.8, lineheight = 0.9),
        plot.subtitle = element_text(size = 6, color = "grey30", lineheight = 1.05))

W_MM <- 100; H_MM <- 82
ggsave(file.path(OUT_DIR, "Fig3J_v2_PT_Gradient_Validated.pdf"), p,
       width = W_MM/25.4, height = H_MM/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3J_v2_PT_Gradient_Validated.png"), p,
       width = W_MM/25.4, height = H_MM/25.4, dpi = 600)
cat(sprintf("Fig3J_v2: %.0f x %.0f mm | headline peak %s Z=%+.2f | NewPT %s\n",
            W_MM, H_MM, peak_lab, peak_z,
            paste(sprintf("%+.2f", newpt_z), collapse="/")))
