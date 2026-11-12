# DRS Tech Showcase — Unified Color Palette & Theme
# Source: shared/palette/color_palette.R

library(ggplot2)

# ============================================================
# 5 Group Colors (reused from main paper)
# ============================================================
group_colors <- c(
  NC   = "#4DBBD5",
  IS   = "#B09C85",
  Day1 = "#E64B35",
  Day3 = "#F39B7F",
  Day7 = "#7E6148"
)

group_order <- c("NC", "IS", "Day1", "Day3", "Day7")

# ============================================================
# 4 Base Modification Colors (Panel A standard)
# ============================================================
mod_colors <- c(
  m6A     = "#CC5A4C",  # terracotta red
  m5C     = "#E8948A",  # light red
  PSU     = "#00A087",  # green
  Inosine = "#7B3F8C"   # purple
)

# ============================================================
# 2'-O-Me Modification Colors (Panel A standard)
# ============================================================
tome_colors <- c(
  Am = "#8DA0CB",  # light blue
  Cm = "#A8B8D8",  # lighter blue
  Gm = "#3B4F8A",  # deep blue
  Um = "#5C6B9A"   # slate blue
)

# ============================================================
# All 8 Modifications Combined (Panel A standard)
# ============================================================
all_8mod_colors <- c(mod_colors, tome_colors)
all_8mod_order <- c("m6A", "m5C", "PSU", "Inosine", "Gm", "Um", "Am", "Cm")

# ============================================================
# Comparison Group Pairs
# ============================================================
comparison_labels <- c(
  "NC_vs_IS"   = "NC vs IS",
  "NC_vs_Day1" = "NC vs Day1",
  "NC_vs_Day3" = "NC vs Day3",
  "NC_vs_Day7" = "NC vs Day7"
)

# ============================================================
# AS Type Colors
# ============================================================
as_colors <- c(
  SE   = "#E64B35",
  MXE  = "#4DBBD5",
  A5SS = "#00A087",
  A3SS = "#F39B7F",
  RI   = "#984EA3",
  AF   = "#FFD92F",
  AL   = "#A6D854"
)

# ============================================================
# Multi-Omics Layer Colors
# ============================================================
omics_colors <- c(
  Expression    = "#E64B35",
  Modification  = "#4DBBD5",
  AS            = "#00A087",
  PolyA         = "#984EA3",
  ATAC          = "#F39B7F",
  Metabolomics  = "#FFD92F"
)

# ============================================================
# Enzyme Category Colors (harmonized with 8-mod palette)
# ============================================================
enzyme_colors <- c(
  "m6A Writer"      = "#CC5A4C",  # terracotta (m6A)
  "m6A Reader"      = "#E8948A",  # light red (warm family)
  "m6A Eraser"      = "#F4A460",  # sandy brown
  "m5C Enzyme"      = "#7B3F8C",  # purple (Inosine family; distinct from reds)
  "2'-O-Me Writer"  = "#3B4F8A",  # deep blue (Gm)
  "PSU Synthase"    = "#00A087",  # green (PSU)
  "Other"           = "#B09C85"   # neutral grey-brown
)

# ============================================================
# Direction Colors
# ============================================================
direction_colors <- c(
  Up    = "#E64B35",
  Down  = "#4DBBD5",
  NS    = "#B09C85"
)

# ============================================================
# Shared ggplot2 Theme
# ============================================================
theme_drs <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 9),
    plot.title = element_text(color = "black", size = 11, face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# ============================================================
# Figure Export Helper
# ============================================================
save_fig <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "pdf"
  )
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "png"
  )
  message("Saved: ", filename, ".{pdf,png}")
}

# ============================================================
# Shared Data Paths
# ============================================================
PATH_E_DRS   <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/"
PATH_E_RNA   <- "/mnt/e/BN2505012AH02S66N1/RNAseq/Result/"
PATH_E_METAB <- "/mnt/e/BN2505012AH02S66N1/Metabolomics/"
PATH_D_ATAC  <- "/mnt/d/desktop/DRS/ATAC-03.align_location/"
PATH_D_MAIN  <- "/mnt/d/desktop/DRS/正刊-DRS-1/new/"
PATH_D_HERE  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/"
PATH_SHARED   <- file.path(PATH_D_HERE, "shared/data/")

message("DRS Tech Showcase palette loaded.")
message("Available: group_colors, mod_colors, tome_colors, all_8mod_colors, as_colors, omics_colors, enzyme_colors")
message("Theme: theme_drs | Helper: save_fig()")
