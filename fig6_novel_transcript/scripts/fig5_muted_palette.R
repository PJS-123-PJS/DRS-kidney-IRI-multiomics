# Fig5 Muted Class Code Palette — Desaturated, colorblind-safe, journal-ready
# Replaces bright class_colors in Fig5 BC EFGHIJ script
# ============================================================

class_colors_muted <- c(
  "k" = "#D55E5E",   # Muted red (was #E64B35)
  "j" = "#8CBFD0",   # Muted sky blue (was #4DBBD5)
  "u" = "#9C78B8",   # Muted lavender (was #984EA3)
  "m" = "#5E8C8C",   # Desaturated teal (was #00A087)
  "x" = "#E69F5B",   # Muted orange (was #F39B7F)
  "i" = "#DCC68B",   # Pale ochre (was #FFD92F — too bright yellow)
  "n" = "#B09C85"    # Warm taupe (was #A6D854 — avoids green overuse)
)

class_labels <- c(
  "k" = "k: contained",
  "j" = "j: multi-exon",
  "u" = "u: intergenic",
  "m" = "m: intron",
  "x" = "x: antisense",
  "i" = "i: intronic",
  "n" = "n: other"
)

# Muted GO/functional colors
go_colors_muted <- c(
  "Stress/IRI"    = "#D55E5E",
  "Cell Death"    = "#E69F5B",
  "Proliferation" = "#8CBFD0",
  "Regulation"    = "#9C78B8",
  "Processing"    = "#5E8C8C",
  "Signaling"     = "#DCC68B",
  "Metabolism"    = "#B09C85",
  "Organelle"     = "#D4C5B9"
)

# Muted heatmap gradient colors
hm_low_muted    <- "#8CBFD0"   # muted blue
hm_mid_muted    <- "white"
hm_high_muted   <- "#D55E5E"   # muted red

# Muted RNA category colors
rna_cat_muted <- c(
  "Known mRNA"   = "#8CBFD0",
  "Known lncRNA" = "#E69F5B",
  "Novel (k)"    = "#D55E5E",
  "Novel (u)"    = "#9C78B8",
  "Novel (other)"= "#D4C5B9"
)

message("Fig5 muted palette loaded. class_colors_muted, go_colors_muted, rna_cat_muted")
