# Fig4 Muted AS Color Palette — Vision-optimized, colorblind-safe, journal-ready
# Replaces overly bright shared as_colors for Fig4 panels
# ============================================================

as_colors_muted <- c(
  SE   = "#D55E5E",   # Desaturated red (was #E64B35)
  MXE  = "#8CBFD0",   # Muted sky blue (was #4DBBD5)
  A5SS = "#5E8C8C",   # Desaturated teal (was #00A087)
  A3SS = "#E69F5B",   # Muted orange (was #F39B7F)
  RI   = "#9C78B8",   # Muted lavender (was #984EA3)
  AF   = "#FDBB84",   # Soft peach (was #FFD92F)
  AL   = "#B09C85"    # Warm brown (colorblind-safe, avoids red-green conflict)
)

# Name mapping: data AS_type -> palette key
as_name_map <- c("SE"="SE", "RI"="RI", "A3"="A3SS", "A5"="A5SS",
                 "AF"="AF", "AL"="AL", "MX"="MXE")

message("Fig4 muted AS palette loaded. 7 colors, desaturated, colorblind-safe.")
