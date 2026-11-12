# Fig3 Panel I v2 — Two Writing Mechanisms, Five Regional Signatures
# ============================================================================
# 2026-07-23 — REPLACES v1 (falsified). v1 asserted a two-tier spatial code
# (base mods -> 3'UTR "45-56%", 2'-O-Me -> CDS "49-58%", Wilcoxon p=0.0286)
# tracing to the Jun-era trans_diff chain. Canonical all-sites analysis
# (regional_signature_allsites.R, Table_regional_signature_allsites.csv)
# shows each modification carries its OWN regional signature — no clean
# base-vs-Nm dichotomy (diff-site Wilcoxon p=0.80). v2 renders the REAL
# per-modification regional spectrum read from the canonical CSV at runtime
# (zero hardcoding) + keeps the DREME writing-mechanism annotation (§4.2,
# unaffected) + validated PT-S3 anchor (GSE139107, m6A_UP Z=+2.03).
# HARD FAIL if the canonical CSV is missing — no fallback numbers.
# ============================================================================

Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

OUT_DIR  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
DATA_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CSV <- file.path(DATA_DIR, "Table_regional_signature_allsites.csv")
if (!file.exists(CSV)) stop("HARD FAIL: canonical regional signature CSV missing: ", CSV)
sig <- fread(CSV)

# ---- merge 8 submod rows -> 5 validated identities (Nm = site-weighted pool) ----
pool <- function(dt) {
  n <- sum(dt$n_sites)
  data.table(n = n,
             u5 = sum(dt$utr5_pct * dt$n_sites) / n,
             cd = sum(dt$cds_pct  * dt$n_sites) / n,
             u3 = sum(dt$utr3_pct * dt$n_sites) / n)
}
ident <- rbindlist(list(
  pool(sig[mod == "m6a"])[, identity := "M6A"],
  pool(sig[mod == "m5c"])[, identity := "M5C"],
  pool(sig[mod == "psu"])[, identity := "PSU"],
  pool(sig[mod == "inosine"])[, identity := "INOSINE"],
  pool(sig[type == "2'-O-Me"])[, identity := "Nm"]
))
ident[, identity := factor(identity, levels = c("M6A","M5C","PSU","INOSINE","Nm"))]
setcolorder(ident, "identity")
print(ident)   # transcript-visible provenance

MOD_COLORS <- c(M6A = "#D6404C", M5C = "#E8753A", PSU = "#3D7EA6",
                INOSINE = "#7B5294", Nm = "#2E7D6F")
REGION_FILL <- c("5'UTR" = "#E8D9B0", "CDS" = "#A8C8A0", "3'UTR" = "#9FC3D8")
MECH <- c(M6A = "sequence motif (DRACH)", M5C = "sequence motif",
          PSU = "sequence motif", INOSINE = "ADAR / dsRNA-guided",
          Nm = "snoRNA-guided (Fbl/Nop56)")

theme_nature <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 7.5, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "grey30"),
        plot.caption = element_text(size = 5.5, color = "grey50", hjust = 0,
                                     margin = margin(t = 4), lineheight = 1.15),
        plot.margin = margin(4, 4, 4, 4), panel.grid = element_blank(),
        legend.position = "none")

# ---- canvas 100 x 84mm; bars stacked vertically, one per modification ----
bar_x0 <- 20; bar_x1 <- 88            # bar track (proportional 0-100%)
bar_h  <- 6.2
y_top  <- 62; y_gap <- 11.5           # first bar centre y, vertical pitch

p <- ggplot() + coord_cartesian(xlim = c(0, 100), ylim = c(0, 78), expand = FALSE)

for (i in seq_len(nrow(ident))) {
  m <- as.character(ident$identity[i]); y <- y_top - (i - 1) * y_gap
  u5 <- ident$u5[i]; cd <- ident$cd[i]; u3 <- ident$u3[i]
  w <- bar_x1 - bar_x0
  x5 <- bar_x0 + w * u5 / 100
  xc <- x5     + w * cd / 100
  x3 <- xc     + w * u3 / 100
  p <- p +
    annotate("point", x = 8, y = y, size = 3, color = MOD_COLORS[m]) +
    annotate("text", x = 11.5, y = y, label = ifelse(m == "M6A", "m6A",
              ifelse(m == "M5C", "m5C", ifelse(m == "PSU", "PSU",
              ifelse(m == "INOSINE", "Inosine", "Nm")))),
             size = 2.5, fontface = "bold", hjust = 0, family = "Arial") +
    annotate("text", x = 11.5, y = y - 3.6, label = MECH[m], size = 1.75,
             color = "grey45", hjust = 0, fontface = "italic", family = "Arial") +
    annotate("rect", xmin = bar_x0, xmax = x5, ymin = y - bar_h/2, ymax = y + bar_h/2,
             fill = REGION_FILL["5'UTR"], color = "white", linewidth = 0.25) +
    annotate("rect", xmin = x5, xmax = xc, ymin = y - bar_h/2, ymax = y + bar_h/2,
             fill = REGION_FILL["CDS"], color = "white", linewidth = 0.25) +
    annotate("rect", xmin = xc, xmax = x3, ymin = y - bar_h/2, ymax = y + bar_h/2,
             fill = REGION_FILL["3'UTR"], color = "white", linewidth = 0.25) +
    annotate("text", x = (bar_x0 + x5)/2, y = y, label = sprintf("%.0f%%", u5),
             size = 1.9, color = "grey20", family = "Arial") +
    annotate("text", x = (x5 + xc)/2, y = y, label = sprintf("%.0f%%", cd),
             size = 1.9, color = "grey20", family = "Arial") +
    annotate("text", x = (xc + x3)/2, y = y, label = sprintf("%.0f%%", u3),
             size = 1.9, color = "grey20", family = "Arial")
}
# region legend header over the bar track
p <- p +
  annotate("text", x = bar_x0 + (bar_x1 - bar_x0) * c(0.06, 0.45, 0.86),
           y = y_top + 5.6, label = c("5'UTR", "CDS", "3'UTR"), size = 2.1,
           fontface = "bold", color = "grey30", family = "Arial") +
  annotate("text", x = bar_x0, y = y_top + 9.8,
           label = sprintf("All-sites regional signature (canonical, n = %.1fM-%.1fM per modification)",
                           min(sig$n_sites)/1e6, max(sig$n_sites)/1e6),
           size = 2.0, color = "grey35", hjust = 0, family = "Arial")

# ---- bottom: validated PT-subtype anchor (GSE139107) ----
p <- p +
  annotate("segment", x = 4, xend = 96, y = 8.5, yend = 8.5,
           linewidth = 0.3, color = "grey75") +
  annotate("text", x = 50, y = 4.6,
           label = "Validated anchor: all direction sets rise along PT-S1→S2→S3 and peak at PT-S3 (GSE139107; m6A↑ Z = +2.03, Nm↑ +2.48, inosine↑ +2.42)",
           size = 1.9, color = "grey30", family = "Arial") +
  theme_nature +
  labs(title = "Two Writing Mechanisms, Five Regional Signatures",
       subtitle = "Each modification carries its own constitutive regional code — no base-mod vs 2'-O-Me spatial dichotomy",
       caption = paste0("Region %: canonical all-sites (Table_regional_signature_allsites.csv, all_*.rds Gene_Struc).\n",
                        "Mechanisms: Fig3C DREME + ADAR/snoRNA literature. PT-S3 Z: validated GSE139107 PT-subtype analysis."))

W_MM <- 100; H_MM <- 84
ggsave(file.path(OUT_DIR, "Fig3I_SpatialSignatures_v2.pdf"), p,
       width = W_MM/25.4, height = H_MM/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3I_SpatialSignatures_v2.png"), p,
       width = W_MM/25.4, height = H_MM/25.4, dpi = 600)
cat(sprintf("Fig3I v2 saved: %.0f x %.0f mm\n", W_MM, H_MM))
