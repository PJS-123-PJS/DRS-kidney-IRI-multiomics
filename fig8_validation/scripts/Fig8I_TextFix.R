# Fig8 Panel I — meRIP-DRS correlation template, RETIRED-GENE text fix (2026-07-31)
# ============================================================================
# Original template listed the ALL-RETIRED five (Park7/Fis1/Gclm/Gpx4/Trp53 —
# fabrication-chain genes, HC50-verified 0 sites, 2026-07-23 reevaluation).
# Text-only template (Experiment Pending); swap to the validated new five:
# Glyat / Ninj1 / Ehhadh / Tomm5 / Ndufv3.
# Old render archived to archive_superseded_20260731/.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
ARCHIVE <- file.path(OUT_DIR, "archive_superseded_20260731")
dir.create(ARCHIVE, showWarnings = FALSE, recursive = TRUE)
for (ext in c("pdf","png")) {
  src <- file.path(OUT_DIR, paste0("Fig8I_meRIP_DRS_Correlation.", ext))
  if (file.exists(src)) file.rename(src, file.path(ARCHIVE, basename(src)))
}

p_i <- ggplot() +
  annotate("text", x=0.5, y=0.5, label=paste0(
    "meRIP-qPCR × DRS Correlation\n\n",
    "Δ%input (meRIP) vs ΔFraction (DRS)\n",
    "5 genes: Glyat/Ninj1/Ehhadh/Tomm5/Ndufv3\n",
    "(validated set, 2026-07-23 reevaluation;\n",
    " replaces retired Park7/Fis1/Gclm/Gpx4/Trp53)\n\n",
    "Experiment Pending\n",
    "Expected: strong positive correlation\n",
    "(R>0.8), validating DRS as quantitative\n",
    "modification measurement"),
    size=4.0, hjust=0.5, lineheight=1.2) +
  xlim(0,1) + ylim(0,1) + theme_void() +
  labs(title=expression(bold("i")~~"meRIP-DRS Correlation")) +
  theme(plot.title=element_text(size=9,face="bold",hjust=0.5),
        plot.border=element_rect(fill=NA,color="grey70",linewidth=1))

ggsave(file.path(OUT_DIR, "Fig8I_meRIP_DRS_Correlation.pdf"), p_i,
       width = 130/25.4, height = 100/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig8I_meRIP_DRS_Correlation.png"), p_i,
       width = 130/25.4, height = 100/25.4, dpi = 600)
cat("Panel I: template text fixed (retired five -> Glyat/Ninj1/Ehhadh/Tomm5/Ndufv3)\n")
