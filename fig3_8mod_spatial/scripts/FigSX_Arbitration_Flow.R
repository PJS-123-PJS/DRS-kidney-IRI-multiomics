# =============================================================================
# Arbitration flow schematic (supplement candidate)
# 2026-09-27 — one-page visual summary of the dual-anchor arbitration framework:
# input -> complete paired sites -> two anchors -> verdict tree -> outputs,
# with the three robustness gates shown as side checkpoints.
# English labels (submission-ready). Arial 7pt-ish, cairo_pdf.
# Output: output/FigSX_Arbitration_Flow.{pdf,png}
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"

# palette
c_in   <- "#4A4A4A"  # input grey (Sham canonical)
c_d1   <- "#D6404C"  # Day1 canonical red
c_gain <- "#D6404C"
c_loss <- "#2E7D6F"
c_stab <- "#8A8A8A"
c_gate <- "#E8753A"  # Day3 orange for side gates
c_bg   <- "white"

box <- function(x0, x1, y0, y1, label, fill, col = "grey30", tc = "white",
                size = 2.2, bold = FALSE) {
  list(
    annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1,
             fill = fill, color = col, linewidth = 0.3),
    annotate("text", x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label,
             color = tc, size = size, fontface = if (bold) "bold" else "plain",
             lineheight = 0.95)
  )
}
arrow <- function(x0, y0, x1, y1, col = "grey40", curved = 0) {
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1, color = col,
           linewidth = 0.4, arrow = grid::arrow(length = unit(0.08, "inches")))
}

p <- ggplot() +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 10), expand = FALSE) +

  # ---- main vertical flow ---------------------------------------------------
  box(3.6, 8.4, 9.1, 9.9,
      "DRS site-level stoichiometry\nSham (n = 3)  vs  Day 1 (n = 3)",
      c_in, bold = TRUE) +
  arrow(6, 9.1, 6, 8.55) +
  box(3.9, 8.1, 7.75, 8.55,
      "Complete paired sites\n(3/3 + 3/3 detection)",
      "grey55") +
  arrow(5.1, 7.75, 3.1, 7.1) +
  arrow(6.9, 7.75, 8.9, 7.1) +

  # two anchors
  box(1.0, 5.2, 6.15, 7.1,
      "SH anchor\nSham high-confidence sites\ntracked into Day 1\n(paired Wilcoxon)",
      "grey65", tc = "grey10") +
  box(6.8, 11.0, 6.15, 7.1,
      "D1 anchor\nDay 1 high-confidence sites\nback-queried in Sham\n(paired Wilcoxon)",
      c_d1, tc = "white") +
  arrow(3.1, 6.15, 3.1, 5.5) +
  arrow(8.9, 6.15, 8.9, 5.5) +

  # deltas converge
  box(1.4, 4.8, 4.75, 5.5, expression(Delta[SH]~"(pp)"), "grey80", tc = "grey10") +
  box(7.2, 10.6, 4.75, 5.5, expression(Delta[D1]~"(pp)"), "#f3c2c6", tc = "grey10") +
  arrow(3.1, 4.75, 4.9, 4.15) +
  arrow(8.9, 4.75, 7.1, 4.15) +

  # verdict tree
  box(4.4, 7.6, 3.45, 4.15,
      "Verdict tree\nsame sign | opposite sign | small shift",
      "grey20", bold = TRUE) +

  arrow(4.9, 3.45, 2.6, 2.75) +
  arrow(6.0, 3.45, 6.0, 2.75) +
  arrow(7.1, 3.45, 9.4, 2.75) +
  box(1.0, 4.2, 1.95, 2.75,
      "Same sign\nReal gain / Real loss",
      c_loss, bold = TRUE) +
  box(4.6, 7.4, 1.95, 2.75,
      "Opposite sign\nThird witness:\nHC70 count (>=15%)",
      c_gate, bold = TRUE) +
  box(7.8, 11.0, 1.95, 2.75,
      "|delta| < 5 pp (or anchors\nat floor/ceiling)\nStable / not callable",
      c_stab, bold = TRUE) +
  arrow(6.0, 1.95, 6.0, 1.35) +

  # outputs
  box(2.6, 9.4, 0.65, 1.35,
      "Validated direction sets: m6A gain (mild) | m5C dominant loss | Um dominant gain |\ninosine gain | Nm subtype split —  m6A 'loss wave' rejected (detection artifact)",
      c_bg, col = "grey40", tc = "grey10", size = 2.2) +

  # ---- side gates (robustness) ----------------------------------------------
  box(0.15, 2.3, 8.9, 9.7, "Negative control:\nm6A-loss residual set\n(no orthogonal support)",
      "white", col = c_gate, tc = "grey20", size = 1.9) +
  box(0.15, 2.3, 4.55, 5.35, "Animal-level gate:\n3/3 concordant,\nzero group overlap",
      "white", col = c_gate, tc = "grey20", size = 1.9) +
  box(0.15, 2.3, 2.9, 3.7, "Caller independence:\nsibling Nm callers flat,\ncross-modification null",
      "white", col = c_gate, tc = "grey20", size = 1.9) +
  annotate("text", x = 1.22, y = 10.0 + 0, label = "", size = 1) +  # spacer
  annotate("segment", x = 2.3, y = 9.3, xend = 3.6, yend = 9.5,
           color = c_gate, linewidth = 0.35, linetype = 2) +
  annotate("segment", x = 2.3, y = 4.95, xend = 4.4, yend = 3.9,
           color = c_gate, linewidth = 0.35, linetype = 2) +
  annotate("segment", x = 2.3, y = 3.3, xend = 4.4, yend = 3.7,
           color = c_gate, linewidth = 0.35, linetype = 2) +

  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

ggsave(file.path(out_dir, "FigSX_Arbitration_Flow.pdf"), p,
       width = 7.2, height = 5.6, device = cairo_pdf)
ggsave(file.path(out_dir, "FigSX_Arbitration_Flow.png"), p,
       width = 7.2, height = 5.6, dpi = 600)
cat("written:", file.path(out_dir, "FigSX_Arbitration_Flow.{pdf,png}"), "\n")
