# =============================================================================
# meRIP-qPCR cell validation (NC vs H/R) — 3-layer readout, CELL-style mini bars
# 7 genes (rows) x 3 layers (cols) = 21 mini panels, per-gene free y (dynamic
# range up to ~500x makes shared linear/log bars impossible/misleading)
# Extreme-ratio cells render as normal free-y cells (user ruling 09-07: axis
# break effect cancelled; slim bars 0.30 kept). ggbreak is layout-incompatible
# with nested patchwork, do NOT reintroduce.
# Source: data/meRIP_HR_cells_summary_tidy_20260907.csv (user-curated, EH corrected)
# Style: soft palette NC #B8C2CB / HR #E7A199, slim bars (0.30), Arial 7pt, cairo
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

font_family <- "Arial"
nc_fill <- "#B8C2CB"; nc_dark <- "#5F6B76"   # soft gray
hr_fill <- "#E7A199"; hr_dark <- "#BF5F52"   # soft coral

dt <- fread("../data/meRIP_HR_cells_summary_tidy_20260907.csv")

# --- per-replicate stoichiometry = IP / Input ---
wide <- dcast(dt, code + gene + group + rep ~ layer, value.var = "value")
wide[, stoich := meRIP_IP / mRNA_input]
ratio_rows <- wide[, .(code, gene, group, rep, value = stoich, layer = "stoich")]
all <- rbind(dt[, .(code, gene, group, rep, layer, value)], ratio_rows)

# --- orders ---
fc <- ratio_rows[, .(fc = mean(value[group == "HR"]) / mean(value[group == "NC"])), by = gene]
gene_levels <- fc[order(-fc)]$gene                       # stoich FC desc
layer_levels <- c("mRNA_input", "meRIP_IP", "stoich")
layer_titles <- c(
  "mRNA_input" = "a  Input mRNA",
  "meRIP_IP"   = "b  m6A-IP\nenrichment",                      # 09-14: wrapped for narrow cells
  "stoich"     = "c  m6A stoichiometry\n(IP / Input)"          # 09-14: two-line after narrowing
)
all[, gene := factor(gene, levels = gene_levels)]
all[, group := factor(group, levels = c("NC", "HR"))]

# --- Welch t-test on log10 values per gene per layer ---
stats <- all[, {
  p <- tryCatch(t.test(log10(value[group == "HR"]), log10(value[group == "NC"]))$p.value,
                error = function(e) NA_real_)
  .(p = p, fc = mean(value[group == "HR"]) / mean(value[group == "NC"]))
}, by = .(layer, gene)]
stats[, sig := fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "")]
fwrite(stats, "../data/meRIP_HR_cells_stats_20260907.csv")

# --- shared geom stack ---
# 09-14 final geometry: bar width 0.30 on CONTINUOUS x (NC=0, HR=0.6) so the
# NC-HR gap is halved (0.60 -> 0.30 of data units); discrete scale cannot
# compress group spacing. Errorbar cap 0.15, jitter 0.05, x margin add 0.32.
BAR_W   <- 0.30
GROUP_X <- c(NC = 0, HR = 0.6)   # group center distance compressed 1.0 -> 0.6
bar_layers <- function(p, d, sm, sig = "", star_y = NULL) {
  d  <- copy(d)[,  x := unname(GROUP_X[as.character(group)])]
  sm <- copy(sm)[, x := unname(GROUP_X[as.character(group)])]
  p <- p +
    geom_col(data = sm, aes(x = x, y = m, fill = group), width = BAR_W) +
    geom_errorbar(data = sm[group == "NC"], aes(x = x, y = m, ymin = pmax(0, m - sd), ymax = m + sd),
                  color = nc_dark, width = 0.15, linewidth = 0.45) +
    geom_errorbar(data = sm[group == "HR"], aes(x = x, y = m, ymin = pmax(0, m - sd), ymax = m + sd),
                  color = hr_dark, width = 0.15, linewidth = 0.45) +
    geom_point(data = d[group == "NC"], aes(x = x), color = nc_dark, size = 0.9, alpha = 0.85,
               position = position_jitter(width = 0.05, height = 0, seed = 1)) +
    geom_point(data = d[group == "HR"], aes(x = x), color = hr_dark, size = 0.9, alpha = 0.85,
               position = position_jitter(width = 0.05, height = 0, seed = 1)) +
    scale_fill_manual(values = c(NC = nc_fill, HR = hr_fill), name = NULL)
  if (!is.null(star_y))
    p <- p + annotate("text", x = unname(GROUP_X["HR"]), y = star_y, label = sig,
                      family = font_family, size = 2.8, color = "#333333")
  p
}

cell_theme <- function(top = FALSE, left = FALSE, bottom = FALSE) {
  theme_bw(base_size = 7, base_family = font_family) +
    theme(
      panel.border = element_rect(color = "grey75", linewidth = 0.3),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = if (bottom) element_text(size = 6.5) else element_blank(),
      axis.ticks.x = if (bottom) element_line(linewidth = 0.25) else element_blank(),
      axis.text.y = element_text(size = 6, color = "grey30"),
      axis.title.y = if (left) element_text(face = "italic", size = 7.5) else element_blank(),
      plot.title = if (top) element_text(hjust = 0, face = "bold", size = 7) else element_blank(),
      plot.subtitle = if (top) element_text(hjust = 0, size = 6.5, color = "grey45") else element_blank(),
      plot.margin = margin(3, 3, 3, 3),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.32, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )
}

# --- normal mini bar cell ---
make_cell <- function(g, lyr, top = FALSE, left = FALSE, bottom = FALSE) {
  d  <- all[gene == g & layer == lyr]
  s  <- stats[layer == lyr & gene == g]
  sm <- d[, .(m = mean(value), sd = sd(value)), by = group]
  setorder(sm, group)
  star_y <- max(d$value, sm$m + sm$sd) * 1.14
  ymax   <- star_y * 1.18

  p <- ggplot(d, aes(y = value))
  p <- bar_layers(p, d, sm, s$sig, star_y)
  p <- p +
    scale_y_continuous(limits = c(0, ymax), breaks = pretty(c(0, star_y), n = 2),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = unname(GROUP_X), labels = names(GROUP_X),
                       expand = expansion(add = 0.32)) +
    cell_theme(top, left, bottom)
  if (top)    p <- p + labs(title = layer_titles[[lyr]], subtitle = "relative, NC = 1")
  if (left)   p <- p + labs(y = as.character(g))
  p
}

# --- assemble 7 x 3 grid ---
cells <- list()
for (g in gene_levels)
  for (lyr in layer_levels) {
    cells[[length(cells) + 1]] <- make_cell(
      g, lyr,
      top    = (g == gene_levels[1]),
      left   = (lyr == layer_levels[1]),
      bottom = (g == gene_levels[length(gene_levels)])
    )
  }

p <- wrap_plots(cells, ncol = 3, byrow = TRUE) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

out_stub <- "../output/meRIP_HR_CellValidation_3layer"
# 09-14: 6.7 -> 5.0 -> 4.0 in (user: shrink cell frames to hug the bars,
# no need to fill the page); height unchanged
ggsave(paste0(out_stub, ".pdf"), p, width = 4.0, height = 6.9, device = cairo_pdf)
ggsave(paste0(out_stub, ".png"), p, width = 4.0, height = 6.9, dpi = 300, type = "cairo")
cat("saved:", out_stub, "\n")
