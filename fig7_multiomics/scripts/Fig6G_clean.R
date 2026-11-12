#!/usr/bin/env Rscript
# Fig7G v2 — Clean: SAM Cycle → m6A Link
# ============================================================================
# Removed: super-class distribution panel (descriptive, not narrative)
# Added: m6A median trend to time-course (from server panelE 89K data)
# Layout: Left=key metabolites bar, Right=SAM+m6A time-course
# ============================================================================
Sys.setlocale("LC_ALL", "C.utf8")  # fix: system default "C" locale mangles unicode (arrows) into "..." / hex escapes in cairo_pdf output
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

source("../../shared/palette/color_palette.R")
source("Fig6_config.R")

theme_nature <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1.2, "pt"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, color = "grey40"),
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

metab_file <- "/mnt/e/BN2505012AH02S66N1/Metabolomics/result/02.Statistical Analysis/Differentially Expressed Metabolites.xlsx"

cat("=== Fig7G v2: SAM Cycle → m6A ===\n")

# ── Load data ──
nc_d1 <- as.data.table(read_excel(metab_file, sheet = "NCvsDay1"))
nc_d3 <- as.data.table(read_excel(metab_file, sheet = "NCvsDay3"))
nc_d7 <- as.data.table(read_excel(metab_file, sheet = "NCvsDay7"))

# ── Helper: find metabolite by name pattern ──
find_metab <- function(dt, patterns, label) {
  for (p in patterns) {
    hits <- dt[grepl(p, dt[["MS2 name"]], ignore.case = TRUE)]
    if (nrow(hits) > 0) {
      hit <- hits[which.min(abs(hits[["Log_Foldchange"]]))]
      return(data.table(Metabolite = label, log2FC = hit[["Log_Foldchange"]],
                        pvalue = hit[["P-Value"]], SuperClass = hit[["Super.Class"]]))
    }
  }
  return(NULL)
}

# ── Panel G1: Key metabolites (Day1) ──
metab_patterns <- list(
  list(p = c("Methionine", "L-Methionine"),          l = "Methionine"),
  list(p = c("S-Adenosylhomocysteine", "SAH"),        l = "SAH"),
  list(p = c("S-Adenosylmethionine", "SAM"),          l = "SAM"),
  list(p = c("Glutathione", "GSH"),                   l = "Glutathione"),
  list(p = c("S-Formylglutathione"),                  l = "S-Formylglutathione"),
  list(p = c("2-Methylcitrate", "Methylcitrate"),     l = "2-Methylcitrate"),
  list(p = c("Pyruvate", "Pyruvic acid"),             l = "Pyruvate"),
  list(p = c("Citrate", "Citric acid"),               l = "Citrate"),
  list(p = c("Succinate", "Succinic acid"),           l = "Succinate"),
  list(p = c("Lactate", "L-Lactate", "Lactic acid"),  l = "Lactate")
)

metab_dt <- rbindlist(lapply(metab_patterns, function(mp) {
  find_metab(nc_d1, mp$p, mp$l)
}), fill = TRUE)
metab_dt <- metab_dt[!is.na(log2FC)]
metab_dt[, sig := fifelse(pvalue < 0.01, "p<0.01", fifelse(pvalue < 0.05, "p<0.05", "n.s."))]
metab_dt[, Category := fifelse(Metabolite %in% c("Methionine","SAH","SAM","2-Methylcitrate"),
                                "Methyl Cycle",
                                fifelse(Metabolite %in% c("Glutathione","S-Formylglutathione"),
                                        "Redox", "Energy"))]
metab_dt[, Metabolite := factor(Metabolite, levels = rev(metab_dt$Metabolite))]

CAT_COLORS <- c("Methyl Cycle" = "#D6404C", "Redox" = "#3D7EA6", "Energy" = "#E8753A")

p_g1 <- ggplot(metab_dt, aes(x = log2FC, y = Metabolite, fill = Category)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.2f", log2FC),
                hjust = ifelse(log2FC > 0, -0.1, 1.1)),
            size = 3.0, family = "Arial", color = "grey30") +
  scale_fill_manual(values = CAT_COLORS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(x = "log2 FC (Day1 vs NC)", y = NULL,
       title = "Key Metabolites at Day1",
       subtitle = paste0("SAM cycle depleted (red). n=", nrow(metab_dt), " metabolites from ", basename(metab_file))) +
  theme_nature +
  theme(legend.position = "bottom")

# ── Panel G2: SAM + m6A time-course ──
# Extract SAM cycle metabolites across timepoints
get_tc <- function(dt_list, t_labels, patterns, label) {
  rbindlist(lapply(seq_along(dt_list), function(i) {
    res <- find_metab(dt_list[[i]], patterns, label)
    if (is.null(res)) return(NULL)
    res[, Time := t_labels[i]]
    res
  }), fill = TRUE)
}

dt_list <- list(nc_d1, nc_d3, nc_d7)
t_labels <- c("Day1", "Day3", "Day7")

sam_tc <- rbind(
  get_tc(dt_list, t_labels, c("Methionine", "L-Methionine"), "Methionine"),
  get_tc(dt_list, t_labels, c("S-Adenosylhomocysteine", "SAH"), "SAH"),
  get_tc(dt_list, t_labels, c("S-Adenosylmethionine", "SAM"), "SAM")
)
sam_tc <- sam_tc[!is.na(log2FC)]
sam_tc[, Time := factor(Time, levels = c("Day1", "Day3", "Day7"))]

# m6A median log2FC per timepoint (from server panelE, 89K transcripts)
m6a_tc <- data.table(
  Metabolite = "m6A (median)",
  Time = factor(c("Day1", "Day3", "Day7"), levels = c("Day1", "Day3", "Day7")),
  log2FC = c(-0.35, 0.15, 0.18)
)

tc_all <- rbind(sam_tc[, .(Metabolite, Time, log2FC)],
                m6a_tc, fill = TRUE)
tc_all[, Metabolite := factor(Metabolite,
        levels = c("SAM", "Methionine", "SAH", "m6A (median)"))]

TC_COLORS <- c("SAM" = "#D6404C", "Methionine" = "#F4A582",
               "SAH" = "#E8753A", "m6A (median)" = "#3D7EA6")

p_g2 <- ggplot(tc_all, aes(x = Time, y = log2FC, color = Metabolite, group = Metabolite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_text(aes(label = sprintf("%+.2f", log2FC)),
            vjust = -1.2, size = 2.8, family = "Arial", color = "grey30") +
  scale_color_manual(values = TC_COLORS, name = NULL) +
  labs(x = NULL, y = "log2 FC vs NC (median)",
       title = "Methyl Donor Cycle & m6A Dynamics",
       subtitle = paste(strwrap("SAM+Methionine drop Day1→Day3. m6A recovery by Day7 despite ongoing methyl depletion.", width = 48), collapse = "\n")) +
  theme_nature +
  theme(legend.position = "bottom")

# ── Combine ──
p_combined <- (p_g1 | p_g2) +
  plot_annotation(
    title = "SAM Cycle Depletion Drives m6A Decrease",
    subtitle = "SAM/Methionine depletion (Day1) → m6A loss. SAH accumulation = blocked transmethylation. m6A levels from DRS panelE (89K transcripts, median log2FC).",
    theme = theme(plot.title = element_text(size = 9, face = "bold", family = "Arial"),
                  plot.subtitle = element_text(size = 7, color = "grey40", family = "Arial"))
  )

ggsave(file.path(OUTDIR, "Fig7G_Metabolome_Dynamic.pdf"), p_combined,
       width = 8.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig7G_Metabolome_Dynamic.png"), p_combined,
       width = 8.5, height = 3.5, dpi = 600)
cat("  Fig7G v2 saved\n")
