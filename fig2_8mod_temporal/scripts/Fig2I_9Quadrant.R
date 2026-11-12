# Fig2 Panel I: Modification × Expression 9-Quadrant (Sham vs Day1)
# 2026-09-23: title NC→Sham (in-vivo convention); data path /mnt/e→/mnt/f
# (DRS data drive); QUAD_COLORS unified to project canonical (#D6404C/#8A8A8A).
# =============================================================================
library(ggplot2); library(dplyr)

source("/mnt/d/desktop/DRS/DRS_TechShowcase/shared/palette/color_palette.R")
source("Fig2_config.R")

data_dir <- "/mnt/f/BN2505012AH02S66N1/DirectRNA/05.methy"
stat_file <- file.path(data_dir, "Um", "trans_diff_Expression", "NC_vs_Day1",
                       "NC_vs_Day1.Um_expression.corr.raw.stat.xls")

quad_data <- NULL
if (file.exists(stat_file)) {
  cat("Loading real 9-quadrant statistics...\n")
  raw <- read.table(stat_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # real schema (2026-06-04 redelivery): Diff_Condition = <mod>_<expr> directions
  if (all(c("Diff_Condition", "Number") %in% colnames(raw))) {
    keep <- raw[grepl("^(Up|Down)_(Up|Down)$", raw$Diff_Condition), ]
    lab_map <- c("Up_Up" = "Hyper-mod\n+Up-expr", "Up_Down" = "Hyper-mod\n+Down-expr",
                 "Down_Up" = "Hypo-mod\n+Up-expr", "Down_Down" = "Hypo-mod\n+Down-expr")
    conc_map <- c("Up_Up" = "Concordant", "Down_Down" = "Concordant",
                  "Up_Down" = "Discordant", "Down_Up" = "Discordant")
    quad_data <- data.frame(
      Category = lab_map[keep$Diff_Condition],
      Count = keep$Number,
      Concordance = conc_map[keep$Diff_Condition],
      stringsAsFactors = FALSE)
    # zero-count cells kept explicit so the quadrant shape is visible
    missing_cells <- setdiff(names(lab_map), keep$Diff_Condition)
    if (length(missing_cells)) {
      quad_data <- rbind(quad_data, data.frame(
        Category = lab_map[missing_cells], Count = 0L,
        Concordance = conc_map[missing_cells], stringsAsFactors = FALSE))
    }
    quad_data$Category <- factor(quad_data$Category,
      levels = c("Hyper-mod\n+Up-expr", "Hypo-mod\n+Down-expr",
                 "Hyper-mod\n+Down-expr", "Hypo-mod\n+Up-expr"))
    n_mod_down <- sum(raw$Number[grepl("^Down_", raw$Diff_Condition)])
    cat(sprintf("  real rows: %d | mod-loss transcripts: %d\n", nrow(keep), n_mod_down))
  }
}
if (is.null(quad_data) || nrow(quad_data) == 0) {
  stop("9-quadrant stat file missing or unparseable — no fallback numbers allowed (2026-07-22 rule)")
}

total <- sum(quad_data$Count)
quad_data$Pct <- round(quad_data$Count / total * 100, 1)
quad_data$Category <- factor(quad_data$Category, levels = quad_data$Category)
n_conc <- sum(quad_data$Count[quad_data$Concordance == "Concordant"])
quad_data$Label <- sprintf("%d\n(%.0f%%)", quad_data$Count, quad_data$Pct)

quad_colors <- QUAD_COLORS

p <- ggplot(quad_data, aes(x = Category, y = Count, fill = Concordance)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = Label, y = Count + max(Count) * 0.03),
            size = 5, fontface = "bold", color = "#333333") +
  scale_fill_manual(values = quad_colors, name = "Mod × Expr Direction") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                     limits = c(0, max(quad_data$Count) * 1.15)) +
  labs(x = NULL, y = "Transcripts with both modification & expression change",
       title = "Um Modification × Expression Quadrant (Sham vs Day1, real calls)",
       subtitle = sprintf("n=%d Um-gained transcripts | %.0f%% also expression-up — legacy calls are expression-coupled, motivating paired-stoichiometry validation",
                          total, round(n_conc / total * 100))) +
  theme_drs +
  theme(axis.text.x = element_text(size = 9, face = "bold"),
        legend.position = "bottom", legend.title = element_text(size = 9),
        plot.subtitle = element_text(size = 8))

ggsave(file.path(OUTDIR, "Fig2I_9Quadrant.pdf"), p, width = 8, height = 5)
ggsave(file.path(OUTDIR, "Fig2I_9Quadrant.png"), p, width = 8, height = 5, dpi = 300)
message("Saved: Fig2I_9Quadrant.{pdf,png}")
