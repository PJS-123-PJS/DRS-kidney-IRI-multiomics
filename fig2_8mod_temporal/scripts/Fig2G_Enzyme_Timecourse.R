# Fig2 Panel B: Combined Enzyme Timecourse — All 5 Categories
# Merges former B2 + C2 into single panel (wave framing retired 2026-07-22)
library(ggplot2)
library(dplyr)
library(tidyr)

source("../../shared/palette/color_palette.R")
source("Fig2_config.R")

# ---- Data: TPM matrix (local copy) ----
tpm_file <- "../data/gene.TPM.matrix.xls"
info_file <- "../data/gene.TPM.matrix.add_info.xls"

if (!file.exists(tpm_file)) {
  stop("TPM matrix not found: ", normalizePath(tpm_file), "\n",
       "Copy gene.TPM.matrix.xls from E drive to fig2_8mod_temporal/data/")
}

drs_tpm  <- read.table(tpm_file, header = TRUE, sep = "\t", check.names = FALSE)
drs_info <- read.table(info_file, header = TRUE, sep = "\t", check.names = FALSE,
                       fill = TRUE, quote = "")

colnames(drs_tpm)[1] <- "gene_id"
drs_mat <- as.matrix(drs_tpm[, -1])
rownames(drs_mat) <- drs_tpm$gene_id

gn_map <- drs_info$gene_name
names(gn_map) <- drs_info$ID
gn_map <- gn_map[!is.na(gn_map) & gn_map != "" & gn_map != "-"]

find_gene_id <- function(sym, gn_map) {
  idx <- which(toupper(gn_map) == toupper(sym))
  if (length(idx) > 0) return(names(gn_map)[idx[1]])
  return(NA)
}

tpm_trajectory <- function(gid, drs_mat) {
  if (!gid %in% rownames(drs_mat)) return(rep(NA, 5))
  vals <- drs_mat[gid, ]
  grps <- gsub("s[0-9]+$", "", names(vals))
  grps[grps == "I"] <- "IS"
  sapply(c("NC","IS","Day1","Day3","Day7"), function(g) mean(vals[grps == g]))
}

# Fallback MGP gene IDs for 3 enzymes whose gene symbols are absent from the
# company's add_info.xls gene_name column (verified present in GTF and TPM matrix).
# Verified 2026-07-07: all 3 IDs confirmed in TPM matrix with reasonable expression.
enzyme_fallbacks <- list(
  "Wtap"   = "MGP_C57BL6NJ_G0023567",   # GTF: ENSMUSG00000060475
  "Nop56"  = "MGP_C57BL6NJ_G0026918",   # GTF: ENSMUSG00000027405
  "Ythdf2" = "MGP_C57BL6NJ_G0029288")   # GTF: ENSMUSG00000040025

# ---- All 5 enzyme categories ----
enzyme_categories <- list(
  "m6A Writers"         = c("Mettl3","Mettl14","Wtap","Rbm15","Zc3h13","Virma"),
  "m6A Readers"         = c("Hnrnpa2b1","Hnrnpc","Ythdc1","Ythdc2","Ythdf1","Ythdf2","Ythdf3"),
  "m5C Writers"         = c("Nsun2","Nsun4","Nsun5","Nsun7","Trdmt1"),
  "Nm Writers"     = c("Fbl","Nop56","Nop58","Ftsj3","Ftsj1"),
  "PSU Synthases"       = c("Dkc1","Pus1","Pus3","Pus7","Pus10")
)

# Category colors (shared palette)
cat_colors <- ENZYME_CAT_COLORS

# ---- Compute per-enzyme trajectories ----
enzyme_tpm_list <- list()
for (cat_name in names(enzyme_categories)) {
  for (sym in enzyme_categories[[cat_name]]) {
    gid <- find_gene_id(sym, gn_map)
    if (is.na(gid) && sym %in% names(enzyme_fallbacks)) {
      gid <- enzyme_fallbacks[[sym]]
    }
    if (!is.na(gid) && gid %in% rownames(drs_mat)) {
      traj <- tpm_trajectory(gid, drs_mat)
      enzyme_tpm_list[[length(enzyme_tpm_list) + 1]] <- data.frame(
        Symbol = sym, Category = cat_name, gene_id = gid,
        NC = traj[1], IS = traj[2], Day1 = traj[3], Day3 = traj[4], Day7 = traj[5])
    } else {
      message(sprintf("  NOT FOUND: %s (%s)", sym, cat_name))
    }
  }
}

enzyme_tpm <- do.call(rbind, enzyme_tpm_list)
enzyme_long <- enzyme_tpm %>%
  pivot_longer(c(NC, IS, Day1, Day3, Day7),
               names_to = "Timepoint", values_to = "TPM") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("NC","IS","Day1","Day3","Day7")))

# Normalize each enzyme to its NC value → fold change relative to NC
enzyme_norm <- enzyme_long %>%
  group_by(Symbol, Category) %>%
  mutate(FC_vs_NC = TPM / TPM[Timepoint == "NC"]) %>%
  ungroup() %>%
  mutate(FC_vs_NC = ifelse(is.na(FC_vs_NC) | is.infinite(FC_vs_NC), NA, FC_vs_NC))

# Per-category mean±SD of fold change
cat_stats <- enzyme_norm %>%
  group_by(Category, Timepoint) %>%
  summarise(MeanFC = mean(FC_vs_NC, na.rm = TRUE),
            SDFC = sd(FC_vs_NC, na.rm = TRUE),
            n = n(), .groups = "drop")

# ---- Plot: Combined timecourse (NC-normalized) ----
hline_df <- data.frame(yint = 1, label = "NC baseline")
p_combined <- ggplot(cat_stats, aes(x = Timepoint, y = MeanFC,
                                     color = Category, group = Category,
                                     fill = Category)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(aes(ymin = pmax(0, MeanFC - SDFC), ymax = MeanFC + SDFC),
              alpha = 0.08, color = NA) +
  geom_errorbar(aes(ymin = pmax(0, MeanFC - SDFC), ymax = MeanFC + SDFC),
                width = 0.15, linewidth = 0.35, alpha = 0.5) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_color_manual(values = cat_colors, name = "Enzyme Category") +
  scale_fill_manual(values = cat_colors, guide = "none") +
  labs(title = "Fig2G: RNA Modification Enzyme Timecourse",
       subtitle = "Fold change relative to NC (mean ± SD across enzymes), colored by modification class",
       x = NULL, y = "Expression fold change (NC = 1)") +
  theme_drs + theme(legend.position = c(0.12, 0.85), legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR, "Fig2G_Enzyme_Timecourse.pdf"),
       p_combined, width = 8, height = 5.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig2G_Enzyme_Timecourse.png"),
       p_combined, width = 8, height = 5.5, dpi = 300)

message("Saved: Fig2G_Enzyme_Timecourse.{pdf,png}")

# ---- Print key FC values ----
cat("\nEnzyme Day1 fold changes vs NC (normalized):\n")
for (cat_name in names(enzyme_categories)) {
  d1_fc <- cat_stats$MeanFC[cat_stats$Category == cat_name & cat_stats$Timepoint == "Day1"]
  d1_sd <- cat_stats$SDFC[cat_stats$Category == cat_name & cat_stats$Timepoint == "Day1"]
  cat(sprintf("  %-20s: Day1=%.2f ± %.2f  (%+.0f%%)\n",
      cat_name, d1_fc, d1_sd, (d1_fc - 1) * 100))
}
cat("\nDone.\n")
