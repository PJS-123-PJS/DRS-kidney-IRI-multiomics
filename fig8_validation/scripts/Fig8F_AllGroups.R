# Fig8 Panel F — AS RT-PCR real-data panel, ALL-GROUP extension (2026-07-31)
# ============================================================================
# Original Fig8F was a TEXT TEMPLATE (Rab28/Ninj1/Ssbp4, "Experiment Pending")
# listing Day1 dpsi only. User directive: extend wet-lab panels to the full
# time course. This rebuilds F from canonical DRS dpsi files:
#   E:result/07.Structure_Analysis/03.AS/02.Diff_AS/<cmp>-diff/<type>_isoform.dpsi.add_info.xls
# Events = the SAME three shown in Fig4D (narrative consistency, 2026-07-14):
#   Rab28 A3 (chr5:42206019-42214303), Ninj1 AF (chr13:49604853...49605135),
#   Ssbp4 SE (chr8:73556037-73556120).
# Per-sample PSI: NC replicates are identical across the 4 comparison files
# (verified 2026-07-31); each file contributes its treatment group's PSI.
# Plot: PSI trajectory NC->IS->Day1->Day3->Day7, n=3 points + mean±SEM.
# Old template renders archived to output/archive_superseded_20260731/.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

AS_DIR  <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/07.Structure_Analysis/03.AS/02.Diff_AS"
OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
ARCHIVE <- file.path(OUT_DIR, "archive_superseded_20260731")
dir.create(ARCHIVE, showWarnings = FALSE, recursive = TRUE)

# archive previous template renders (kept, not deleted)
for (ext in c("pdf","png")) {
  src <- file.path(OUT_DIR, paste0("Fig8F_AS_RT_PCR_Template.", ext))
  if (file.exists(src)) file.rename(src, file.path(ARCHIVE, basename(src)))
}

EVENTS <- list(
  list(gene = "Rab28", type = "A3", key = "A3:chr5:42206019-42214303:42205245-42214303:-",
       label = "Rab28 — A3 (alt. 3' splice site)"),
  list(gene = "Ninj1", type = "AF", key = "AF:chr13:49604853:49605005-49611104:49605135:49605224-49611104:+",
       label = "Ninj1 — AF (alt. first exon)"),
  list(gene = "Ssbp4", type = "SE", key = "SE:chr8:73556037-73556120:73556179-73556676:-",
       label = "Ssbp4 — SE (exon skipping)")
)
CMPS   <- c("NC_vs_I" = "IS", "NC_vs_Day1" = "Day1", "NC_vs_Day3" = "Day3", "NC_vs_Day7" = "Day7")
GROUPS <- c("NC", "IS", "Day1", "Day3", "Day7")
G_COLORS <- c("NC" = "#4A4A4A", "IS" = "#8A8A8A", "Day1" = "#D6404C",
              "Day3" = "#E8753A", "Day7" = "#2E7D6F")

# ── extract per-sample PSI + dpsi/p per event × comparison ──
psi_long <- list(); annot <- list()
for (ev in EVENTS) {
  for (cmp in names(CMPS)) {
    f <- file.path(AS_DIR, paste0(cmp, "-diff"),
                   paste0(ev$type, "_isoform.dpsi.add_info.xls"))
    dt <- fread(f)
    hit <- dt[grepl(ev$key, event_id, fixed = TRUE)]
    stopifnot(nrow(hit) == 1)
    grp <- CMPS[[cmp]]
    # per-sample PSI columns: NCs1-PSI.. and <grp>s1-PSI.. (IS = Is1..)
    nc_cols  <- paste0("NCs", 1:3, "-PSI")
    g_prefix <- switch(grp, IS = "Is", paste0(grp, "s"))
    g_cols   <- paste0(g_prefix, 1:3, "-PSI")
    psi_long[[length(psi_long)+1]] <- data.table(
      Gene = ev$gene, Label = ev$label, Group = "NC",
      Rep = 1:3, PSI = as.numeric(hit[1, nc_cols, with = FALSE]))
    psi_long[[length(psi_long)+1]] <- data.table(
      Gene = ev$gene, Label = ev$label, Group = grp,
      Rep = 1:3, PSI = as.numeric(hit[1, g_cols, with = FALSE]))
    annot[[length(annot)+1]] <- data.table(
      Gene = ev$gene, Label = ev$label, Group = grp,
      dpsi = hit$dpsi, p = hit$`diffAS-pvalue`)
  }
}
psi <- rbindlist(psi_long)[, Group := factor(Group, levels = GROUPS)]
ann <- rbindlist(annot)
# compact labels: value+stars only ("Δψ vs Sham" stated in subtitle);
# avoid ugly "-0.00" for tiny negative dpsi
ann[, dv := gsub("-0\\.00", "0.00", sprintf("%+.2f", dpsi))]
ann[, lab := paste0(dv,
                    fifelse(p < 0.001, "***", fifelse(p < 0.01, "**",
                    fifelse(p < 0.05, "*", " ns"))))]
cat("=== per-comparison dpsi / significance ===\n")
print(dcast(ann, Gene ~ Group, value.var = "lab"))

# per-group mean±SEM
stat <- psi[, .(m = mean(PSI), se = sd(PSI)/sqrt(.N)), by = .(Label, Group)]

p_f <- ggplot() +
  geom_line(data = stat, aes(x = Group, y = m, group = 1),
            color = "grey60", linewidth = 0.4) +
  geom_errorbar(data = stat, aes(x = Group, ymin = m - se, ymax = m + se),
                width = 0.25, linewidth = 0.3, color = "grey40") +
  geom_point(data = psi, aes(x = Group, y = PSI, color = Group),
             size = 1.6, alpha = 0.9) +
  geom_point(data = stat, aes(x = Group, y = m, color = Group),
             size = 2.6, shape = 18) +
  geom_text(data = ann, aes(x = Group, y = -Inf, label = lab),
            vjust = -0.6, size = 1.9, family = "Arial", color = "grey25") +
  facet_wrap(~ Label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = G_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0.06, 0.02))) +
  labs(x = NULL, y = "PSI (percent spliced-in)",
       title = expression(bold("f")~~"AS RT-PCR: Validated Events — Full Time Course"),
       subtitle = "DRS per-sample PSI (n=3/group, diamond=mean±SEM); Δψ vs Sham with diffAS p (*<0.05 **<0.01 ***<0.001)",
       caption = "Same 3 events as Fig4D. Day1 peak exon exclusion (Rab28 −0.29, Ssbp4 −0.26, Ninj1 −0.14), resolving by Day7.") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.3, color = "grey30"),
        plot.caption = element_text(size = 5.8, color = "grey50", hjust = 0),
        strip.text = element_text(size = 7, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(2, "mm"))

ggsave(file.path(OUT_DIR, "Fig8F_AS_RT_PCR_AllGroups.pdf"), p_f,
       width = 183/25.4, height = 80/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig8F_AS_RT_PCR_AllGroups.png"), p_f,
       width = 183/25.4, height = 80/25.4, dpi = 600)

fwrite(psi,  file.path(OUT_DIR, "Table_Fig8F_per_sample_PSI.csv"))
fwrite(ann,  file.path(OUT_DIR, "Table_Fig8F_dpsi_all_comparisons.csv"))
cat("=== Fig8F all-group real-data panel saved ===\n")
