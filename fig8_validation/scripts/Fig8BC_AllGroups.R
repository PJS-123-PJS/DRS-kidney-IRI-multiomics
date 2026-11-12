# Fig8 Panels B & C — ALL-GROUP extension (2026-07-31)
# ============================================================================
# User directive: wet-lab panels B (qRT-PCR enzymes) and C (WB targets) showed
# only Sham vs Day1 — extend to the full time course (IS/Day1/Day3/Day7 vs NC).
# Data: local DRS TPM matrix (fig2_8mod_temporal/data/gene.TPM.matrix.xls —
# byte-identical company delivery copy) + add_info symbols + the 2026-07-19
# MGP ID rescue (Wtap/Ythdf2/Nop56, verified against uniprot column).
# Same log2FC convention as the original Fig7_all_panels.R:
#   log2((group_mean+0.1)/(NC_mean+0.1)), per-sample TPM, group means.
# Old Sham/Day1-only renders archived to output/archive_superseded_20260731/.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

TPM_F    <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig2_8mod_temporal/data/gene.TPM.matrix.xls"
ANNOT_F  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig2_8mod_temporal/data/gene.TPM.matrix.add_info.xls"
OUT_DIR  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig8_validation/output"
ARCHIVE  <- file.path(OUT_DIR, "archive_superseded_20260731")
dir.create(ARCHIVE, showWarnings = FALSE, recursive = TRUE)

# archive previous Sham/Day1-only renders (kept, not deleted)
for (f in c("Fig8B_qRT_PCR_Enzymes","Fig8C_Western_Blot_Targets")) {
  for (ext in c("pdf","png")) {
    src <- file.path(OUT_DIR, paste0(f, ".", ext))
    if (file.exists(src)) file.rename(src, file.path(ARCHIVE, paste0(f, ".", ext)))
  }
}

# ── load TPM + symbols (with rescue) ──
tpm <- fread(TPM_F)
annot <- fread(ANNOT_F)
if ("gene_name" %in% names(annot)) setnames(annot, "gene_name", "Gene")
tpm <- merge(tpm, annot[, .(ID, Gene)], by = "ID", all.x = TRUE)
tpm[, Gene := ifelse(is.na(Gene) | Gene == "-", ID, Gene)]
gene_id_rescue <- data.table(
  ID   = c("MGP_C57BL6NJ_G0023567","MGP_C57BL6NJ_G0029288","MGP_C57BL6NJ_G0026918"),
  Gene = c("Wtap","Ythdf2","Nop56"))
tpm[gene_id_rescue, Gene := i.Gene, on = "ID"]

GROUPS   <- c("IS","Day1","Day3","Day7")
SAMPLES  <- list(NC = c("NCs1","NCs2","NCs3"), IS = c("Is1","Is2","Is3"),
                 Day1 = c("Day1s1","Day1s2","Day1s3"), Day3 = c("Day3s1","Day3s2","Day3s3"),
                 Day7 = c("Day7s1","Day7s2","Day7s3"))
TP_COLORS <- c("IS"="#8A8A8A","Day1"="#D6404C","Day3"="#E8753A","Day7"="#2E7D6F")

tpm[, NC_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = SAMPLES$NC]
for (g in GROUPS) {
  tpm[, paste0(g, "_mean") := rowMeans(.SD, na.rm = TRUE), .SDcols = SAMPLES[[g]]]
  tpm[, paste0(g, "_log2FC") := log2((get(paste0(g, "_mean")) + 0.1) / (NC_mean + 0.1))]
}

enzyme_genes <- data.table(
  Gene = c("Mettl3","Wtap","Hnrnpa2b1","Ythdf2","Fto","Alkbh5",
           "Fbl","Nop56","Dkc1","Tet2"),
  Category = c(rep("m6A Writer/Reader/Eraser", 6), rep("2-O-Me/PSU/Bridge", 4)))

expr <- tpm[toupper(Gene) %in% toupper(enzyme_genes$Gene)]
expr[, Gene_up := toupper(Gene)]
enzyme_genes[, Gene_up := toupper(Gene)]
em <- merge(enzyme_genes,
            expr[, .(Gene_up, NC_mean,
                     IS_log2FC, Day1_log2FC, Day3_log2FC, Day7_log2FC)],
            by = "Gene_up", all.x = TRUE)
em <- unique(em[!is.na(Day1_log2FC)], by = "Gene")
em <- em[order(Day1_log2FC)]
em[, Gene := factor(Gene, levels = Gene)]

long <- melt(em, id.vars = c("Gene","Category"),
             measure.vars = paste0(GROUPS, "_log2FC"),
             variable.name = "Timepoint", value.name = "log2FC")
long[, Timepoint := factor(gsub("_log2FC","",Timepoint), levels = GROUPS)]
cat("=== Enzyme log2FC matrix (vs NC) ===\n")
print(dcast(long, Gene ~ Timepoint, value.var = "log2FC"))

save_n7 <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(OUT_DIR, paste0("Fig8", name, ".pdf")), p,
         width = w_mm/25.4, height = h_mm/25.4, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0("Fig8", name, ".png")), p,
         width = w_mm/25.4, height = h_mm/25.4, dpi = 600)
}

theme_fig8 <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.5, color = "grey30"),
        strip.text = element_text(size = 6.5, face = "bold"),
        strip.background = element_blank())

# ── PANEL B: enzyme × timepoint heatmap ──
p_b <- ggplot(long, aes(x = Timepoint, y = Gene, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.2f", log2FC)), size = 2.4, family = "Arial") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3D7EA6", mid = "white", high = "#D6404C",
                       midpoint = 0, limits = c(-1.2, 1.2), oob = scales::squish,
                       name = "log2FC\n(vs Sham)") +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = NULL, y = NULL,
       title = expression(bold("b")~~"qRT-PCR: Modification Enzyme Expression — Full Time Course"),
       subtitle = sprintf("DRS TPM-derived log2FC vs Sham (n=3/group); %d enzymes, ordered by Day1 response", nrow(em))) +
  theme_fig8 +
  theme(legend.key.size = unit(8, "pt"), legend.title = element_text(size = 6),
        legend.text = element_text(size = 5.5),
        panel.spacing = unit(1.5, "mm"))
save_n7(p_b, "B_qRT_PCR_Enzymes", 110, 120)
cat("Panel B (all-group heatmap) saved\n")

# ── PANEL C: 5 WB targets × timepoint dodged bars ──
wb <- em[Gene_up %in% toupper(c("Mettl3","Hnrnpa2b1","Fbl","Nop56","Alkbh5"))]
wb <- unique(wb, by = "Gene_up")
wb_long <- melt(wb, id.vars = c("Gene","Category"),
                measure.vars = paste0(GROUPS, "_log2FC"),
                variable.name = "Timepoint", value.name = "log2FC")
wb_long[, Timepoint := factor(gsub("_log2FC","",Timepoint), levels = GROUPS)]

p_c <- ggplot(wb_long, aes(x = Gene, y = log2FC, fill = Timepoint)) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_text(aes(label = sprintf("%+.2f", log2FC), group = Timepoint),
            position = position_dodge(0.75),
            vjust = ifelse(wb_long$log2FC > 0, -0.4, 1.4), size = 1.9, family = "Arial") +
  scale_fill_manual(values = TP_COLORS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  labs(x = NULL, y = "log2FC (vs Sham)",
       title = expression(bold("c")~~"Western Blot: Key Modification Proteins — Full Time Course"),
       subtitle = sprintf("%d proteins, DRS TPM-derived; m6A machinery down, 2'-O-Me writers up (Day1), trajectories to Day7", nrow(wb))) +
  theme_fig8 + theme(legend.position = "top", legend.margin = margin(b = -2))
save_n7(p_c, "C_Western_Blot_Targets", 150, 100)
cat("Panel C (all-group dodged bars) saved\n")

fwrite(dcast(long, Gene + Category ~ Timepoint, value.var = "log2FC"),
       file.path(OUT_DIR, "Table_Fig8BC_enzyme_log2FC_allgroups.csv"))
cat("=== B/C all-group extension done ===\n")
