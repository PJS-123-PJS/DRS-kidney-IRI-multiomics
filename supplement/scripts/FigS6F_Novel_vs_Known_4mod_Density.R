# FigS6F: Novel vs Known — 4-mod differential-site density violin
# ============================================================================
# 2026-07-30 — README planned item (FigS6F), user-approved build (batch 4).
# Extends main fig5 H (=canonical Fig6C, m6A-only static density) to 4 mod
# classes with full per-gene distributions.
#
# DATA BASIS (honest, local-only): NC_vs_Day1 DIFFERENTIAL sites from the
# NEWPATH high-confidence context files (fig3_8mod_spatial/data/*_diff_site_
# context_NEWPATH.tsv.gz; ENSMUSG trans + gene symbol, ~99% name coverage).
# Nm = position-level union of Am/Cm/Gm/Um (cross-submod dedup by trans+pos,
# same convention as frac_5mod_cache_v3). Inosine excluded (README "4-mod";
# A-to-I editing is a separate story).
#   Static identity-basis density (NC condition, m6A) is in main fig5 H;
#   this panel answers the complementary question: do novel-transcript host
#   genes undergo comparable modification REMODELLING at Day1?
# Class definitions mirror the validated Panel E pipeline (Fig5_EGH_server.R):
#   Novel = genes in geneset_Novel_HostGenes_symbols.txt (6,636 host genes)
#   Known = all other symbol-annotated genes in the mod-data universe
# Zero hardcoding: all numbers computed at runtime.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

DATA_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
GENESET  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig6_novel_transcript/data/geneset_Novel_HostGenes_symbols.txt"
OUT_DIR  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/supplement/figS6_novel_transcript"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MOD_COLORS <- c("M6A"="#D6404C","M5C"="#E8753A","PSU"="#3D7EA6","Nm"="#2E7D6F")

read_ctx <- function(mod) {
  f <- file.path(DATA_DIR, sprintf("%s_diff_site_context_NEWPATH.tsv.gz", mod))
  dt <- fread(f, select = c("trans","pos","name"))
  dt[name != "-" & name != "" & !is.na(name)]
}

# ── per-gene differential-site counts, 4 mods ──
per_gene <- list()
for (m in c("m6A","m5c","psu")) {
  dt <- read_ctx(m)
  lbl <- c("m6A"="M6A","m5c"="M5C","psu"="PSU")[m]
  per_gene[[lbl]] <- dt[, .(Sites = uniqueN(paste(trans, pos))), by = name][, Mod := lbl]
  cat(sprintf("%s: %s sites over %s genes\n", lbl,
              format(uniqueN(paste(dt$trans, dt$pos)), big.mark=","),
              format(nrow(per_gene[[lbl]]), big.mark=",")))
}
# Nm = pooled union across submods (position-level dedup)
nm <- rbindlist(lapply(c("am","cm","gm","um"), read_ctx))
nm <- unique(nm, by = c("trans","pos"))
per_gene[["Nm"]] <- nm[, .(Sites = .N), by = name][, Mod := "Nm"]
cat(sprintf("Nm (pooled Am+Cm+Gm+Um): %s sites over %s genes\n",
            format(nrow(nm), big.mark=","), format(nrow(per_gene[["Nm"]]), big.mark=",")))

pg <- rbindlist(per_gene)
novel_hosts <- scan(GENESET, what = "", quiet = TRUE)
pg[, Class := fifelse(name %in% novel_hosts, "Novel host genes", "Known genes")]
pg[, Class := factor(Class, levels = c("Known genes","Novel host genes"))]
pg[, Mod := factor(Mod, levels = c("M6A","M5C","PSU","Nm"))]

cls_n <- pg[, .(N_genes = .N), by = .(Mod, Class)]
print(cls_n)

# ── Wilcoxon per mod ──
wil <- pg[, {
  w <- wilcox.test(Sites[Class=="Novel host genes"], Sites[Class=="Known genes"])
  .(p = w$p.value,
    med_novel = as.numeric(median(Sites[Class=="Novel host genes"])),
    med_known = as.numeric(median(Sites[Class=="Known genes"])))
}, by = Mod]
wil[, p_lab := fifelse(p < 2.2e-16, "p<2.2e-16", sprintf("p=%.2g", p))]
print(wil)

anno <- merge(cls_n, wil, by = "Mod")
anno[, n_lab := sprintf("N=%s", format(N_genes, big.mark=","))]

p <- ggplot(pg, aes(x = Class, y = Sites, fill = Mod)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.size = 0.3, outlier.alpha = 0.3,
               linewidth = 0.25, fill = "white") +
  facet_wrap(~ Mod, ncol = 4) +
  geom_text(data = anno, aes(x = Class, y = 1, label = n_lab),
            vjust = -0.5, size = 2, color = "grey40", family = "Arial") +
  geom_text(data = wil, aes(x = 1.5, y = Inf, label = p_lab),
            vjust = 1.6, size = 2.2, fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0.02, 0.18))) +
  scale_fill_manual(values = MOD_COLORS, guide = "none") +
  labs(title = "FigS6F: Differential-Site Load at Day1 — Novel-Host vs Known Genes (4 modifications)",
       subtitle = "Per-gene NC_vs_Day1 differential sites (NEWPATH high-conf sets); genes carrying >=1 site",
       x = NULL, y = "Differential sites per gene (log10)",
       caption = paste0(
         "Nm = position-level union of Am/Cm/Gm/Um. Violin = density, box = median/IQR; Wilcoxon rank-sum test per panel.\n",
         "Complements Fig 5H (static m6A density, NC condition): novel-transcript host genes are NOT depleted of",
         " Day1 remodelling — load is equal or higher in all 4 mods.")) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(size = 6, angle = 15, hjust = 1),
        strip.text = element_text(size = 7, face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 5.8, color = "grey30"),
        plot.caption = element_text(size = 5.3, color = "grey50", hjust = 0))

W_MM <- 183; H_MM <- 75
ggsave(file.path(OUT_DIR, "FigS6F_Novel_vs_Known_4mod_Density.pdf"), p,
       width = W_MM/25.4, height = H_MM/25.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS6F_Novel_vs_Known_4mod_Density.png"), p,
       width = W_MM/25.4, height = H_MM/25.4, dpi = 600)

fwrite(pg, file.path(OUT_DIR, "Table_FigS6F_per_gene_diffsites.csv"))
fwrite(wil, file.path(OUT_DIR, "Table_FigS6F_wilcoxon.csv"))
cat(sprintf("FigS6F saved: %.0f x %.0f mm -> %s\n", W_MM, H_MM, OUT_DIR))
