# FigS6F v2: Novel vs Known — 4-mod density violin, IDENTITY basis (NC condition)
# ============================================================================
# 2026-07-30 — replaces v1 (NEWPATH NC_vs_Day1 differential/response basis)
# with the static identity basis matching main fig5H (Fig6C): distinct
# NC-detected sites per gene. Requires /mnt/e (05.methy) mounted.
#
# Pipeline (memory-safe, one mod at a time, 7GB box):
#   A) trans→symbol map: site_diff NC_vs_Day1 per_site.add_name files
#      (full tested universe — NOT restricted to differential sites),
#      select trans+name only.
#   B) per mod: all_<mod>.rds → in-place slim (Trans/Pos/group) → NC filter →
#      unique (Trans,Pos) → per-transcript site counts. Nm = pooled
#      Am/Cm/Gm/Um with position-level union (cache-v3 convention).
#   C) per-gene = sum of per-transcript counts per symbol (Panel E precedent:
#      server m6a_gene_*_summary was per-gene site counts).
# Class definitions mirror Fig5_EGH_server.R:
#   Novel = geneset_Novel_HostGenes_symbols.txt (6,636 host genes)
#   Known = all other symbol-annotated genes in the mod-data universe
# Zero hardcoding: all numbers computed at runtime.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

SITE_DIR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
GENESET  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig6_novel_transcript/data/geneset_Novel_HostGenes_symbols.txt"
OUT_DIR  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/supplement/figS6_novel_transcript"
CACHE    <- "/mnt/d/desktop/DRS/DRS_TechShowcase/supplement/data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

stopifnot(dir.exists(SITE_DIR))
MOD_COLORS <- c("M6A"="#D6404C","M5C"="#E8753A","PSU"="#3D7EA6","Nm"="#2E7D6F")
BASE_MODS  <- c(M6A="m6a", M5C="m5c", PSU="psu")   # label -> dir
NM_SUBS    <- c("Am","Cm","Gm","Um")

# ── A) trans→symbol map from site_diff per_site files (full tested universe) ──
map_rds <- file.path(CACHE, "trans2symbol_map.rds")
if (file.exists(map_rds)) {
  t2s <- readRDS(map_rds)
} else {
  maps <- list()
  for (md in c("m6a","m5c","psu","inosine","Am","Cm","Gm","Um")) {
    f <- file.path(SITE_DIR, md, "site_diff", "NC_vs_Day1",
                   sprintf("NC_vs_Day1.%s.per_site.add_name.xls.gz", md))
    if (!file.exists(f)) { cat(sprintf("  map source missing: %s\n", f)); next }
    m <- fread(f, select = c("trans","name"))
    m <- unique(m[name != "-" & name != "" & !is.na(name)])
    maps[[md]] <- m
    cat(sprintf("  map %s: %s trans rows\n", md, format(nrow(m), big.mark=",")))
  }
  t2s <- unique(rbindlist(maps), by = "trans")   # keep-first across mods
  saveRDS(t2s, map_rds)
}
cat(sprintf("trans→symbol map: %s transcripts\n", format(nrow(t2s), big.mark=",")))
sym_of <- setNames(t2s$name, t2s$trans)

# ── B) per-mod NC site sets / per-transcript counts (cached) ──
nc_siteset <- function(md) {
  cf <- file.path(CACHE, sprintf("nc_siteset_%s.rds", md))
  if (file.exists(cf)) return(readRDS(cf))
  r <- readRDS(file.path(SITE_DIR, md, sprintf("all_%s.rds", md))); setDT(r)
  keep <- intersect(names(r), c("Trans","Pos","group","sample"))
  r[, setdiff(names(r), keep) := NULL]          # in-place slim
  gc(verbose = FALSE)
  # all_<mod>.rds quirk: group values are Day1/Day3/Day7/IS/ctrl — Sham = "ctrl"
  # (NOT "NC"). Filter by sample names to be explicit and version-proof.
  r <- r[sample %in% c("NCs1","NCs2","NCs3")]
  r <- unique(r[, .(Trans, Pos)])
  cat(sprintf("  %s: %s NC sites over %s transcripts\n", md,
              format(nrow(r), big.mark=","), format(uniqueN(r$Trans), big.mark=",")))
  saveRDS(r, cf)
  r
}

tx <- list()
for (lbl in names(BASE_MODS)) {
  s <- nc_siteset(BASE_MODS[[lbl]])
  tx[[lbl]] <- s[, .(Sites = .N), by = Trans]
  rm(s); gc(verbose = FALSE)
}
# Nm = position-level union across submods (cache-v3 convention), then count
nm_set <- unique(rbindlist(lapply(NM_SUBS, nc_siteset)), by = c("Trans","Pos"))
tx[["Nm"]] <- nm_set[, .(Sites = .N), by = Trans]
rm(nm_set); gc(verbose = FALSE)

# ── C) per-gene (symbol) counts + class ──
pg <- rbindlist(lapply(names(tx), function(lbl) {
  d <- tx[[lbl]]
  d[, Symbol := sym_of[Trans]]
  d <- d[!is.na(Symbol)]
  d[, .(Sites = sum(Sites)), by = Symbol][, Mod := lbl]
}))
cat(sprintf("per-gene table: %s rows; symbol coverage by mod:\n", format(nrow(pg), big.mark=",")))
print(pg[, .(genes = .N), by = Mod])

novel_hosts <- scan(GENESET, what = "", quiet = TRUE)
pg[, Class := fifelse(Symbol %in% novel_hosts, "Novel host genes", "Known genes")]
pg[, Class := factor(Class, levels = c("Known genes","Novel host genes"))]
pg[, Mod := factor(Mod, levels = c("M6A","M5C","PSU","Nm"))]

cls_n <- pg[, .(N_genes = .N), by = .(Mod, Class)]
print(cls_n)

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
  labs(title = "FigS6F: Modification-Site Load — Novel-Host vs Known Genes (4 modifications)",
       subtitle = "Per-gene distinct sites detected in NC (Sham) samples, all transcript regions; genes carrying >=1 site",
       x = NULL, y = "Sites per gene (log10)",
       caption = paste0(
         "Nm = pooled Am/Cm/Gm/Um. Violin = density, box = median/IQR; Wilcoxon rank-sum test per panel.\n",
         "Identity basis (NC condition), same convention as Fig 5H: novel-transcript host genes are NOT",
         " modification-depleted.")) +
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

fwrite(pg, file.path(OUT_DIR, "Table_FigS6F_per_gene_NCsites.csv"))
fwrite(wil, file.path(OUT_DIR, "Table_FigS6F_wilcoxon.csv"))
cat(sprintf("FigS6F v2 (identity basis) saved: %.0f x %.0f mm -> %s\n", W_MM, H_MM, OUT_DIR))
