# =============================================================================
# A2 step-1 (gene level): REDIportal mm10 x inosine validated-set overlap
# 2026-10-22 — orthogonal external validation of the inosine verdict arm.
#   Our set  : INOSINE_UP validated genes (HC50, D1-anchor gain; 245 genes,
#              fig3_8mod_spatial/data/upload_5mod_validated/INOSINE_UP.txt)
#   External : REDIportal mouse atlas TABLE1_mm10_v3 (107,094 A-to-I sites,
#              9,358 edited genes, GENCODE VM16 symbols) — downloaded
#              2026-10-22 from rediportal.cloud.ba.infn.it/download/
#   Universe : bg_genes.txt (10,841 DRS-detected background genes)
#   Controls : INOSINE_DOWN (same-mod direction control), M6A_UP / PSU_UP
#              (cross-modification controls), 1,000 random draws (empirical P)
# A2 plan ref: plans/A2_orthogonal_validation_plan.md §1.2 (site-level
# overlap + A->G mismatch signature = step 2, pending E: drive).
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "../data"
out_dir  <- file.path(data_dir, "rediportal_20261022")

# ── Load sets ────────────────────────────────────────────────────────────────
rd_lines <- function(p) { g <- unique(readLines(p)); g[!g %in% c("", "-")] }
ino_up   <- rd_lines(file.path(data_dir, "upload_5mod_validated/INOSINE_UP.txt"))
ino_dn   <- rd_lines(file.path(data_dir, "upload_5mod_validated/INOSINE_DOWN.txt"))
m6a_up   <- rd_lines(file.path(data_dir, "upload_5mod_validated/M6A_UP.txt"))
psu_up   <- rd_lines(file.path(data_dir, "upload_5mod_validated/PSU_UP.txt"))
m5c_dn   <- rd_lines(file.path(data_dir, "upload_5mod_validated/M5C_DOWN.txt"))
universe <- rd_lines(file.path(data_dir, "bg_genes.txt"))

rp <- fread(file.path(out_dir, "TABLE1_mm10_v3.txt.gz"))
rp_genes_all <- unique(rp$Gene.wgEncodeGencodeBasicVM16)
rp_genes_all <- rp_genes_all[!is.na(rp_genes_all) & rp_genes_all != "-"]
# Stringent subset: genes with >=1 transcript-body site (exonic/UTR3/UTR5/splicing)
body_func <- c("exonic", "UTR3", "UTR5", "splicing")
rp_genes_body <- unique(rp[Func.wgEncodeGencodeBasicVM16 %in% body_func]$Gene.wgEncodeGencodeBasicVM16)
rp_genes_body <- rp_genes_body[!is.na(rp_genes_body) & rp_genes_body != "-"]

cat(sprintf("REDIportal edited genes: %d (all) | %d (transcript-body sites)\n",
            length(rp_genes_all), length(rp_genes_body)))
cat(sprintf("Universe: %d | INOSINE_UP: %d | INOSINE_DOWN: %d\n",
            length(universe), length(ino_up), length(ino_dn)))

# ── Overlap engine ───────────────────────────────────────────────────────────
overlap_test <- function(set, rp_genes, universe, label, n_draw = 1000) {
  set_u <- intersect(set, universe)
  rp_u  <- intersect(rp_genes, universe)
  a <- length(intersect(set_u, rp_u))
  b <- length(set_u) - a
  c <- length(rp_u) - a
  d <- length(universe) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  # empirical P: random draws of |set_u| from universe
  set.seed(42)
  draws <- replicate(n_draw, length(intersect(sample(universe, length(set_u)), rp_u)))
  emp_p <- (sum(draws >= a) + 1) / (n_draw + 1)
  res <- data.table(
    set = label, set_n = length(set), set_in_universe = length(set_u),
    rp_in_universe = length(rp_u), overlap = a,
    pct_of_set = round(100 * a / length(set_u), 1),
    expected_by_chance = round(length(set_u) * length(rp_u) / length(universe), 1),
    fold_enrichment = round(a / (length(set_u) * length(rp_u) / length(universe)), 2),
    fisher_OR = round(unname(ft$estimate), 2), fisher_P = signif(ft$p.value, 3),
    empirical_P = signif(emp_p, 3)
  )
  list(res = res, hits = intersect(set_u, rp_u))
}

run_panel <- function(rp_genes, panel_label) {
  cat(sprintf("\n== Panel: %s ==\n", panel_label))
  out <- rbindlist(list(
    overlap_test(ino_up, rp_genes, universe, "INOSINE_UP (verdict arm)")$res,
    overlap_test(ino_dn, rp_genes, universe, "INOSINE_DOWN (direction control)")$res,
    overlap_test(m6a_up, rp_genes, universe, "M6A_UP (cross-mod control)")$res,
    overlap_test(psu_up, rp_genes, universe, "PSU_UP (cross-mod control)")$res,
    overlap_test(m5c_dn, rp_genes, universe, "M5C_DOWN (cross-mod control)")$res
  ))
  print(out)
  out
}

res_all  <- run_panel(rp_genes_all,  "REDIportal all edited genes")
res_body <- run_panel(rp_genes_body, "REDIportal transcript-body edited genes")

hits <- overlap_test(ino_up, rp_genes_all, universe, "INOSINE_UP")$hits
cat(sprintf("\nINOSINE_UP ∩ REDIportal: %d genes\n", length(hits)))
cat(paste(sort(hits), collapse = ", "), "\n")

# canonical edited-gene anecdotes (well-characterized mammalian targets)
anecdotes <- c("Blcap", "Cyfip2", "Copa", "Cog3", "Flnb", "Flna", "Ncoa3",
               "Gabra3", "Gria2", "Htr2c", "Kcna1", "Neil1", "Sgip1", "Zdhhc6")
cat("\nAnecdote check (in INOSINE_UP / in REDIportal):\n")
for (g in anecdotes) {
  cat(sprintf("  %-8s inosine_up=%s  rediportal=%s\n", g,
              g %in% ino_up, g %in% rp_genes_all))
}

fwrite(rbind(res_all, res_body), file.path(out_dir, "REDIportal_overlap_stats.tsv"), sep = "\t")
fwrite(data.table(gene = sort(hits)), file.path(out_dir, "REDIportal_overlap_genes_INOSINE_UP.tsv"), sep = "\t")
cat("\nSaved: REDIportal_overlap_stats.tsv | REDIportal_overlap_genes_INOSINE_UP.tsv\n")
