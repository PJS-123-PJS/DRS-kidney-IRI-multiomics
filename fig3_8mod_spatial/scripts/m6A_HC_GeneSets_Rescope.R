# =============================================================================
# Re-scope m6A gene sets to the REAL high-confidence signal (new canonical path)
# 2026-07-22 — replaces the noise-dominated 530-gene m6A-down set used in
# Fig2J / Fig3J / Fig4 enrichments.
#
# Inputs (all real, produced by m6A_NewPath_Reliability.R):
#   data/m6A_diff_site_context_NEWPATH.tsv.gz   matched site_diff + is_drach
#   data/geneset_m6A↓.txt                       legacy 530-gene set (for overlap)
# Outputs:
#   data/geneset_m6A_UP_HC50_NEWPATH.txt        genes with >=1 UP HC50 DRACH+ site
#   data/geneset_m6A_UP_HC70_NEWPATH.txt
#   data/geneset_m6A_DOWN_HC50_NEWPATH.txt      (kept for honesty; noise-rich)
#   data/m6A_geneset_rescope_summary.tsv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "../data"
dj <- fread(file.path(data_dir, "m6A_diff_site_context_NEWPATH.tsv.gz"))
cat(sprintf("matched diff sites: %s | genes: %s\n",
            format(nrow(dj), big.mark = ","),
            format(length(unique(dj$name)), big.mark = ",")))

mk <- function(d, id) {
  g <- sort(unique(d$name))
  data.table(set = id, n_sites = nrow(d), n_genes = length(g),
             n_genes_drach = length(unique(d[is_drach == TRUE]$name)))
}
summ <- rbind(
  mk(dj[diffType == "Up"],                                  "UP_all"),
  mk(dj[diffType == "Up"   & b_pct >= 0.5],                 "UP_HC50"),
  mk(dj[diffType == "Up"   & b_pct >= 0.7],                 "UP_HC70"),
  mk(dj[diffType == "Up"   & b_pct >= 0.5 & is_drach],      "UP_HC50_DRACH"),
  mk(dj[diffType == "Down"],                                "DOWN_all"),
  mk(dj[diffType == "Down" & a_pct >= 0.5],                 "DOWN_HC50"),
  mk(dj[diffType == "Down" & a_pct >= 0.5 & is_drach],      "DOWN_HC50_DRACH"))
print(summ)
fwrite(summ, file.path(data_dir, "m6A_geneset_rescope_summary.tsv"), sep = "\t")

g_up50  <- sort(unique(dj[diffType == "Up"   & b_pct >= 0.5 & is_drach]$name))
g_up70  <- sort(unique(dj[diffType == "Up"   & b_pct >= 0.7 & is_drach]$name))
g_dn50  <- sort(unique(dj[diffType == "Down" & a_pct >= 0.5]$name))
writeLines(g_up50, file.path(data_dir, "geneset_m6A_UP_HC50_NEWPATH.txt"))
writeLines(g_up70, file.path(data_dir, "geneset_m6A_UP_HC70_NEWPATH.txt"))
writeLines(g_dn50, file.path(data_dir, "geneset_m6A_DOWN_HC50_NEWPATH.txt"))

# ---- overlap with the legacy 530-gene m6A-down set ---------------------------
legacy_file <- file.path(data_dir, "geneset_m6A↓.txt")
if (file.exists(legacy_file)) {
  legacy <- readLines(legacy_file)
  legacy <- trimws(legacy[nzchar(trimws(legacy))])
  tab <- data.table(
    comparison = c("legacy530 vs DOWN_all", "legacy530 vs DOWN_HC50",
                   "legacy530 vs UP_HC50_DRACH", "legacy530 vs UP_all"),
    n_overlap = c(length(intersect(legacy, unique(dj[diffType == "Down"]$name))),
                  length(intersect(legacy, g_dn50)),
                  length(intersect(legacy, g_up50)),
                  length(intersect(legacy, unique(dj[diffType == "Up"]$name)))),
    pct_of_legacy = 0)
  tab[, pct_of_legacy := round(100 * n_overlap / length(legacy), 1)]
  cat(sprintf("legacy set: %d genes\n", length(legacy)))
  print(tab)
  fwrite(tab, file.path(data_dir, "m6A_geneset_legacy_overlap.tsv"), sep = "\t")
}
cat("== DONE ==\n")
