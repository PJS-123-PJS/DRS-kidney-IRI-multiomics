# build_Fig3F_region_sets.R — region-stratified m6A-UP gene sets (validated basis)
# ============================================================================
# 2026-07-23 — for Fig3F v2 (replaces wave-era Fig3F_WaveRegion_CellType).
# Validated M6A_UP set (2,829 genes, upload_5mod_validated/M6A_UP.txt)
# stratified by WHERE their HC UP sites sit on the transcript (canonical
# all_m6a.rds Gene_Struc via Trans+Pos mult="first" join — memory-safe, no
# pre-dedup copy; 7GB-box lessons 07-22/23).
# Gene assignment: majority region of its HC UP sites; ties -> "mixed"
# (excluded from the three strata, count disclosed).
# Outputs: data/fig3f_region_sets/{m6A_UP_CDS,m6A_UP_3UTR,m6A_UP_5UTR}.txt
#          + data/fig3f_region_sets/_assignment_summary.csv
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
set.seed(20260723)

SITE <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/m6a"
DATA <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
OUTD <- file.path(DATA, "fig3f_region_sets"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

up <- fread(file.path(DATA, "upload_5mod_validated/M6A_UP.txt"), header = FALSE)$V1
cat(sprintf("Validated M6A_UP set: %d genes\n", length(up)))

# ---- HC UP sites in validated genes (same filter as cache v2) ----
dt <- fread(cmd = paste("zcat", file.path(SITE, "site_diff/NC_vs_Day1/NC_vs_Day1.m6a.per_site.add_name.xls.gz")),
            select = c("trans","name","pos","a_pct_modified","b_pct_modified","map_pvalue","effect_size"))
dt <- dt[name %in% up]
dt <- dt[pmax(a_pct_modified, b_pct_modified) > 0.7 &
         effect_size >= 0.1 & map_pvalue < 0.05]
cat(sprintf("HC UP sites in validated genes: %s (genes: %d)\n",
            format(nrow(dt), big.mark = ","), uniqueN(dt$name)))

# ---- canonical structure join ----
r <- readRDS(file.path(SITE, "all_m6a.rds")); setDT(r)
r[, c("Gene","Depth","Fraction","Chr","Motif","sample","group") := NULL]   # in-place
gc(verbose = FALSE)
setkey(r, Trans, Pos)
dt2 <- r[dt, on = .(Trans = trans, Pos = pos), mult = "first", nomatch = 0]
dt2 <- dt2[Gene_Struc %in% c("five_prime_UTR","CDS","three_prime_UTR")]
cat(sprintf("structure-matched: %s sites (%.1f%%)\n",
            format(nrow(dt2), big.mark = ","), 100 * nrow(dt2) / nrow(dt)))
rm(r); gc(verbose = FALSE)

# ---- majority region per gene ----
cnt <- dt2[, .N, by = .(name, Gene_Struc)]
cnt[, tot := sum(N), by = name]
cnt[, frac := N / tot]
maj <- cnt[, {
  top <- which.max(frac)
  tied <- sum(frac == max(frac)) > 1
  .(region = ifelse(tied, "mixed", Gene_Struc[top]), n_sites = tot[1], top_frac = frac[top])
}, by = name]

cat("\n=== assignment ===\n"); print(table(maj$region))
cat(sprintf("sites per gene: median %d, max %d\n", median(maj$n_sites), max(maj$n_sites)))

sets <- list(m6A_UP_CDS  = maj[region == "CDS", name],
             m6A_UP_3UTR = maj[region == "three_prime_UTR", name],
             m6A_UP_5UTR = maj[region == "five_prime_UTR", name])
for (nm in names(sets)) {
  writeLines(sets[[nm]], file.path(OUTD, paste0(nm, ".txt")))
  cat(sprintf("  %-12s: %d genes\n", nm, length(sets[[nm]])))
}
fwrite(maj, file.path(OUTD, "_assignment_summary.csv"))
cat("=== Fig3F region sets done ===\n")
