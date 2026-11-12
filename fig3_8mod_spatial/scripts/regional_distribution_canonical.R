# regional_distribution_canonical.R — per-mod region% + 3'UTR preference test
# ============================================================================
# 2026-07-23 — canonical replacement for Table_8mod_regional_distribution.csv
# (Jun-17, old trans_diff basis, m6a row tagged "Down (746)" = falsified set).
# Method: site_diff differential filter (same as Fig3A cache v2) joined to
# all_<mod>.rds structure via mult="first" (memory-safe: no pre-dedup copy).
# Outputs:
#   data/Table_regional_distribution_canonical.csv  (8 mods, region %)
#   prints Wilcoxon p for base-mod vs 2'-O-Me 3'UTR preference index
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })

SITE_DIR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
DATA_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"

MODS <- list(m6a="m6a", m5c="m5c", psu="psu", inosine="inosine",
             Am="Am", Cm="Cm", Gm="Gm", Um="Um")
TYPE <- c(m6a="Base", m5c="Base", psu="Base", inosine="Base",
          Am="2'-O-Me", Cm="2'-O-Me", Gm="2'-O-Me", Um="2'-O-Me")

rows <- list()
for (m in names(MODS)) {
  cat(sprintf("\n=== %s ===\n", m))
  r <- readRDS(file.path(SITE_DIR, MODS[[m]], paste0("all_", m, ".rds")))
  setDT(r)
  r[, c("Gene","Depth","Fraction","Chr","Motif","sample","group") := NULL]
  gc(verbose = FALSE)
  setkey(r, Trans, Pos)          # in-place ordering, no dedup copy
  for (cmp in c("NC_vs_Day1","NC_vs_Day3","NC_vs_Day7")) {
    f <- file.path(SITE_DIR, MODS[[m]], "site_diff", cmp,
                   paste0(cmp, ".", m, ".per_site.add_name.xls.gz"))
    if (!file.exists(f)) { cat("  missing:", cmp, "\n"); next }
    dt <- fread(f, select = c("trans","pos","a_pct_modified","b_pct_modified",
                              "effect_size","map_pvalue"))
    dt <- dt[pmax(a_pct_modified, b_pct_modified) > 0.7 &
             abs(effect_size) >= 0.1 & map_pvalue < 0.05]
    if (!nrow(dt)) next
    dt2 <- r[dt, on = .(Trans = trans, Pos = pos), mult = "first", nomatch = 0]
    dt2 <- dt2[Gene_Struc %in% c("five_prime_UTR","CDS","three_prime_UTR")]
    tb <- table(factor(dt2$Gene_Struc, levels = c("five_prime_UTR","CDS","three_prime_UTR")))
    tot <- sum(tb)                              # zero-count regions -> 0, never NA
    cat(sprintf("  %s: n=%d | 5'UTR %.1f%% CDS %.1f%% 3'UTR %.1f%%\n", cmp, tot,
                100*tb["five_prime_UTR"]/tot, 100*tb["CDS"]/tot, 100*tb["three_prime_UTR"]/tot))
    rows[[paste(m, cmp)]] <- data.table(
      mod = m, type = TYPE[[m]], comparison = cmp, n_sites = tot,
      utr5_pct = 100*tb["five_prime_UTR"]/tot,
      cds_pct  = 100*tb["CDS"]/tot,
      utr3_pct = 100*tb["three_prime_UTR"]/tot)
    rm(dt, dt2); gc(verbose = FALSE)
  }
  rm(r); gc(verbose = FALSE)
}

out <- rbindlist(rows)
# pooled across comparisons (site-count weighted)
pooled <- out[, .(n_sites = sum(n_sites),
                  utr5_pct = sum(utr5_pct*n_sites)/sum(n_sites),
                  cds_pct  = sum(cds_pct*n_sites)/sum(n_sites),
                  utr3_pct = sum(utr3_pct*n_sites)/sum(n_sites)), by = .(mod, type)]
pooled[, pref_index := utr3_pct / utr5_pct]
fwrite(out,    file.path(DATA_DIR, "Table_regional_distribution_canonical_bycmp.csv"))
fwrite(pooled, file.path(DATA_DIR, "Table_regional_distribution_canonical.csv"))

cat("\n=== POOLED (canonical) ===\n"); print(pooled)
b <- pooled[type == "Base", pref_index]; n <- pooled[type == "2'-O-Me", pref_index]
w <- suppressWarnings(wilcox.test(b, n))
cat(sprintf("\n3'UTR preference index: Base %.2f-%.2f vs 2'-O-Me %.2f-%.2f | Wilcoxon p = %.4f\n",
            min(b), max(b), min(n), max(n), w$p.value))
cat(sprintf("Base 3'UTR%% range: %.1f-%.1f | 2'-O-Me CDS%% range: %.1f-%.1f\n",
            min(pooled[type=="Base", utr3_pct]), max(pooled[type=="Base", utr3_pct]),
            min(pooled[type=="2'-O-Me", cds_pct]), max(pooled[type=="2'-O-Me", cds_pct])))
