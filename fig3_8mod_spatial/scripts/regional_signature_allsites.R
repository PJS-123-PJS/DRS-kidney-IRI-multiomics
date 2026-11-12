# regional_signature_allsites.R — ALL-SITES regional signature per modification
# ============================================================================
# 2026-07-23 — §4.1's "regional coding" identity claims need the all-sites
# basis (n=millions, stable), NOT the differential-site subset (Nm n=1-34,
# underpowered). all_*.rds Gene_Struc is the company's own annotation, self-
# consistent with pos coordinates (same source as cache v2 join).
# Memory-safe: in-place slim to Gene_Struc only (7GB box lessons 07-22/23).
# Outputs:
#   data/Table_regional_signature_allsites.csv   (8 mods, region %, n)
#   prints Wilcoxon base vs 2'-O-Me 3'UTR-preference index
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })

SITE_DIR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
DATA_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"

MODS <- c("m6a","m5c","psu","inosine","Am","Cm","Gm","Um")
TYPE <- c(m6a="Base", m5c="Base", psu="Base", inosine="Base",
          Am="2'-O-Me", Cm="2'-O-Me", Gm="2'-O-Me", Um="2'-O-Me")

rows <- list()
for (m in MODS) {
  cat(sprintf("=== %s ===\n", m))
  r <- readRDS(file.path(SITE_DIR, m, paste0("all_", m, ".rds"))); setDT(r)
  keep <- c("Trans","Pos","Gene_Struc")
  drop <- setdiff(names(r), keep)
  r[, (drop) := NULL]                       # in-place, no copy
  r <- unique(r, by = c("Trans","Pos"))     # site-level dedup (small key cols only)
  tb <- table(factor(r$Gene_Struc, levels = c("five_prime_UTR","CDS","three_prime_UTR")))
  other <- sum(!r$Gene_Struc %in% c("five_prime_UTR","CDS","three_prime_UTR"))
  tot <- sum(tb)
  cat(sprintf("  sites=%s (excl. non-3region: %s) | 5'UTR %.1f%% CDS %.1f%% 3'UTR %.1f%%\n",
              format(tot, big.mark=","), format(other, big.mark=","),
              100*tb[1]/tot, 100*tb[2]/tot, 100*tb[3]/tot))
  rows[[m]] <- data.table(mod = m, type = TYPE[[m]], n_sites = tot,
                          utr5_pct = as.numeric(100*tb[1]/tot),
                          cds_pct  = as.numeric(100*tb[2]/tot),
                          utr3_pct = as.numeric(100*tb[3]/tot))
  rm(r); gc(verbose = FALSE)
}

out <- rbindlist(rows)
out[, pref_index := utr3_pct / utr5_pct]
fwrite(out, file.path(DATA_DIR, "Table_regional_signature_allsites.csv"))

cat("\n=== ALL-SITES regional signature ===\n"); print(out)
b <- out[type == "Base", pref_index]; n <- out[type == "2'-O-Me", pref_index]
w <- suppressWarnings(wilcox.test(b, n))
cat(sprintf("\n3'UTR preference index (all sites): Base %.2f-%.2f vs 2'-O-Me %.2f-%.2f | Wilcoxon p = %.4f\n",
            min(b), max(b), min(n), max(n), w$p.value))
