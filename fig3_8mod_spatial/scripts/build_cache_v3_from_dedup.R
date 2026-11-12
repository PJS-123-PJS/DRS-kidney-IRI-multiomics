# build_cache_v3_from_dedup.R — derive frac_5mod_cache_v3.rds
# ============================================================================
# 2026-07-29 — cache v3 = the fidelity-gated NC-deduplicated, true-Nm-pooled
# rows from build_Fig3B_dedup_rows.R (fig3b_5mod_dedup_rows.rds).
#   * fixes cache v2 silent bug: "Nm" was Um-only (all_frac key lacked submod
#     name -> each submod overwrote the previous; gate-verified Um-only)
#   * NC/Sham = union across NC_vs_Day1/3/7 (no cross-comparison duplicates;
#     v2 counted recurring sites up to 3x, inflating Sham x1.14-1.61)
#   * Day1/3/7 rows identical to v2 (fidelity gate passed)
# Schema: superset of v2 (adds Trans/Pos for future audits). Consumers:
# Fig3A_5mod_metagene.R, Fig3B_3UTR_timecourse.R (both switched 2026-07-29).
# Fig3E_v4_RealData_Rebuild.R still reads v2 — panel handled manually in
# typesetting by user; re-point only if Fig3E is ever re-rendered.
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })
OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"

res <- readRDS(file.path(OUT_DIR, "fig3b_5mod_dedup_rows.rds")); setDT(res)
saveRDS(res, file.path(OUT_DIR, "frac_5mod_cache_v3.rds"))

cat(sprintf("cache v3 written: %s rows\n", format(nrow(res), big.mark = ",")))
cat("\ncount matrix:\n")
print(dcast(res[, .N, by = .(Modification, Group)], Modification ~ Group, value.var = "N"))
cat("\nregion % by mod:\n")
for (m in levels(res$Modification)) {
  tb <- round(prop.table(table(res[Modification == m]$Gene_Struc)) * 100, 1)
  cat(sprintf("  %-8s 5'UTR=%.1f%% CDS=%.1f%% 3'UTR=%.1f%%\n",
              m, tb["five_prime_UTR"], tb["CDS"], tb["three_prime_UTR"]))
}
