# Day1-anchored paired stoichiometry — the fair arbiter for UP-wave claims
# Reads existing site_ctx_chunks_<mod>/ctx_*.tsv.gz (6 NC/Day1 samples)
# For each mod: among fully-paired sites, anchor on Day1-HC (mean_frac_Day1>=t)
# and test whether Sham was lower (UP real) or equal (diff-call noise).
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
nc_samples <- c("NCs1", "NCs2", "NCs3")
d1_samples <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc_samples, d1_samples)

run_test <- function(dt, label) {
  if (nrow(dt) < 10) return(data.table(set = label, n = nrow(dt)))
  d <- dt$mean_frac_Day1 - dt$mean_frac_Sham
  wt <- suppressWarnings(wilcox.test(dt$mean_frac_Day1, dt$mean_frac_Sham, paired = TRUE))
  data.table(set = label, n = nrow(dt),
             median_Sham = median(dt$mean_frac_Sham),
             median_Day1 = median(dt$mean_frac_Day1),
             median_delta = median(d),
             pct_up = 100 * mean(d > 0),
             pct_up_gt10pp = 100 * mean(d > 0.10),
             pct_down_gt10pp = 100 * mean(d < -0.10),
             wilcox_p = wt$p.value)
}

out_all <- list()
for (mod in c("um", "psu", "inosine", "am", "cm", "gm")) {
  ctx_dir <- file.path(data_dir, paste0("site_ctx_chunks_", mod))
  if (!dir.exists(ctx_dir)) { cat("skip", mod, "(no ctx dir)\n"); next }
  # gate motifs learned by the battery (descriptive)
  me <- fread(file.path(data_dir, paste0(mod, "_motif_enrichment_NEWPATH.tsv")))
  gate <- me[enrich >= 3 & n_hc70 >= 100]$motif

  ctx <- rbindlist(lapply(six, function(s) {
    dt <- fread(file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")))
    dt[, sample := s]
    dt[, grp := ifelse(s %chin% nc_samples, "Sham", "Day1")]
    dt
  }))
  site_ctx <- ctx[, .(
    motif = motif[1],
    n_det_NC = sum(grp == "Sham"),
    n_det_Day1 = sum(grp == "Day1"),
    mean_frac_Sham = mean(fraction[grp == "Sham"]),
    mean_frac_Day1 = mean(fraction[grp == "Day1"])
  ), by = .(trans, pos)]
  site_ctx[, is_gate := motif %chin% gate]
  rm(ctx); gc(verbose = FALSE)

  paired <- site_ctx[n_det_NC == 3 & n_det_Day1 == 3]
  res <- rbind(
    run_test(paired[mean_frac_Day1 >= 0.7],            "D1-HC70 all"),
    run_test(paired[mean_frac_Day1 >= 0.7 & is_gate],  "D1-HC70 gate+"),
    run_test(paired[mean_frac_Day1 >= 0.5],            "D1-HC50 all"),
    run_test(paired[mean_frac_Day1 >= 0.5 & is_gate],  "D1-HC50 gate+"),
    # symmetric anchor used before, for direct comparison
    run_test(paired[mean_frac_Sham >= 0.5],            "SH-HC50 all"),
    fill = TRUE)
  res[, mod := mod]
  cat("\n=====", toupper(mod), "=====\n")
  print(res, digits = 3)
  out_all[[mod]] <- res
  rm(site_ctx, paired); gc(verbose = FALSE)
}
final <- rbindlist(out_all, fill = TRUE)
fwrite(final, file.path(data_dir, "day1anchored_paired_all_mods_NEWPATH.tsv"), sep = "\t")
cat("\nALL DAY1-ANCHORED DONE\n")
