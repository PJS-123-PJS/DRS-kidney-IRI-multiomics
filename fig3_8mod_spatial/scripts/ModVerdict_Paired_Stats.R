# =============================================================================
# 5-mod paired-arbitration stats with CIs (replaces all "wave" claims)
# 2026-07-22 — per-mod caller-independent verdicts:
#   m6A | m5C | PSU | inosine | Nm (Um/Am/Cm/Gm reported per base)
# Both anchors per caller: D1-HC50/70 (mean_frac_Day1>=t) and SH-HC50/70.
# wilcox.test(conf.int=TRUE) gives pseudomedian CI.
# Output: data/modverdict_paired_stats.tsv  (figure script reads this)
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
nc_samples <- c("NCs1", "NCs2", "NCs3")
d1_samples <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc_samples, d1_samples)

run_test <- function(dt, label) {
  if (nrow(dt) < 10) return(data.table(set = label, n = nrow(dt)))
  d <- dt$mean_frac_Day1 - dt$mean_frac_Sham
  wt <- suppressWarnings(wilcox.test(dt$mean_frac_Day1, dt$mean_frac_Sham,
                                   paired = TRUE, conf.int = TRUE))
  data.table(set = label, n = nrow(dt),
             median_Sham = median(dt$mean_frac_Sham),
             median_Day1 = median(dt$mean_frac_Day1),
             median_delta_pp = 100 * median(d),
             est_pp = 100 * wt$estimate,          # pseudomedian of differences
             ci_lo_pp = 100 * wt$conf.int[1],
             ci_hi_pp = 100 * wt$conf.int[2],
             pct_up = 100 * mean(d > 0),
             pct_up_gt10pp = 100 * mean(d > 0.10),
             pct_down_gt10pp = 100 * mean(d < -0.10),
             wilcox_p = wt$p.value)
}

paired_from_ctx <- function(site_ctx, gate_motifs = NULL) {
  if (!is.null(gate_motifs)) site_ctx[, is_gate := motif %chin% gate_motifs]
  else site_ctx[, is_gate := FALSE]
  site_ctx[n_det_NC == 3 & n_det_Day1 == 3]
}

anchors <- function(paired, tag) {
  rbind(
    run_test(paired[mean_frac_Day1 >= 0.7], paste0(tag, "|D1-HC70")),
    run_test(paired[mean_frac_Day1 >= 0.5], paste0(tag, "|D1-HC50")),
    run_test(paired[mean_frac_Sham >= 0.7], paste0(tag, "|SH-HC70")),
    run_test(paired[mean_frac_Sham >= 0.5], paste0(tag, "|SH-HC50")),
    fill = TRUE)
}

out <- list()

# ---- m6A: saved site_ctx with is_drach ----
cat("== m6A ==\n")
sc <- fread(file.path(data_dir, "site_ctx_nc_day1_NEWPATH.tsv.gz"))
paired <- sc[n_det_NC == 3 & n_det_Day1 == 3]
out[["m6A"]] <- anchors(paired, "m6A")
out[["m6A-DRACH"]] <- rbind(
  run_test(paired[mean_frac_Sham >= 0.7 & is_drach], "m6A DRACH+|SH-HC70"),
  run_test(paired[mean_frac_Day1 >= 0.7 & is_drach], "m6A DRACH+|D1-HC70"),
  fill = TRUE)
rm(sc, paired); gc(verbose = FALSE)

# ---- m5C: saved site_ctx with motif gate ----
cat("== m5C ==\n")
me <- fread(file.path(data_dir, "m5c_motif_enrichment_NEWPATH.tsv"))
gate5 <- me[enrich >= 3 & n_hc70 >= 100]$motif
sc <- fread(file.path(data_dir, "site_ctx_nc_day1_m5c_NEWPATH.tsv.gz"))
sc[, is_gate := motif %chin% gate5]
paired <- sc[n_det_NC == 3 & n_det_Day1 == 3]
out[["m5C"]] <- anchors(paired, "m5C")
out[["m5C-gate"]] <- rbind(
  run_test(paired[mean_frac_Sham >= 0.7 & is_gate], "m5C gate+|SH-HC70"),
  run_test(paired[mean_frac_Day1 >= 0.7 & is_gate], "m5C gate+|D1-HC70"),
  fill = TRUE)
rm(sc, paired); gc(verbose = FALSE)

# ---- 6 other callers from ctx chunks ----
for (mod in c("psu", "inosine", "um", "am", "cm", "gm")) {
  cat("==", mod, "==\n")
  ctx_dir <- file.path(data_dir, paste0("site_ctx_chunks_", mod))
  ctx <- rbindlist(lapply(six, function(s) {
    dt <- fread(file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")))
    dt[, sample := s]
    dt[, grp := ifelse(s %chin% nc_samples, "Sham", "Day1")]
    dt
  }))
  site_ctx <- ctx[, .(
    n_det_NC = sum(grp == "Sham"),
    n_det_Day1 = sum(grp == "Day1"),
    mean_frac_Sham = mean(fraction[grp == "Sham"]),
    mean_frac_Day1 = mean(fraction[grp == "Day1"])
  ), by = .(trans, pos)]
  rm(ctx); gc(verbose = FALSE)
  tag <- switch(mod, psu = "PSU", inosine = "inosine",
                um = "Nm-Um", am = "Nm-Am", cm = "Nm-Cm", gm = "Nm-Gm")
  paired <- site_ctx[n_det_NC == 3 & n_det_Day1 == 3]
  out[[mod]] <- anchors(paired, tag)
  rm(site_ctx, paired); gc(verbose = FALSE)
}

final <- rbindlist(out, fill = TRUE)
final[, c("mod", "anchor") := tstrsplit(set, "|", fixed = TRUE)]
setcolorder(final, c("mod", "anchor", "n", "median_Sham", "median_Day1",
                     "median_delta_pp", "est_pp", "ci_lo_pp", "ci_hi_pp",
                     "pct_up", "pct_up_gt10pp", "pct_down_gt10pp", "wilcox_p", "set"))
fwrite(final, file.path(data_dir, "modverdict_paired_stats.tsv"), sep = "\t")
print(final[, .(mod, anchor, n, median_delta_pp, ci_lo_pp, ci_hi_pp, pct_up, wilcox_p)], digits = 3)
cat("\nMODVERDICT STATS DONE\n")
