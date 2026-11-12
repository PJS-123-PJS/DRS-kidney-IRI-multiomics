# =============================================================================
# m6A reliability battery on the CANONICAL new path — chunked v4 (7GB-safe)
# 2026-07-22
#
# v1/v2/v3 OOM-killed: rds (49.1M rows) loads at ~4.5GB, any subset copy kills
# the box. v3 proved depth>=20 share = 1.0000 and frac==0 share = 0.0000 (new
# path is pre-filtered, modified positions only). v4: immediately after load,
# compute is_drach and DROP motif/gene/chr/rel_location/group IN PLACE (no
# frame copy), then per-sample loop on the 6-column slim frame.
#
# Products (all _NEWPATH):
#   data/m6a_threshold_motif_counts_NEWPATH.tsv   (Part C, all 15 samples)
#   data/m6A_highconf_counts_NEWPATH.tsv          (Part D)
#   data/site_ctx_nc_day1_NEWPATH.tsv.gz          (per-site context, 6 samples)
#   data/m6A_highconf_tier_summary_NEWPATH.tsv    (Part A)
#   data/m6A_diff_site_context_NEWPATH.tsv.gz
#   data/m6A_highconf_star_genes_NEWPATH.tsv
#   data/m6A_paired_stoichiometry_NEWPATH.tsv     (Part B)
#   data/m6A_presence_drop_NEWPATH.tsv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

RDS <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/m6a/all_m6a.rds"
DIFF <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/m6a/site_diff/NC_vs_Day1/NC_vs_Day1.m6a.per_site.diff.add_name.xls"
data_dir <- "../data"
out_dir  <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
DRACH_RE <- "^[AGT][AG]AC[ACT]$"
nc_samples <- c("NCs1", "NCs2", "NCs3")
d1_samples <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc_samples, d1_samples)

cat("== loading all_m6a.rds (one big allocation) ==\n")
t0 <- Sys.time()
r <- readRDS(RDS)
setDT(r)
setnames(r, tolower(names(r)))
cat(sprintf("  %s rows in %.1f min | frac==0: %.4f | depth>=20: %.4f\n",
            format(nrow(r), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            mean(r$fraction == 0), mean(r$depth >= 20)))

# slim in place BEFORE any subsetting (7GB box: full frame ~4.5GB, copies kill)
r[, is_drach := grepl(DRACH_RE, motif)]
r[, motif := NULL]
r[, c("gene", "chr", "rel_location", "group") := NULL]
gc(verbose = FALSE)

all_samples <- unique(r$sample)
cat("  samples:", paste(all_samples, collapse = ", "), "\n")

# ---- chunked per-sample pass --------------------------------------------------
bin_out <- file.path(data_dir, "m6a_threshold_motif_counts_NEWPATH.tsv")
hc_out  <- file.path(data_dir, "m6A_highconf_counts_NEWPATH.tsv")
ctx_dir <- file.path(data_dir, "site_ctx_chunks")
dir.create(ctx_dir, showWarnings = FALSE, recursive = TRUE)
first <- TRUE
for (s in all_samples) {
  sub <- r[sample == s,
           .(trans, pos, fraction, gene_struc, is_drach)]
  grp <- ifelse(grepl("^NC", s), "Sham",
         ifelse(grepl("^I", s), "IS", sub("s[0-9]+$", "", s)))
  # Part C bins
  b <- sub[, .(N = .N, DRACH = sum(is_drach)),
           by = .(bin = pmin(9L, as.integer(floor(fraction * 10))))]
  b[, `:=`(sample = s, group = grp)]
  fwrite(b, bin_out, sep = "\t", append = !first)
  # Part D HC counts
  hc <- data.table(sample = s, group = grp,
                   N_hc50 = sub[fraction >= 0.5, .N],
                   N_hc70 = sub[fraction >= 0.7, .N])
  fwrite(hc, hc_out, sep = "\t", append = !first)
  # per-site chunk for Parts A/B (only NC/Day1 samples)
  if (s %chin% six) {
    fwrite(sub, file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")),
           sep = "\t", compress = "gzip")
  }
  cat(sprintf("  %s done (%s rows)\n", s, format(nrow(sub), big.mark = ",")))
  rm(sub, b, hc); gc(verbose = FALSE)
  first <- FALSE
}
rm(r); gc(verbose = FALSE)

# Part C pooled summary
bins <- fread(bin_out)
pooled <- bins[, .(N = sum(N), DRACH = sum(DRACH)), by = bin]
setorder(pooled, bin)
pooled[, `:=`(cumN = rev(cumsum(rev(N))), cumDRACH = rev(cumsum(rev(DRACH))))]
pooled[, pct := 100 * cumDRACH / cumN]
cat("== Part C: cumulative DRACH% by threshold (new path, all 15) ==\n")
print(pooled[, .(thr = bin / 10, cumN, pct = round(pct, 1))])

# ---- combine 6 per-site chunks into site context ------------------------------
cat("== building site context from chunks ==\n")
ctx_list <- lapply(six, function(s) {
  dt <- fread(file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")))
  dt[, sample := s]
  dt[, grp := ifelse(s %chin% nc_samples, "Sham", "Day1")]
  dt
})
ctx <- rbindlist(ctx_list); rm(ctx_list); gc(verbose = FALSE)
site_ctx <- ctx[, .(
  is_drach = any(is_drach),
  gene_struc = gene_struc[1],
  n_det_NC = sum(grp == "Sham"),
  n_det_Day1 = sum(grp == "Day1"),
  mean_frac_Sham = mean(fraction[grp == "Sham"]),
  mean_frac_Day1 = mean(fraction[grp == "Day1"])
), by = .(trans, pos)]
rm(ctx); gc(verbose = FALSE)
cat(sprintf("  unique sites: %s\n", format(nrow(site_ctx), big.mark = ",")))
fwrite(site_ctx, file.path(data_dir, "site_ctx_nc_day1_NEWPATH.tsv.gz"),
       sep = "\t", compress = "gzip")

# ============ Part A: site_diff tiers ==========================================
cat("== Part A: site_diff tiers ==\n")
diff_dt <- fread(DIFF)
setnames(diff_dt, c("trans", "name", "pos", "score", "a_pct", "b_pct",
                    "map_pvalue", "effect_size", "diffType"))
cat(sprintf("  diff sites: %s | Up=%s Down=%s\n",
            format(nrow(diff_dt), big.mark = ","),
            format(sum(diff_dt$diffType == "Up"), big.mark = ","),
            format(sum(diff_dt$diffType == "Down"), big.mark = ",")))
dj <- merge(diff_dt, site_ctx, by = c("trans", "pos"))
cat(sprintf("  matched %s / %s (%.1f%%)\n", format(nrow(dj), big.mark = ","),
            format(nrow(diff_dt), big.mark = ","), 100 * nrow(dj) / nrow(diff_dt)))
mk_tier <- function(d, id) {
  n <- nrow(d); k <- sum(d$is_drach)
  ci <- if (n > 0) binom.test(k, n)$conf.int else c(NA_real_, NA_real_)
  data.table(tier = id, N = n, DRACH = k, DRACH_pct = 100 * k / n,
             CI_lo = 100 * ci[1], CI_hi = 100 * ci[2])
}
tier_dt <- rbind(
  mk_tier(dj[diffType == "Down"],                "DOWN_all"),
  mk_tier(dj[diffType == "Down" & a_pct >= 0.5], "DOWN_HC50"),
  mk_tier(dj[diffType == "Down" & a_pct >= 0.7], "DOWN_HC70"),
  mk_tier(dj[diffType == "Up"],                  "UP_all"),
  mk_tier(dj[diffType == "Up"   & b_pct >= 0.5], "UP_HC50"),
  mk_tier(dj[diffType == "Up"   & b_pct >= 0.7], "UP_HC70"))
print(tier_dt)
fwrite(tier_dt, file.path(data_dir, "m6A_highconf_tier_summary_NEWPATH.tsv"), sep = "\t")
fwrite(dj[, .(trans, pos, name, diffType, a_pct, b_pct, effect_size,
              is_drach, gene_struc)],
       file.path(data_dir, "m6A_diff_site_context_NEWPATH.tsv.gz"), sep = "\t", compress = "gzip")

star <- c("Park7", "Fis1", "Gclm", "Gpx4", "Trp53", "Spp1", "Kap", "Gas5", "Brap")
star_dt <- dj[name %in% star, .(
  n_DOWN_all = sum(diffType == "Down"),
  n_DOWN_HC50 = sum(diffType == "Down" & a_pct >= 0.5),
  n_DOWN_HC70 = sum(diffType == "Down" & a_pct >= 0.7),
  n_DOWN_HC50_DRACH = sum(diffType == "Down" & a_pct >= 0.5 & is_drach),
  n_UP_all = sum(diffType == "Up"),
  n_UP_HC50 = sum(diffType == "Up" & b_pct >= 0.5)), by = name]
print(star_dt)
fwrite(star_dt, file.path(data_dir, "m6A_highconf_star_genes_NEWPATH.tsv"), sep = "\t")
rm(dj, diff_dt); gc(verbose = FALSE)

# ============ Part B: paired stoichiometry =====================================
cat("== Part B: paired stoichiometry ==\n")
run_test <- function(dt, label) {
  wt <- suppressWarnings(wilcox.test(dt$mean_frac_Day1, dt$mean_frac_Sham, paired = TRUE))
  data.table(set = label, n = nrow(dt),
             median_NC = median(dt$mean_frac_Sham),
             median_Day1 = median(dt$mean_frac_Day1),
             median_delta = median(dt$mean_frac_Day1 - dt$mean_frac_Sham),
             pct_down = 100 * mean(dt$mean_frac_Day1 < dt$mean_frac_Sham),
             pct_down_gt10pp = 100 * mean((dt$mean_frac_Sham - dt$mean_frac_Day1) > 0.10),
             pct_up_gt10pp = 100 * mean((dt$mean_frac_Day1 - dt$mean_frac_Sham) > 0.10),
             wilcox_p = wt$p.value)
}
paired <- site_ctx[n_det_NC == 3 & n_det_Day1 == 3]
cat(sprintf("  fully paired sites: %s\n", format(nrow(paired), big.mark = ",")))
res_b <- rbind(
  run_test(paired[mean_frac_Sham >= 0.7],             "HC70 all"),
  run_test(paired[mean_frac_Sham >= 0.7 & is_drach],  "HC70 DRACH+"),
  run_test(paired[mean_frac_Sham >= 0.7 & !is_drach], "HC70 DRACH-"),
  run_test(paired[mean_frac_Sham >= 0.5],             "HC50 all"),
  run_test(paired[mean_frac_Sham >= 0.5 & is_drach],  "HC50 DRACH+"),
  run_test(paired[mean_frac_Sham >= 0.5 & !is_drach], "HC50 DRACH-"))
print(res_b, digits = 3)
fwrite(res_b, file.path(data_dir, "m6A_paired_stoichiometry_NEWPATH.tsv"), sep = "\t")

hc_sham <- site_ctx[n_det_NC >= 2 & mean_frac_Sham >= 0.7]
drop_test <- hc_sham[, .(
  n = .N,
  n_absent_all_Day1 = sum(n_det_Day1 == 0),
  pct_absent = 100 * mean(n_det_Day1 == 0),
  mean_Day1_when_present = mean(mean_frac_Day1[n_det_Day1 > 0], na.rm = TRUE)
), by = is_drach]
cat("  presence-drop test (Sham HC70, >=2/3 NC):\n")
print(drop_test)
fwrite(drop_test, file.path(data_dir, "m6A_presence_drop_NEWPATH.tsv"), sep = "\t")

cat("\n== ALL DONE ==\n")
