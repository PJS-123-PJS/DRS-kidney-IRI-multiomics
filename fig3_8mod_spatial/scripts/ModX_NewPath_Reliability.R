# =============================================================================
# Generic per-mod reliability battery (new canonical path) — v4-pattern
# 2026-07-22 — usage: Rscript ModX_NewPath_Reliability.R <mod>
#   mod in: Um PSU inosine Am Cm Gm   (m6a/m5c have dedicated scripts)
#
# None of these mods has a DRACH-class strict consensus, so the evidence
# hierarchy is: (1) site_diff direction symmetry, (2) data-driven 5-mer gate
# tiers (descriptive), (3) PAIRED STOICHIOMETRY (arbiter), (4) presence-drop,
# (5) HC counts. Produces HC gene sets for the co-modification panel.
#
# Products (data/, <mod> prefix):
#   <mod>_threshold_bins_NEWPATH.tsv  <mod>_highconf_counts_NEWPATH.tsv
#   <mod>_motif_enrichment_NEWPATH.tsv  <mod>_highconf_tier_summary_NEWPATH.tsv
#   <mod>_diff_site_context_NEWPATH.tsv.gz  <mod>_paired_stoichiometry_NEWPATH.tsv
#   <mod>_presence_drop_NEWPATH.tsv  geneset_<MOD>_{UP,DOWN}_HC50_NEWPATH.txt
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
mod <- args[1]
base <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
rds_cands <- file.path(base, mod, c(paste0("all_", mod, ".rds"),
                                    paste0("all_", tolower(mod), ".rds")))
RDS <- rds_cands[file.exists(rds_cands)][1]
stopifnot(!is.na(RDS))
diff_dir <- file.path(base, mod, "site_diff", "NC_vs_Day1")
DIFF <- list.files(diff_dir, pattern = "per_site\\.diff\\.add_name\\.xls$",
                   full.names = TRUE)[1]
stopifnot(!is.na(DIFF))
cat("mod:", mod, "\nRDS:", RDS, "\nDIFF:", DIFF, "\n")

data_dir <- "../data"
nc_samples <- c("NCs1", "NCs2", "NCs3")
d1_samples <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc_samples, d1_samples)
PFX <- tolower(mod)

cat("== loading rds ==\n")
t0 <- Sys.time()
r <- readRDS(RDS)
setDT(r)
setnames(r, tolower(names(r)))
cat(sprintf("  %s rows in %.1f min | frac==0: %.4f | depth>=20: %.4f\n",
            format(nrow(r), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            mean(r$fraction == 0), mean(r$depth >= 20)))
stopifnot(all(c("trans", "pos", "fraction", "motif", "sample") %in% names(r)))
r[, c("depth", "group", "chr", "gene", "rel_location") := NULL]
gc(verbose = FALSE)
all_samples <- unique(r$sample)

bin_out   <- file.path(data_dir, paste0(PFX, "_threshold_bins_NEWPATH.tsv"))
hc_out    <- file.path(data_dir, paste0(PFX, "_highconf_counts_NEWPATH.tsv"))
motif_out <- file.path(data_dir, paste0(PFX, "_motif_persample_NEWPATH.tsv"))
ctx_dir   <- file.path(data_dir, paste0("site_ctx_chunks_", PFX))
dir.create(ctx_dir, showWarnings = FALSE, recursive = TRUE)
first <- TRUE
for (s in all_samples) {
  sub <- r[sample == s, .(trans, pos, fraction, gene_struc, motif)]
  grp <- ifelse(grepl("^NC", s), "Sham",
         ifelse(grepl("^I", s), "IS", sub("s[0-9]+$", "", s)))
  b <- sub[, .(N = .N), by = .(bin = pmin(9L, as.integer(floor(fraction * 10))))]
  b[, `:=`(sample = s, group = grp)]
  fwrite(b, bin_out, sep = "\t", append = !first)
  hc <- data.table(sample = s, group = grp,
                   N_hc50 = sub[fraction >= 0.5, .N],
                   N_hc70 = sub[fraction >= 0.7, .N])
  fwrite(hc, hc_out, sep = "\t", append = !first)
  mc <- sub[, .(n_hc70 = sum(fraction >= 0.7),
                n_lo   = sum(fraction <  0.2),
                n_all  = .N), by = motif]
  mc[, sample := s]
  fwrite(mc, motif_out, sep = "\t", append = !first)
  if (s %chin% six) {
    fwrite(sub, file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")),
           sep = "\t", compress = "gzip")
  }
  cat(sprintf("  %s done (%s rows)\n", s, format(nrow(sub), big.mark = ",")))
  rm(sub, b, hc, mc); gc(verbose = FALSE)
  first <- FALSE
}
rm(r); gc(verbose = FALSE)

# data-driven gate (descriptive only for non-m6A mods)
mc <- fread(motif_out)
mp <- mc[, .(n_hc70 = sum(n_hc70), n_lo = sum(n_lo), n_all = sum(n_all)), by = motif]
tot_hc <- sum(mp$n_hc70); tot_lo <- sum(mp$n_lo)
mp[, enrich := (n_hc70 / tot_hc) / (n_lo / tot_lo)]
setorder(mp, -enrich)
gate <- mp[enrich >= 3 & n_hc70 >= 100]$motif
cat(sprintf("  gate: %d motifs (of %d seen)\n", length(gate), nrow(mp)))
fwrite(mp, file.path(data_dir, paste0(PFX, "_motif_enrichment_NEWPATH.tsv")), sep = "\t")

ctx_list <- lapply(six, function(s) {
  dt <- fread(file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")))
  dt[, sample := s]
  dt[, grp := ifelse(s %chin% nc_samples, "Sham", "Day1")]
  dt
})
ctx <- rbindlist(ctx_list); rm(ctx_list); gc(verbose = FALSE)
site_ctx <- ctx[, .(
  motif = motif[1],
  gene_struc = gene_struc[1],
  n_det_NC = sum(grp == "Sham"),
  n_det_Day1 = sum(grp == "Day1"),
  mean_frac_Sham = mean(fraction[grp == "Sham"]),
  mean_frac_Day1 = mean(fraction[grp == "Day1"])
), by = .(trans, pos)]
site_ctx[, is_gate := motif %chin% gate]
rm(ctx); gc(verbose = FALSE)
cat(sprintf("  unique sites: %s | gate+: %.1f%%\n",
            format(nrow(site_ctx), big.mark = ","), 100 * mean(site_ctx$is_gate)))

# Part A: tiers
diff_dt <- fread(DIFF)
stopifnot(ncol(diff_dt) == 9)
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
  n <- nrow(d); k <- sum(d$is_gate)
  ci <- if (n > 0) binom.test(k, n)$conf.int else c(NA_real_, NA_real_)
  data.table(tier = id, N = n, gate_pos = k, gate_pct = 100 * k / n,
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
fwrite(tier_dt, file.path(data_dir, paste0(PFX, "_highconf_tier_summary_NEWPATH.tsv")), sep = "\t")
fwrite(dj[, .(trans, pos, name, diffType, a_pct, b_pct, effect_size,
              motif, is_gate, gene_struc)],
       file.path(data_dir, paste0(PFX, "_diff_site_context_NEWPATH.tsv.gz")),
       sep = "\t", compress = "gzip")

g_up50 <- sort(unique(dj[diffType == "Up"   & b_pct >= 0.5 & is_gate]$name))
g_dn50 <- sort(unique(dj[diffType == "Down" & a_pct >= 0.5 & is_gate]$name))
g_up50 <- g_up50[g_up50 != "-"]; g_dn50 <- g_dn50[g_dn50 != "-"]
# fallback: if gate captures nothing meaningful, HC50 sets without gate
if (length(g_up50) < 20) {
  cat("  NOTE: gate+ UP set <20 genes — writing ungated HC50 sets instead\n")
  g_up50 <- sort(unique(dj[diffType == "Up"   & b_pct >= 0.5]$name))
  g_dn50 <- sort(unique(dj[diffType == "Down" & a_pct >= 0.5]$name))
  g_up50 <- g_up50[g_up50 != "-"]; g_dn50 <- g_dn50[g_dn50 != "-"]
}
writeLines(g_up50, file.path(data_dir, paste0("geneset_", toupper(mod), "_UP_HC50_NEWPATH.txt")))
writeLines(g_dn50, file.path(data_dir, paste0("geneset_", toupper(mod), "_DOWN_HC50_NEWPATH.txt")))
cat(sprintf("  UP_HC50 genes: %d | DOWN_HC50 genes: %d\n", length(g_up50), length(g_dn50)))
rm(dj, diff_dt); gc(verbose = FALSE)

# Part B: paired stoichiometry
run_test <- function(dt, label) {
  if (nrow(dt) < 10) return(data.table(set = label, n = nrow(dt)))
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
  run_test(paired[mean_frac_Sham >= 0.7],              "HC70 all"),
  run_test(paired[mean_frac_Sham >= 0.7 & is_gate],    "HC70 gate+"),
  run_test(paired[mean_frac_Sham >= 0.7 & !is_gate],   "HC70 gate-"),
  run_test(paired[mean_frac_Sham >= 0.5],              "HC50 all"),
  run_test(paired[mean_frac_Sham >= 0.5 & is_gate],    "HC50 gate+"),
  run_test(paired[mean_frac_Sham >= 0.5 & !is_gate],   "HC50 gate-"),
  fill = TRUE)
print(res_b, digits = 3)
fwrite(res_b, file.path(data_dir, paste0(PFX, "_paired_stoichiometry_NEWPATH.tsv")), sep = "\t")

hc_sham <- site_ctx[n_det_NC >= 2 & mean_frac_Sham >= 0.7]
drop_test <- hc_sham[, .(
  n = .N,
  n_absent_all_Day1 = sum(n_det_Day1 == 0),
  pct_absent = 100 * mean(n_det_Day1 == 0),
  mean_Day1_when_present = mean(mean_frac_Day1[n_det_Day1 > 0], na.rm = TRUE),
  median_delta_when_present = median(mean_frac_Day1[n_det_Day1 > 0] -
                                     mean_frac_Sham[n_det_Day1 > 0], na.rm = TRUE)
), by = is_gate]
cat("  presence-drop (Sham HC70, >=2/3 NC):\n")
print(drop_test)
fwrite(drop_test, file.path(data_dir, paste0(PFX, "_presence_drop_NEWPATH.tsv")), sep = "\t")

cat("\n== ", mod, "BATTERY DONE ==\n")
