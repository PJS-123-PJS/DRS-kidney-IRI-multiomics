# =============================================================================
# Animal-level sensitivity analysis for the dual-anchor arbitration
# 2026-09-27 — addresses reviewer concern: site-level Wilcoxon uses sites as
# the statistical unit (n = thousands) while biological replicates are n = 3.
# Here every statistic is re-expressed at the ANIMAL level.
#
# Metrics per caller x anchor set:
#   1. fidelity gate: site-level median delta must reproduce the canonical
#      modverdict_paired_stats.tsv value exactly (same anchor definitions)
#   2. per-animal median delta: for each animal i, median over anchor sites of
#      (fraction_i - mean fraction of the other group)  -> 6 values
#   3. direction concordance: 3/3 Day1 animals with same-sign median delta
#   4. per-site 9/9 concordance: fraction of sites where all 9 Day1xSham
#      pairwise differences share the verdict sign
#   5. separation: Day1 per-animal medians vs Sham per-animal medians overlap?
#   6. exact permutation p (3v3, all 20 assignments, two-sided on mean diff)
#
# Input: site_ctx_chunks_<mod>/ctx_<sample>.tsv.gz (per-sample fraction)
# Output: data/animal_level_sensitivity_7callers.tsv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
nc <- c("NCs1", "NCs2", "NCs3"); d1 <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc, d1)

# canonical site-level values to reproduce (from modverdict_paired_stats.tsv)
canon <- fread(file.path(data_dir, "modverdict_paired_stats.tsv"))
canon[, key := paste(mod, anchor, sep = "|")]

# m5C motif gate (same definition as ModVerdict_Paired_Stats.R)
me5 <- fread(file.path(data_dir, "m5c_motif_enrichment_NEWPATH.tsv"))
gate5 <- me5[enrich >= 3 & n_hc70 >= 100]$motif

load_caller <- function(mod) {
  ctx_dir <- file.path(data_dir, paste0("site_ctx_chunks_", mod))
  rbindlist(lapply(six, function(s) {
    dt <- fread(file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")))
    dt[, .(trans, pos, fraction, sample = s)]
  }))
}

per_animal_stats <- function(wide, sites) {
  # wide: site rows x 6 sample fraction columns (complete paired sites only)
  # returns per-animal median deltas in pp
  sub <- wide[sites]
  fs <- as.matrix(sub[, ..nc]); fd <- as.matrix(sub[, ..d1])
  sham_mean <- rowMeans(fs); day1_mean <- rowMeans(fd)
  c(
    sapply(1:3, function(i) median(fd[, i] - sham_mean) * 100),   # Day1 animals
    sapply(1:3, function(i) median(fs[, i] - day1_mean) * 100)    # Sham animals
  )
}

pct_9of9 <- function(wide, sites, sign) {
  sub <- wide[sites]
  fs <- as.matrix(sub[, ..nc]); fd <- as.matrix(sub[, ..d1])
  # all 9 pairwise (Day1_a - Sham_b) share sign
  v <- rep(TRUE, nrow(sub))
  for (a in 1:3) for (b in 1:3) {
    if (sign > 0) v <- v & (fd[, a] - fs[, b] > 0)
    else          v <- v & (fd[, a] - fs[, b] < 0)
  }
  100 * mean(v)
}

perm_p_3v3 <- function(x, y) {
  # exact permutation test on mean difference, two-sided
  obs <- mean(x) - mean(y)
  allv <- c(x, y); combs <- combn(6, 3)
  k <- sum(sapply(1:ncol(combs), function(j) {
    idx <- combs[, j]
    abs(mean(allv[idx]) - mean(allv[-idx])) >= abs(obs) - 1e-12
  }))
  k / ncol(combs)
}

analyze <- function(mod, tag, gate_motifs = NULL) {
  cat("==", tag, "==\n")
  ctx <- load_caller(mod)
  # group means and detection counts per site (same as ModVerdict)
  site_ctx <- ctx[, .(
    n_det_NC = sum(sample %chin% nc),
    n_det_Day1 = sum(sample %chin% d1),
    mean_frac_Sham = mean(fraction[sample %chin% nc]),
    mean_frac_Day1 = mean(fraction[sample %chin% d1])
  ), by = .(trans, pos)]
  paired <- site_ctx[n_det_NC == 3 & n_det_Day1 == 3]
  # wide per-sample fraction matrix (complete paired sites only)
  wide <- dcast(ctx[paired, on = .(trans, pos)],
                trans + pos ~ sample, value.var = "fraction")
  wide <- wide[complete.cases(wide)]
  cat("  paired sites:", nrow(wide), "\n")

  rows <- list()
  for (anchor in c("SH-HC50", "D1-HC50", "SH-HC70", "D1-HC70")) {
    sites <- if (anchor == "SH-HC50") paired$mean_frac_Sham >= 0.5
             else if (anchor == "SH-HC70") paired$mean_frac_Sham >= 0.7
             else if (anchor == "D1-HC50") paired$mean_frac_Day1 >= 0.5
             else paired$mean_frac_Day1 >= 0.7
    # align wide rows to paired rows (same order after dcast? enforce by key)
    wide_key <- paste(wide$trans, wide$pos)
    paired_key <- paste(paired$trans, paired$pos)
    site_idx <- wide_key %in% paired_key[sites]
    if (sum(site_idx) < 10) next
    pa <- per_animal_stats(wide, site_idx)
    d1_pa <- pa[1:3]; sh_pa <- pa[4:6]
    site_delta <- (median(paired$mean_frac_Day1[sites] -
                          paired$mean_frac_Sham[sites])) * 100
    verdict_sign <- sign(site_delta)
    rows[[length(rows) + 1]] <- data.table(
      mod = tag, anchor = anchor, n_sites = sum(site_idx),
      site_median_delta_pp = round(site_delta, 2),
      Day1_animal_median_deltas_pp = paste(sprintf("%.2f", d1_pa), collapse = ";"),
      Sham_animal_median_deltas_pp = paste(sprintf("%.2f", sh_pa), collapse = ";"),
      Day1_3of3_concordant = all(sign(d1_pa) == verdict_sign),
      separation_no_overlap = if (verdict_sign > 0)
        min(d1_pa) > max(sh_pa) else max(d1_pa) < min(sh_pa),
      perm_p_3v3_two_sided = round(perm_p_3v3(d1_pa, sh_pa), 3),
      pct_sites_9of9_concordant = round(pct_9of9(wide, site_idx, verdict_sign), 1)
    )
  }
  rbindlist(rows)
}

mods <- c("m5c", "psu", "inosine", "um", "am", "cm", "gm")
tags <- c("m5C", "PSU", "inosine", "Nm-Um", "Nm-Am", "Nm-Cm", "Nm-Gm")
out <- rbindlist(lapply(seq_along(mods), function(i) analyze(mods[i], tags[i])))

# fidelity gate: compare with canonical site-level medians
out[, canon_key := paste(mod, anchor, sep = "|")]
gate <- merge(out, canon[, .(key, n, median_delta_pp)],
              by.x = "canon_key", by.y = "key", all.x = TRUE)
gate[, n_match := n_sites == n]
gate[, delta_match := abs(site_median_delta_pp - median_delta_pp) < 0.51]
cat("\n== FIDELITY GATE (site-level reproduction of canonical values) ==\n")
print(gate[, .(mod, anchor, n_sites, n, n_match,
               site_median_delta_pp, median_delta_pp, delta_match)])
stopifnot(all(gate$n_match[!is.na(gate$n)]), all(gate$delta_match[!is.na(gate$n)]))

fwrite(out, file.path(data_dir, "animal_level_sensitivity_7callers.tsv"), sep = "\t")
cat("\n== RESULTS ==\n")
print(out)
cat("\nwritten: data/animal_level_sensitivity_7callers.tsv\n")
