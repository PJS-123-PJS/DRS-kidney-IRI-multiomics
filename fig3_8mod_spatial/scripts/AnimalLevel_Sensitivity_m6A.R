# =============================================================================
# Animal-level sensitivity analysis — m6A (the headline verdict)
# 2026-09-27 — same metrics as AnimalLevel_Sensitivity_7callers.R, plus the
# DRACH+ HC70 subset (the canonical n = 5,174 set quoted in the main text).
# Input: all_m6a.rds (per-sample per-site fraction, company pipeline) +
#        site_ctx_nc_day1_NEWPATH.tsv.gz (paired-site definitions, is_drach)
# Output: data/animal_level_sensitivity_m6a.tsv
# Fidelity gate: site-level medians must reproduce modverdict_paired_stats.tsv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
nc <- c("NCs1", "NCs2", "NCs3"); d1 <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc, d1)

cat("== load paired-site definitions ==\n")
sc <- fread(file.path(data_dir, "site_ctx_nc_day1_NEWPATH.tsv.gz"))
paired <- sc[n_det_NC == 3 & n_det_Day1 == 3,
             .(trans, pos, is_drach, mean_frac_Sham, mean_frac_Day1)]
cat("paired sites:", nrow(paired), "\n")
rm(sc); gc(verbose = FALSE)

cat("== load all_m6a.rds (797MB, 6 samples only) ==\n")
all <- readRDS("/mnt/f/BN2505012AH02S66N1/DirectRNA/05.methy/m6a/all_m6a.rds")
setDT(all)
sub <- all[sample %chin% six, .(Trans, Pos, Fraction, sample)]
rm(all); gc(verbose = FALSE)
cat("rows after sample filter:", nrow(sub), "\n")

m <- sub[paired, on = .(Trans = trans, Pos = pos), nomatch = 0]
rm(sub); gc(verbose = FALSE)
cat("merged rows:", nrow(m), "\n")
dup <- m[, .N, by = .(Trans, Pos, sample)][N > 1]
cat("duplicate site-sample rows:", nrow(dup), "\n")

wide <- dcast(m, Trans + Pos ~ sample, value.var = "Fraction", fun.aggregate = mean)
wide <- wide[complete.cases(wide)]
wide_key <- paste(wide$Trans, wide$Pos)
paired_key <- paste(paired$trans, paired$pos)
idx_in_paired <- wide_key %in% paired_key
wide <- wide[idx_in_paired]
pm <- paired[match(wide_key[idx_in_paired], paired_key)]
stopifnot(nrow(wide) == nrow(pm))
cat("wide matrix:", nrow(wide), "sites x 6 animals\n")

per_animal_stats <- function(sel) {
  fs <- as.matrix(wide[sel, ..nc]); fd <- as.matrix(wide[sel, ..d1])
  sham_mean <- rowMeans(fs); day1_mean <- rowMeans(fd)
  c(sapply(1:3, function(i) median(fd[, i] - sham_mean) * 100),
    sapply(1:3, function(i) median(fs[, i] - day1_mean) * 100))
}
pct_9of9 <- function(sel, sg) {
  fs <- as.matrix(wide[sel, ..nc]); fd <- as.matrix(wide[sel, ..d1])
  v <- rep(TRUE, nrow(fs))
  for (a in 1:3) for (b in 1:3)
    v <- if (sg > 0) v & (fd[, a] - fs[, b] > 0) else v & (fd[, a] - fs[, b] < 0)
  100 * mean(v)
}
perm_p_3v3 <- function(x, y) {
  obs <- mean(x) - mean(y); allv <- c(x, y); combs <- combn(6, 3)
  sum(sapply(1:20, function(j) {
    idx <- combs[, j]
    abs(mean(allv[idx]) - mean(allv[-idx])) >= abs(obs) - 1e-12
  })) / 20
}

rows <- list()
add_row <- function(tag, anchor, sel) {
  if (sum(sel) < 10) return(NULL)
  pa <- per_animal_stats(sel); d1_pa <- pa[1:3]; sh_pa <- pa[4:6]
  sd_ <- median(pm$mean_frac_Day1[sel] - pm$mean_frac_Sham[sel]) * 100
  sg <- sign(sd_)
  rows[[length(rows) + 1]] <<- data.table(
    mod = tag, anchor = anchor, n_sites = sum(sel),
    site_median_delta_pp = round(sd_, 2),
    Day1_animal_median_deltas_pp = paste(sprintf("%.2f", d1_pa), collapse = ";"),
    Sham_animal_median_deltas_pp = paste(sprintf("%.2f", sh_pa), collapse = ";"),
    Day1_3of3_concordant = all(sign(d1_pa) == sg),
    separation_no_overlap = if (sg > 0) min(d1_pa) > max(sh_pa) else max(d1_pa) < min(sh_pa),
    perm_p_3v3_two_sided = round(perm_p_3v3(d1_pa, sh_pa), 3),
    pct_sites_9of9_concordant = round(pct_9of9(sel, sg), 1))
}

add_row("m6A", "SH-HC50", pm$mean_frac_Sham >= 0.5)
add_row("m6A", "D1-HC50", pm$mean_frac_Day1 >= 0.5)
add_row("m6A", "SH-HC70", pm$mean_frac_Sham >= 0.7)
add_row("m6A", "D1-HC70", pm$mean_frac_Day1 >= 0.7)
add_row("m6A DRACH+", "SH-HC70", pm$mean_frac_Sham >= 0.7 & pm$is_drach)
add_row("m6A DRACH+", "D1-HC70", pm$mean_frac_Day1 >= 0.7 & pm$is_drach)

out <- rbindlist(rows)
fwrite(out, file.path(data_dir, "animal_level_sensitivity_m6a.tsv"), sep = "\t")
cat("\n== RESULTS m6A ==\n"); print(out)

# fidelity gate vs canonical
canon <- fread(file.path(data_dir, "modverdict_paired_stats.tsv"))
canon[, key := paste(mod, anchor, sep = "|")]
g <- merge(out, canon[, .(key, n, median_delta_pp)],
#            by.x = "paste(mod,anchor,sep='|')", by.y = "key")
# # (data.table NSE workaround)
cat("\n== FIDELITY GATE ==\n")
for (i in seq_len(nrow(out))) {
  k <- paste(out$mod[i], out$anchor[i], sep = "|")
  cz <- canon[key == k]
  if (!nrow(cz)) { cat(k, ": no canonical row\n"); next }
  ok_n <- out$n_sites[i] == cz$n
  ok_d <- abs(out$site_median_delta_pp[i] - cz$median_delta_pp) < 0.51
  cat(sprintf("%-22s n %d vs %d [%s]  delta %.2f vs %.2f [%s]\n", k,
              out$n_sites[i], cz$n, ifelse(ok_n, "PASS", "FAIL"),
              out$site_median_delta_pp[i], cz$median_delta_pp,
              ifelse(ok_d, "PASS", "FAIL")))
  stopifnot(ok_n, ok_d)
}
cat("\nALL GATES PASS\n")
