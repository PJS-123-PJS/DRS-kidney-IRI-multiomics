# =============================================================================
# Um caller-independence check (reviewer concern: Nm callers cross-confuse,
# so the headline Um +41.9pp gain needs evidence it is not a caller artifact)
# 2026-09-27
#
# Test 1: Um-anchored gain sites (D1-HC50 n=246, D1-HC70 n=51) queried in the
#         Am, Cm and Gm caller outputs. If the gain is a caller artifact, the
#         same sites would drift in ALL Nm callers; if Um-specific, they stay
#         flat in the sibling callers.
# Test 2: m6A DRACH+ HC70 sites (n = 5,174, canonical) queried in the Um
#         caller. A globally biased Um caller would also report drift at
#         sites carrying a chemically different modification.
#
# Input: site_ctx_chunks_{um,am,cm,gm}/ctx_<sample>.tsv.gz,
#        site_ctx_nc_day1_NEWPATH.tsv.gz (m6A anchor set)
# Output: data/um_caller_independence.tsv
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
nc <- c("NCs1", "NCs2", "NCs3"); d1 <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc, d1)

load_caller <- function(mod) {
  ctx_dir <- file.path(data_dir, paste0("site_ctx_chunks_", mod))
  rbindlist(lapply(six, function(s) {
    dt <- fread(file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz")))
    dt[, .(trans, pos, fraction, sample = s)]
  }))
}
grp_stats <- function(ctx) {
  ctx[, .(n_det_NC = sum(sample %chin% nc),
          n_det_Day1 = sum(sample %chin% d1),
          mean_frac_Sham = mean(fraction[sample %chin% nc]),
          mean_frac_Day1 = mean(fraction[sample %chin% d1])),
      by = .(trans, pos)]
}

# ---- build Um anchor sets ---------------------------------------------------
cat("== build Um anchors ==\n")
um <- load_caller("um")
um_ctx <- grp_stats(um)
um_paired <- um_ctx[n_det_NC == 3 & n_det_Day1 == 3]
um_d1hc50 <- um_paired[mean_frac_Day1 >= 0.5, .(trans, pos)]
um_d1hc70 <- um_paired[mean_frac_Day1 >= 0.7, .(trans, pos)]
cat("Um D1-HC50:", nrow(um_d1hc50), " D1-HC70:", nrow(um_d1hc70), "\n")

# ---- m6A DRACH+ HC70 anchor set ---------------------------------------------
cat("== m6A anchor set ==\n")
sc <- fread(file.path(data_dir, "site_ctx_nc_day1_NEWPATH.tsv.gz"))
m6a_set <- sc[n_det_NC == 3 & n_det_Day1 == 3 & mean_frac_Sham >= 0.7 & is_drach == TRUE,
              .(trans, pos)]
cat("m6A DRACH+ SH-HC70 sites:", nrow(m6a_set), "\n")
rm(sc); gc(verbose = FALSE)

# ---- query function ---------------------------------------------------------
query_set_in_caller <- function(qset, mod) {
  ctx <- load_caller(mod)
  hit <- ctx[qset, on = .(trans, pos), nomatch = 0]
  empty <- data.table(caller = mod, n_sites_detected_any = 0L,
                      n_sites_both_groups = 0L, median_frac_Sham = NA_real_,
                      median_frac_Day1 = NA_real_, median_delta_pp = NA_real_,
                      wilcox_p = NA_real_)
  if (nrow(hit) < 6) return(empty)
  gs <- hit[, .(n_det = .N,
                mean_frac = mean(fraction)), by = .(trans, pos,
              grp = ifelse(sample %chin% nc, "Sham", "Day1"))]
  wide <- dcast(gs, trans + pos ~ grp, value.var = c("n_det", "mean_frac"))
  wide[is.na(wide)] <- NA
  both <- wide[!is.na(mean_frac_Sham) & !is.na(mean_frac_Day1)]
  data.table(
    caller = mod,
    n_sites_detected_any = uniqueN(hit$trans) ,
    n_sites_both_groups = nrow(both),
    median_frac_Sham = round(median(both$mean_frac_Sham), 4),
    median_frac_Day1 = round(median(both$mean_frac_Day1), 4),
    median_delta_pp = round(100 * median(both$mean_frac_Day1 - both$mean_frac_Sham), 2),
    wilcox_p = if (nrow(both) >= 10)
      signif(wilcox.test(both$mean_frac_Day1, both$mean_frac_Sham, paired = TRUE)$p.value, 3)
      else NA_real_)
}

out <- rbindlist(list(
  cbind(query_set = "Um D1-HC50 gain sites (n=246)", query_set_in_caller(um_d1hc50, "am")),
  cbind(query_set = "Um D1-HC50 gain sites (n=246)", query_set_in_caller(um_d1hc50, "cm")),
  cbind(query_set = "Um D1-HC50 gain sites (n=246)", query_set_in_caller(um_d1hc50, "gm")),
  cbind(query_set = "Um D1-HC70 gain sites (n=51)",  query_set_in_caller(um_d1hc70, "am")),
  cbind(query_set = "Um D1-HC70 gain sites (n=51)",  query_set_in_caller(um_d1hc70, "cm")),
  cbind(query_set = "Um D1-HC70 gain sites (n=51)",  query_set_in_caller(um_d1hc70, "gm")),
  cbind(query_set = "Um D1-HC50 gain sites (self-check, um)", query_set_in_caller(um_d1hc50, "um")),
  cbind(query_set = "m6A DRACH+ HC70 sites (n=5174) in um caller", query_set_in_caller(m6a_set, "um"))
))

# detection-rate context: what fraction of query sites each sibling caller
# even reports (baseline detectability)
det <- rbindlist(lapply(c("am", "cm", "gm"), function(md) {
  ctx <- load_caller(md)
  sites <- unique(ctx[, .(trans, pos)])
  data.table(caller = md,
             n_sites_total = nrow(sites),
             overlap_UmD1HC50 = nrow(fintersect(sites, um_d1hc50)))
}))
cat("\n== sibling-caller total site counts and overlap with Um anchor ==\n")
print(det)

fwrite(out, file.path(data_dir, "um_caller_independence.tsv"), sep = "\t")
cat("\n== RESULTS ==\n"); print(out)
cat("\nwritten: data/um_caller_independence.tsv\n")
