# =============================================================================
# Transcript-level verdict for the m6A-down "wave" (new canonical path)
# 2026-07-22
#
# The paper's m6A-suppressive-wave story came from trans_diff: ~746-813
# transcripts "DOWN" at Day1. Site-level battery showed DOWN sites are
# non-DRACH noise. Here we quantify the transcript level directly:
#   (1) stoichiometry distribution of DOWN transcripts (Group1_pct = Sham)
#   (2) expression-dilution cross-tab (all_diff from trans_diff_Expression)
#   (3) how many DOWN transcripts contain ANY high-confidence (>=0.5) DRACH
#       site in the per-site context
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

TD   <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/m6a/trans_diff/NC_vs_Day1/NC_vs_Day1.m6a.per_trans.diff.add_name.xls"
CORR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/m6a/trans_diff_Expression/NC_vs_Day1/NC_vs_Day1.m6a_expression.corr.new.xls.gz"
data_dir <- "../data"

td <- fread(TD)
setnames(td, tolower(names(td)))
# counts arrive as "a:82" strings -> numeric
td[, g1_counts := as.integer(sub(".*:", "", group1_counts))]
td[, g2_counts := as.integer(sub(".*:", "", group2_counts))]
cat(sprintf("trans_diff rows: %s | Up=%s Down=%s\n",
            format(nrow(td), big.mark = ","),
            sum(td$difftype == "Up"), sum(td$difftype == "Down")))

dn <- td[difftype == "Down"]
up <- td[difftype == "Up"]
summ <- data.table(
  set = c("DOWN", "UP"),
  n = c(nrow(dn), nrow(up)),
  median_pct_G1 = c(median(dn$group1_pct_modified), median(up$group1_pct_modified)),
  median_pct_G2 = c(median(dn$group2_pct_modified), median(up$group2_pct_modified)),
  pct_ge10_G1 = c(100 * mean(dn$group1_pct_modified >= 0.10),
                  100 * mean(up$group1_pct_modified >= 0.10)),
  pct_ge5_G1  = c(100 * mean(dn$group1_pct_modified >= 0.05),
                  100 * mean(up$group1_pct_modified >= 0.05)),
  pct_ge50_G1 = c(100 * mean(dn$group1_pct_modified >= 0.50),
                  100 * mean(up$group1_pct_modified >= 0.50)))
print(summ, digits = 3)
fwrite(summ, file.path(data_dir, "m6A_transdiff_stoich_summary_NEWPATH.tsv"), sep = "\t")

# ---- expression-dilution cross-tab --------------------------------------------
corr <- fread(CORR)
setnames(corr, tolower(names(corr)))
corr[, g1_counts := as.integer(sub(".*:", "", group1_counts))]
corr[, g2_counts := as.integer(sub(".*:", "", group2_counts))]
cat("corr cols:", paste(names(corr), collapse = ", "), "\n")
md <- corr[methyl_diff == "Down"]
cat(sprintf("\nDOWN transcripts with expression data: %d\n", nrow(md)))
xt <- md[, .N, by = .(all_diff)]
xt[, pct := round(100 * N / sum(N), 1)]
print(xt[order(-N)])

# stoichiometry tier x dilution status
md[, stoi_tier := cut(group1_pct_modified,
                      breaks = c(-Inf, 0.02, 0.05, 0.10, 0.50, Inf),
                      labels = c("<2%", "2-5%", "5-10%", "10-50%", ">=50%"))]
xt2 <- md[, .N, by = .(stoi_tier, all_diff)]
print(dcast(xt2, stoi_tier ~ all_diff, value.var = "N", fill = 0))
fwrite(xt2, file.path(data_dir, "m6A_transdiff_stoich_x_dilution_NEWPATH.tsv"), sep = "\t")

# dilution math check: for Down_Up, expected fraction change from totals alone
md[, dilution_fc := (group2_total / group1_total)]
cat("\nDown_Up transcripts: median expression fold = %.2f (fraction mechanically diluted by same factor)\n")
print(md[all_diff == "Down_Up",
         .(n = .N, median_total_fc = median(dilution_fc),
           median_pct_fc = median(group2_pct_modified / group1_pct_modified))])
cat("\nDown_Down: count-drop follows expression?\n")
print(md[all_diff == "Down_Down",
         .(n = .N, median_total_fc = median(dilution_fc),
           median_modified_fc = median(g2_counts / pmax(g1_counts, 1),
                                       na.rm = TRUE))])
cat("== DONE ==\n")
