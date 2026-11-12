# =============================================================================
# Day1 poly(A)-shortened genes: where do they go at Day3 / Day7?
# 2026-09-27 — functional trajectory readout for the structural layer.
# If the Day1 shortening were passive degradation, shortened genes should drift
# back to baseline; if it is a directed program, the SAME genes should show the
# Day3 overshoot (the global trajectory dips at Day1 then overshoots at Day3).
#
# Sets: Day1-shortened / Day1-lengthened (NC_vs_Day1, FDR<0.05 & |eff|>0.3)
# Trajectory: each gene's Diff_median vs Sham at Day1, Day3, Day7 (all vs NC).
# Control: genes unchanged at Day1 (FDR>=0.05 or |eff|<=0.3, both-directions).
# Output: data/polyA_trajectory_after_Day1.tsv + printed summary
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

base <- "/mnt/f/BN2505012AH02S66N1/DirectRNA/result/06.polyA/length/02.diff"
out_d <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig4_as_landscape/data"

read_diff <- function(cmp) {
  d <- fread(file.path(base, cmp, paste0(cmp, ".polya_diff.all.xls")))
  d[, .(Transcript, Name, Diff_median, Effect_size, FDR,
        ctrl_med = Control_median, treat_med = Treatment_median)]
}
d1 <- read_diff("NC_vs_Day1")
d3 <- read_diff("NC_vs_Day3")
d7 <- read_diff("NC_vs_Day7")
setnames(d1, c("Diff_median", "Effect_size", "FDR", "ctrl_med", "treat_med"),
         paste0(c("Diff_median", "Effect_size", "FDR", "ctrl_med", "treat_med"), "_D1"))
setnames(d3, c("Diff_median", "Effect_size", "FDR", "ctrl_med", "treat_med"),
         paste0(c("Diff_median", "Effect_size", "FDR", "ctrl_med", "treat_med"), "_D3"))
setnames(d7, c("Diff_median", "Effect_size", "FDR", "ctrl_med", "treat_med"),
         paste0(c("Diff_median", "Effect_size", "FDR", "ctrl_med", "treat_med"), "_D7"))

m <- merge(merge(d1, d3, by = c("Transcript", "Name"), all = FALSE),
           d7, by = c("Transcript", "Name"), all = FALSE)
cat("transcripts with D1/D3/D7 all measured:", nrow(m), "\n")

sig_short <- m[FDR_D1 < 0.05 & Effect_size_D1 < -0.3]
sig_long  <- m[FDR_D1 < 0.05 & Effect_size_D1 > 0.3]
unchanged <- m[FDR_D1 >= 0.05 | abs(Effect_size_D1) <= 0.3]
cat("Day1-shortened:", nrow(sig_short), " Day1-lengthened:", nrow(sig_long),
    " unchanged:", nrow(unchanged), "\n")

summ <- function(dt, tag) {
  cat(sprintf(
    "[%s] n=%d\n  vs Sham median-diff (nt): Day1 %.1f -> Day3 %+.1f -> Day7 %+.1f\n  fraction above Sham baseline: Day1 %.1f%% -> Day3 %.1f%% -> Day7 %.1f%%\n",
    tag, nrow(dt),
    median(dt$Diff_median_D1), median(dt$Diff_median_D3), median(dt$Diff_median_D7),
    100 * mean(dt$Diff_median_D1 > 0), 100 * mean(dt$Diff_median_D3 > 0),
    100 * mean(dt$Diff_median_D7 > 0)))
  dt
}
s <- summ(sig_short, "Day1-shortened genes")
l <- summ(sig_long, "Day1-lengthened genes")
u <- summ(unchanged, "Day1-unchanged (control)")

# overshoot test: Day1-shortened genes at Day3 vs unchanged genes at Day3
wt <- wilcox.test(s$Diff_median_D3, u$Diff_median_D3)
cat(sprintf("\nDay3 distribution shift (shortened vs unchanged): W test p = %.2e\n",
            wt$p.value))
cat(sprintf("Day1-shortened genes above baseline at Day3: %.1f%% (control %.1f%%)\n",
            100 * mean(s$Diff_median_D3 > 0), 100 * mean(u$Diff_median_D3 > 0)))
tab <- matrix(c(sum(s$Diff_median_D3 > 0), sum(s$Diff_median_D3 <= 0),
                sum(u$Diff_median_D3 > 0), sum(u$Diff_median_D3 <= 0)), 2, byrow = TRUE)
print(fisher.test(tab))

# same-direction persistence check D1->D7
cat(sprintf("\nDay1-shortened still shortened at Day7: %.1f%%\n",
            100 * mean(s$Diff_median_D7 < 0)))
cat(sprintf("Day1-lengthened still lengthened at Day7: %.1f%%\n",
            100 * mean(l$Diff_median_D7 > 0)))

out <- rbind(
  cbind(set = "Day1_shortened", s[, .(Transcript, Name, Diff_median_D1,
      Diff_median_D3, Diff_median_D7)]),
  cbind(set = "Day1_lengthened", l[, .(Transcript, Name, Diff_median_D1,
      Diff_median_D3, Diff_median_D7)]))
fwrite(out, file.path(out_d, "polyA_trajectory_after_Day1.tsv"), sep = "\t")
cat("\nwritten: fig4_as_landscape/data/polyA_trajectory_after_Day1.tsv\n")
