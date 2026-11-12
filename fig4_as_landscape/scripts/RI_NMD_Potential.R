# =============================================================================
# RI -> PTC/NMD potential analysis (functional readout for the AS layer)
# 2026-09-27 — for significant Day1 vs Sham retained-intron events:
#   is the retained intron inside the CDS, and would its retention introduce
#   a premature termination codon (frameshift by length % 3 != 0, or an
#   in-frame stop codon)? Non-significant RI events serve as the control set.
# Simplified NMD rule (documented): CDS-internal RI with frameshift or
# in-frame stop = high NMD potential (the 50-nt EJC rule is not evaluated).
# Input: RI_isoform.dpsi.add_info.xls (NC_vs_Day1), Reference GTF + genome fa
# Output: data/RI_NMD_potential.tsv + console summary
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

base <- "/mnt/f/BN2505012AH02S66N1/DirectRNA/result"
ri_file <- file.path(base, "07.Structure_Analysis/03.AS/02.Diff_AS/NC_vs_Day1-diff/RI_isoform.dpsi.add_info.xls")
gtf_file <- file.path(base, "Reference/Mus_musculus.gtf")
fa_file  <- file.path(base, "Reference/Mus_musculus.fa")
tmp <- "/tmp/ri_nmd"; dir.create(tmp, showWarnings = FALSE)

ri <- fread(ri_file)
ri[, dpsi_abs := abs(dpsi)]
ri[, sig := dpsi_abs > 0.05 & `diffAS-pvalue` < 0.05]
ev <- unique(ri[sig == TRUE], by = "event_id")
ctl <- unique(ri[sig == FALSE], by = "event_id")
cat("significant RI events:", nrow(ev), " control RI events:", nrow(ctl), "\n")
cat("sig direction: dpsi>0 (more retention Day1):", sum(ev$dpsi > 0),
    " dpsi<0:", sum(ev$dpsi < 0), "\n")

# parse event coordinates: gene;RI:chr:upEnd:intronStart-intronEnd:downStart:strand
parse_ev <- function(eid) {
  core <- sub("^[^;]*;RI:", "", eid)
  strand <- substr(core, nchar(core), nchar(core))
  core2 <- substr(core, 1, nchar(core) - 2)   # drop ":strand" (strand "-"/"+" )
  p <- strsplit(core2, "[:-]")[[1]]
  list(chr = p[1], i_start = as.integer(p[3]), i_end = as.integer(p[4]),
       strand = strand)
}
parse_all <- function(dt) {
  m <- t(sapply(dt$event_id, function(x) unlist(parse_ev(x))))
  dt[, `:=`(chr = m[, 1], i_start = as.integer(m[, 2]),
            i_end = as.integer(m[, 3]), strand = m[, 4])]
  dt[, intron_len := i_end - i_start + 1]
  dt
}
ev <- parse_all(ev); ctl <- parse_all(ctl)
set.seed(7); if (nrow(ctl) > 400) ctl <- ctl[sample(.N, 400)]

# ---- subset GTF to relevant genes -------------------------------------------
genes <- unique(c(ev$gene_id, ctl$gene_id))
writeLines(genes, file.path(tmp, "genes.txt"))
cmd <- sprintf("grep -F -f %s %s > %s", file.path(tmp, "genes.txt"),
               gtf_file, file.path(tmp, "sub.gtf"))
cat("subsetting GTF ...\n"); system(cmd, intern = TRUE)
gtf <- fread(file.path(tmp, "sub.gtf"), header = FALSE, sep = "\t",
             col.names = c("chr", "src", "feature", "start", "end", "score",
                           "strand", "phase", "attr"))
cds <- gtf[feature == "CDS"]
cds[, gene := sub(".*gene_id \"?([^\";]+)\"?.*", "\\1", attr)]
cds[, trans := sub(".*transcript_id \"?([^\";]+)\"?.*", "\\1", attr)]
cat("CDS rows:", nrow(cds), "\n")

# ---- NMD potential per event -------------------------------------------------
# CDS-internal: intron overlaps any CDS exon span of its gene
# frameshift: intron_len %% 3 != 0 (retention shifts downstream frame)
# in-frame stop: scan retained intron sequence in the upstream reading frame
cds_by_gene <- split(cds, cds$gene)

classify <- function(chr_, s_, e_, strand_, gene_) {
  g <- cds_by_gene[[gene_]]
  if (is.null(g)) return(list(in_cds = NA, phase0 = NA_real_))
  ov <- g[chr == chr_ & start <= e_ & end >= s_]
  if (!nrow(ov)) return(list(in_cds = FALSE, phase0 = NA_real_))
  # upstream CDS length before intron start (for frame), same strand gene model
  up <- g[chr == chr_ & end < s_]
  up_len <- sum(up$end - up$start + 1)
  list(in_cds = TRUE, phase0 = as.double(up_len %% 3))
}

stops <- c("TAA", "TAG", "TGA")
scan_stops <- function(seq, phase0) {
  # seq = intron sequence in coding orientation (5'->3' of mRNA)
  n <- nchar(seq); if (n < 3) return(FALSE)
  for (i in seq(1 + phase0, n - 2, by = 3))
    if (substr(seq, i, i + 2) %in% stops) return(TRUE)
  FALSE
}
revcomp <- function(s) {
  s <- chartr("ACGTacgt", "TGCAtgca", s)
  paste(rev(strsplit(s, "")[[1]]), collapse = "")
}

fetch_seqs <- function(dt, tag) {
  regs <- unique(dt[, .(chr, i_start, i_end)])
  rf <- file.path(tmp, paste0("regions_", tag, ".txt"))
  fwrite(regs[, .(paste0(chr, ":", i_start, "-", i_end))], rf, col.names = FALSE)
  out <- file.path(tmp, paste0("seqs_", tag, ".fa"))
  system(sprintf("samtools faidx %s -r %s -o %s", fa_file, rf, out))
  fa <- readLines(out)
  ids <- grep("^>", fa)
  seqs <- list()
  for (i in seq_along(ids)) {
    nm <- sub("^>", "", fa[ids[i]])
    e <- if (i < length(ids)) ids[i + 1] - 1 else length(fa)
    seqs[[nm]] <- paste(fa[(ids[i] + 1):e], collapse = "")
  }
  seqs
}

run_set <- function(dt, tag) {
  cat("==", tag, "==\n")
  seqs <- fetch_seqs(dt, tag)
  res <- dt[, {
    cl <- classify(chr, i_start, i_end, strand, gene_id)
    key <- paste0(chr, ":", i_start, "-", i_end)
    sq <- seqs[[key]]
    if (!is.null(sq) && strand == "-") sq <- revcomp(sq)
    fs <- (intron_len %% 3) != 0
    stop_hit <- if (!is.null(sq) && isTRUE(cl$in_cds)) scan_stops(sq, cl$phase0) else NA
    .(in_cds = as.logical(cl$in_cds), phase0 = as.double(cl$phase0),
      frameshift = as.logical(fs), in_frame_stop = as.logical(stop_hit))
  }, by = .(event_id)]
  merge(dt[, .(event_id, dpsi, intron_len, gene_id)], res, by = "event_id")
}

r_ev <- run_set(ev, "sig")
r_ct <- run_set(ctl, "ctl")

summ <- function(r, tag) {
  ok <- r[!is.na(in_cds)]
  cat(sprintf("[%s] n=%d | CDS-internal: %d (%.1f%%) | of CDS-internal: frameshift %.1f%%, in-frame-stop %.1f%%, either(NMD-potential) %.1f%%\n",
      tag, nrow(r), sum(ok$in_cds), 100 * mean(ok$in_cds),
      100 * mean(ok$frameshift[ok$in_cds], na.rm = TRUE),
      100 * mean(ok$in_frame_stop[ok$in_cds], na.rm = TRUE),
      100 * mean(ok$in_cds & (ok$frameshift | ok$in_frame_stop), na.rm = TRUE) /
        mean(ok$in_cds) * 0 + 100 * mean((ok$frameshift | ok$in_frame_stop)[ok$in_cds], na.rm = TRUE)))
  ok
}
s1 <- summ(r_ev, "Day1-significant RI")
s2 <- summ(r_ct, "control RI")
# fisher: NMD-potential among CDS-internal
tab <- matrix(c(sum(s1$in_cds & (s1$frameshift | s1$in_frame_stop), na.rm = TRUE),
                sum(s1$in_cds & !(s1$frameshift | s1$in_frame_stop), na.rm = TRUE),
                sum(s2$in_cds & (s2$frameshift | s2$in_frame_stop), na.rm = TRUE),
                sum(s2$in_cds & !(s2$frameshift | s2$in_frame_stop), na.rm = TRUE)),
              nrow = 2, byrow = TRUE)
print(tab); print(fisher.test(tab))

fwrite(r_ev, "/mnt/d/desktop/DRS/DRS_TechShowcase/fig4_as_landscape/data/RI_NMD_potential_sig.tsv", sep = "\t")
fwrite(r_ct, "/mnt/d/desktop/DRS/DRS_TechShowcase/fig4_as_landscape/data/RI_NMD_potential_ctl.tsv", sep = "\t")
cat("written: fig4_as_landscape/data/RI_NMD_potential_{sig,ctl}.tsv\n")
