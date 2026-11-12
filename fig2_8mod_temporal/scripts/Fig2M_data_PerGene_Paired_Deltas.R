# =============================================================================
# Fig2F F3 (validated): per-gene paired site-level modification deltas
# Replaces the retired "Bidirectional modification wave" panel (Park7/Fis1/
# Malat1 had ZERO validated-set support; star-gene DOWN-HC sites were 0/9 DRACH+).
# Evidence: caller-independent paired sites (3/3 Sham + 3/3 Day1) from
# site_ctx_chunks_<mod>/, median per-gene frac delta (Day1 - Sham, pp).
# Genes: validated-set stars (Tmem176b/Ggt1/Dag1/Sec31a/Tgfbr2).
# Output: data/Fig2F3_per_gene_paired_deltas.tsv  (consumed by Fig2F_star_genes.R)
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages(library(data.table))

data_dir <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
out_f    <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig2_8mod_temporal/data/Fig2F3_per_gene_paired_deltas.tsv"

TARGETS <- c("Tmem176b", "Ggt1", "Dag1", "Sec31a", "Tgfbr2")

# ---- ENSMUSG -> gene name map (union of all per_trans add_name files) --------
map_list <- list()
for (m in c("m6a", "m5c", "psu", "inosine", "um", "am", "cm", "gm")) {
  name_f <- sprintf("/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy/%s/trans_diff/NC_vs_Day1/NC_vs_Day1.%s.per_trans.add_name.xls.gz", m, m)
  if (!file.exists(name_f)) next
  nt <- fread(cmd = paste("zcat", shQuote(name_f)), select = c("Trans", "name"))
  nt <- nt[grepl("^ENSMUSG|^ENSMUST", Trans)]
  nt[, ens := sub("\\.t[0-9]+$", "", Trans)]
  nt[, Trans_key := Trans]
  map_list[[m]] <- unique(nt[, .(Trans_key, ens, name)])
}
map_all <- unique(rbindlist(map_list))
# direct transcript-level keys cover ENSMUST + ENSMUSG.tN; gene-level keys cover ENSMUSG strips
tx2name  <- map_all[!duplicated(Trans_key)][, setNames(name, Trans_key)]
ens2name <- map_all[!duplicated(ens), setNames(name, ens)]
cat(sprintf("ENSMUSG->name map: %d entries\n", length(ens2name)))
stopifnot(all(TARGETS %in% ens2name))

nc_samples  <- c("NCs1", "NCs2", "NCs3")
d1_samples  <- c("Day1s1", "Day1s2", "Day1s3")
six <- c(nc_samples, d1_samples)

chunk_dirs <- c(m6A = "site_ctx_chunks", m5C = "site_ctx_chunks_m5c",
                PSU = "site_ctx_chunks_psu", inosine = "site_ctx_chunks_inosine",
                Um = "site_ctx_chunks_um", Am = "site_ctx_chunks_am",
                Cm = "site_ctx_chunks_cm", Gm = "site_ctx_chunks_gm")

res <- list()
for (mod in names(chunk_dirs)) {
  ctx_dir <- file.path(data_dir, chunk_dirs[mod])
  if (!dir.exists(ctx_dir)) { cat("skip", mod, "\n"); next }
  ctx <- rbindlist(lapply(six, function(s) {
    f <- file.path(ctx_dir, paste0("ctx_", s, ".tsv.gz"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    dt[, sample := s]
    dt[, grp := ifelse(s %chin% nc_samples, "Sham", "Day1")]
    dt
  }))
  ctx[, gene := tx2name[trans]]
  ctx[is.na(gene), gene := ens2name[sub("\\.t[0-9]+$", "", trans)]]
  ctx <- ctx[gene %chin% TARGETS]
  if (!nrow(ctx)) { cat(mod, ": no target-gene sites\n"); next }
  site_ctx <- ctx[, .(n_det_NC = sum(grp == "Sham"), n_det_Day1 = sum(grp == "Day1"),
                      mean_frac_Sham = mean(fraction[grp == "Sham"]),
                      mean_frac_Day1 = mean(fraction[grp == "Day1"])),
                  by = .(gene, trans, pos)]
  paired <- site_ctx[n_det_NC == 3 & n_det_Day1 == 3]
  per_gene <- paired[, .(n_paired = .N,
                         median_frac_Sham = median(mean_frac_Sham),
                         median_frac_Day1 = median(mean_frac_Day1),
                         median_delta_pp = 100 * median(mean_frac_Day1 - mean_frac_Sham),
                         pct_up = round(100 * mean(mean_frac_Day1 > mean_frac_Sham), 1)),
                     by = gene]
  per_gene[, mod := mod]
  res[[mod]] <- per_gene
  cat(sprintf("%-8s: %d paired target-gene sites across %d genes\n",
              mod, sum(per_gene$n_paired), nrow(per_gene)))
}

out <- rbindlist(res, fill = TRUE)
fwrite(out, out_f, sep = "\t")
print(out[order(gene, mod)])
cat("FIG2F3 PER-GENE DELTAS DONE\n")
