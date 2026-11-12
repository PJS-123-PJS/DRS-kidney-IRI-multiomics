# Fig2E v3 — 5×5 partial correlation from per_site files (Nm merged, no OOM)
# ============================================================================
# Reads per_site.add_name.xls.gz (40MB) instead of all_{mod}.rds (760MB).
# Pairwise merge preserves more transcripts per correlation pair.
# Nm = Am/Cm/Gm/Um pooled. Groups: NC, Day1, Day3, Day7.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

SITE_DIR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
TPM_FILE <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/result/03.expression/transcript/isoforms.TPM.matrix.xls"
OUTDIR   <- "../outputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

BASE_MODS <- c("m6a", "m5c", "psu", "inosine")
NM_MODS   <- c("Am", "Cm", "Gm", "Um")
GROUPS <- c("NC", "IS", "Day1", "Day3", "Day7")
DISPLAY  <- c("m6A", "m5C", "PSU", "Inosine", "Nm")
ALL_KEYS <- c(BASE_MODS, "Nm")
SAMPLE_N <- Inf  # no sampling, use all sites

# ═══════════════════════════════════
# Helper: per-transcript mean fraction from per_site file
# ═══════════════════════════════════
agg_trans <- function(mod, cmp, grp_a, grp_b) {
  f <- file.path(SITE_DIR, mod, "site_diff", cmp,
                 paste0(cmp, ".", mod, ".per_site.add_name.xls.gz"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, select = c("trans", "a_pct_modified", "b_pct_modified"))
  if (nrow(dt) > SAMPLE_N) dt <- dt[sample(.N, SAMPLE_N)]
  rbind(
    dt[, .(Trans = trans, Group = grp_a, Frac = mean(a_pct_modified, na.rm = TRUE)), by = trans],
    dt[, .(Trans = trans, Group = grp_b, Frac = mean(b_pct_modified, na.rm = TRUE)), by = trans]
  )
}

# ═══════════════════════════════════
# Step 1: per-transcript per-group fraction (base mods)
# ═══════════════════════════════════
cat("=== Step 1: Base mods ===\n")
mod_trans <- list()
for (mod in BASE_MODS) {
  cat(sprintf("  %s:", mod))
  parts <- list(
    agg_trans(mod, "NC_vs_Day1", "NC", "Day1"),
    agg_trans(mod, "NC_vs_Day3", "NC", "Day3")[Group == "Day3"],
    agg_trans(mod, "NC_vs_Day7", "NC", "Day7")[Group == "Day7"],
    agg_trans(mod, "NC_vs_I",   "NC", "IS")[Group == "IS"]
  )
  all <- rbindlist(parts[!sapply(parts, is.null)])
  wide <- dcast(all, Trans ~ Group, value.var = "Frac", fun.aggregate = mean, fill = NA)
  avail <- intersect(GROUPS, names(wide))
  setnames(wide, old = avail, new = paste0(mod, "_", avail))
  mod_trans[[mod]] <- wide
  cat(sprintf(" %d tx\n", nrow(wide)))
  rm(parts, all); gc()
}

# ═══════════════════════════════════
# Step 2: Nm merged
# ═══════════════════════════════════
cat("\n=== Step 2: Nm merged ===\n")
nm_parts <- list()
for (mod in NM_MODS) {
  cat(sprintf("  %s", mod))
  nm_parts <- c(nm_parts, list(agg_trans(mod, "NC_vs_Day1", "NC", "Day1")))
  nm_parts <- c(nm_parts, list(agg_trans(mod, "NC_vs_Day3", "NC", "Day3")[Group == "Day3"]))
  nm_parts <- c(nm_parts, list(agg_trans(mod, "NC_vs_Day7", "NC", "Day7")[Group == "Day7"]))
  nm_parts <- c(nm_parts, list(agg_trans(mod, "NC_vs_I",   "NC", "IS")[Group == "IS"]))
}
nm_all <- rbindlist(nm_parts[!sapply(nm_parts, is.null)])
nm_pooled <- nm_all[, .(Frac = mean(Frac, na.rm = TRUE)), by = .(Trans, Group)]
nm_wide <- dcast(nm_pooled, Trans ~ Group, value.var = "Frac", fun.aggregate = mean, fill = NA)
avail <- intersect(GROUPS, names(nm_wide))
setnames(nm_wide, old = avail, new = paste0("Nm_", avail))
mod_trans[["Nm"]] <- nm_wide
cat(sprintf("  -> %d tx\n", nrow(nm_wide)))
rm(nm_parts, nm_all, nm_pooled); gc()

# ═══════════════════════════════════
# Step 3: TPM
# ═══════════════════════════════════
cat("\n=== Step 3: TPM ===\n")
tpm <- fread(TPM_FILE); setnames(tpm, "ID", "Trans")
scols <- grep("^(NCs|Is|Day[137]s)[0-9]$", names(tpm), value = TRUE)
tpm_long <- melt(tpm, id.vars = "Trans", measure.vars = scols,
                 variable.name = "sample", value.name = "TPM")
tpm_long[grepl("^NCs", sample), Group := "NC"]
tpm_long[grepl("^Is", sample), Group := "IS"]
tpm_long[grepl("^Day1s", sample), Group := "Day1"]
tpm_long[grepl("^Day3s", sample), Group := "Day3"]
tpm_long[grepl("^Day7s", sample), Group := "Day7"]
tpm_grp <- tpm_long[!is.na(Group), .(TPM = mean(TPM, na.rm = TRUE)), by = .(Trans, Group)]
tpm_wide <- dcast(tpm_grp, Trans ~ Group, value.var = "TPM")
for (g in intersect(GROUPS, names(tpm_wide)))
  setnames(tpm_wide, g, paste0("TPM_", g))
tpms <- grep("^TPM_", names(tpm_wide), value = TRUE)
tpm_wide[, max_tpm := apply(.SD, 1, max, na.rm = TRUE), .SDcols = tpms]
tpm_top <- tpm_wide[max_tpm > 1][, max_tpm := NULL]
cat(sprintf("  TPM>1: %d tx\n", nrow(tpm_top)))

# ═══════════════════════════════════
# Step 4: Pairwise partial correlation
# ═══════════════════════════════════
cat("\n=== Step 4: Pairwise partial correlation ===\n")
# 2026-07-31 fix (audit of FigS2G): two plotting bugs found —
#   (1) old triangle filter `as.integer(Mod1) < as.integer(factor(Mod2, rev))`
#       kept only i+j<6 → showed 4/10 pairs, ALL Nm pairs + m5C-Inosine +
#       PSU-Inosine silently dropped while the title claimed 5×5.
#   (2) scale limits ±0.5 without oob=squish → |R|>0.5 censored to NA grey
#       (strongest correlations rendered as "missing").
#   Also caption claimed "300K sites sampled" but SAMPLE_N=Inf (all sites).
# Fix: compute all matrices first, dedupe pairs by original index (i<j),
# dynamic symmetric limits from the data, honest caption.
mat_list <- list(); n_list <- list()

for (g in GROUPS) {
  tpm_col <- paste0("TPM_", g)
  if (!tpm_col %in% names(tpm_top)) next
  cat(sprintf("  %s:", g))

  pcor_mat <- matrix(NA, 5, 5, dimnames = list(DISPLAY, DISPLAY))
  n_mat    <- matrix(NA, 5, 5, dimnames = list(DISPLAY, DISPLAY))

  for (i in 1:4) for (j in (i+1):5) {
    ki <- ALL_KEYS[i]; kj <- ALL_KEYS[j]
    ci <- paste0(ki, "_", g); cj <- paste0(kj, "_", g)
    if (is.null(mod_trans[[ki]]) || is.null(mod_trans[[kj]])) next
    if (!ci %in% names(mod_trans[[ki]]) || !cj %in% names(mod_trans[[kj]])) next

    pair <- merge(mod_trans[[ki]][, c("Trans", ci), with = FALSE],
                  mod_trans[[kj]][, c("Trans", cj), with = FALSE], by = "Trans")
    pair <- merge(pair, tpm_top[, c("Trans", tpm_col), with = FALSE], by = "Trans")
    pair <- na.omit(pair)
    if (nrow(pair) < 20) next

    ri <- residuals(lm(pair[[ci]] ~ pair[[tpm_col]]))
    rj <- residuals(lm(pair[[cj]] ~ pair[[tpm_col]]))
    pcor_mat[i,j] <- pcor_mat[j,i] <- cor(ri, rj, method = "spearman")
    n_mat[i,j] <- n_mat[j,i] <- nrow(pair)
  }

  cat(sprintf(" n=%d-%d", min(n_mat, na.rm = TRUE), max(n_mat, na.rm = TRUE)))
  if (!is.na(pcor_mat["m6A","m5C"]))
    cat(sprintf(" m6A-m5C:R=%.3f", pcor_mat["m6A","m5C"]))
  cat("\n")
  mat_list[[g]] <- pcor_mat; n_list[[g]] <- n_mat
}

# export full matrices for auditability (missing pairs' values now on disk)
for (g in names(mat_list))
  fwrite(as.data.table(mat_list[[g]], keep.rownames = "Mod"),
         file.path(OUTDIR, sprintf("Table_PartialCor_%s.csv", g)))

# dynamic symmetric color limits across ALL groups (nothing censored)
r_lim <- ceiling(max(abs(unlist(mat_list)), na.rm = TRUE) * 10) / 10
cat(sprintf("Global |R| max -> scale limit ±%.1f\n", r_lim))

plot_list <- list()
for (g in names(mat_list)) {
  pcor_dt <- as.data.table(as.table(mat_list[[g]]))[!is.na(N)]
  setnames(pcor_dt, c("Mod1","Mod2","R"))
  # keep each pair ONCE by original DISPLAY index (upper triangle, incl. Nm)
  pcor_dt <- pcor_dt[match(Mod1, DISPLAY) < match(Mod2, DISPLAY)]
  pcor_dt[, Mod1 := factor(Mod1, levels = DISPLAY)]
  pcor_dt[, Mod2 := factor(Mod2, levels = rev(DISPLAY))]

  p <- ggplot(pcor_dt, aes(x = Mod1, y = Mod2, fill = R)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", R)), size = 3.8, fontface = "bold", color = "grey20") +
    scale_fill_gradient2(low = "#3B4F8A", mid = "white", high = "#CC5A4C",
                         midpoint = 0, limits = c(-r_lim, r_lim),
                         oob = scales::squish,
                         na.value = "grey90", name = "Partial R") +
    labs(title = g) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8), axis.title = element_blank(),
          panel.grid = element_blank(), legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
  plot_list[[g]] <- p
}

# ═══════════════════════════════════
# Step 5: Export
# ═══════════════════════════════════
composite <- wrap_plots(plot_list, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Cross-modification partial correlation (TPM-controlled, Nm merged)",
    subtitle = "5×5: m6A/m5C/PSU/Inosine + Nm (Am/Cm/Gm/Um pooled). Spearman partial R. Pairwise merge.",
    caption = paste0(
      "Nm merged due to caller cross-confusion (R=0.65-0.75). ",
      "per_site.add_name.xls.gz (all tested sites, no sampling) + TPM. ",
      "Partial R = Spearman cor(residuals of lm(mod ~ TPM)), per-transcript mean fraction, TPM>1, pairwise complete, n>=20. ",
      "Groups: NC,IS,Day1,Day3,Day7. Upper triangle shown (all 10 pairs incl. Nm)."
    ),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7, color = "grey40"),
                  plot.caption = element_text(size = 6, color = "grey55", hjust = 0)))

ggsave(file.path(OUTDIR, "Fig2F_PartialCor_v3_Light.pdf"), composite,
       width = 400, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig2F_PartialCor_v3_Light.png"), composite,
       width = 400, height = 100, units = "mm", dpi = 300)
cat("\nDone! 5×5 matrix saved.\n")
