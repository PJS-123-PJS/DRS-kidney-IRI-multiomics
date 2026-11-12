# build_Fig3B_dedup_rows.R — NC-deduplicated row-level data for Fig3B
# ============================================================================
# 2026-07-29 — Fig3B audit follow-up. User decision: B1 site counts = UNIQUE
# sites; Sham (NC) = union across NC_vs_Day1/Day3/Day7 (no cross-comparison
# double-counting). Day groups unaffected (each appears in one comparison).
#
# Pipeline EXACTLY mirrors build_5mod_fraction_cache_v2.R (filter -> canonical
# all_<mod>.rds join (mult="first") -> valid region -> group pct in (0,1]),
# minus sampling (all groups << 30k, v2 never truncated).
#
# FIDELITY GATE: Day1/3/7 counts per mod MUST equal frac_5mod_cache_v2.rds
# exactly — proves pipeline identity; NC deduped counts must be <= cache NC.
#
# ★ 2026-07-29 GATE CAUGHT A CACHE BUG: cache v2 "Nm" is Um-ONLY — the v2
# builder loops submods Am/Cm/Gm/Um but keys all_frac[[paste("Nm",cmp,grp)]]
# WITHOUT the submod name, so each submod OVERWRITES the previous (last = Um).
# This builder's Nm = true pooled Am+Cm+Gm+Um (deduped). Gate therefore
# validates Nm by comparing the Um-only subset vs cache Nm (must match).
#
# Output: output/fig3b_5mod_dedup_rows.rds (Trans/Pos retained for dedup)
# Slim structure maps cached in data/ (one-time E-drive read, reuse later).
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

OUT_DIR  <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
DATA_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
SITE_DIR <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"

BASE_MODS <- list(M6A = c("m6a","m6a"), M5C = c("m5c","m5c"),
                  PSU = c("psu","psu"), INOSINE = c("inosine","inosine"))
NM_SUBS <- list(Am = "Am", Cm = "Cm", Gm = "Gm", Um = "Um")
CMPS <- list(NC_vs_Day1 = c("NC","Day1"), NC_vs_Day3 = c("NC","Day3"), NC_vs_Day7 = c("NC","Day7"))

# ── slim structure map with D-side cache (OOM-safe in-place slimming) ──
load_struc_map <- function(mod_dir, mod_file) {
  slim <- file.path(DATA_DIR, sprintf("struc_map_slim_%s.rds", mod_file))
  if (file.exists(slim)) {
    cat(sprintf("  [map cache hit] %s\n", basename(slim)))
    r <- readRDS(slim); setDT(r); setkey(r, Trans, Pos); return(r)
  }
  cat(sprintf("  [map build] reading all_%s.rds from E: ...\n", mod_file))
  r <- readRDS(file.path(SITE_DIR, mod_dir, paste0("all_", mod_file, ".rds"))); setDT(r)
  r[, c("Gene","Depth","Fraction","Chr","Motif","sample","group") := NULL]  # in-place
  gc(verbose = FALSE)
  saveRDS(r, slim, compress = FALSE)   # one-time; reuse for future Fig3E etc.
  setkey(r, Trans, Pos)
  r
}

# ── per-mod processing: returns list(NC = deduped rows, Day1/3/7 = rows) ──
process_mod <- function(mod_label, mod_dir, mod_file) {
  sm <- load_struc_map(mod_dir, mod_file)
  out <- list(NC = NULL)
  frac_check <- list()   # a_pct consistency across comparisons for dup NC sites
  for (cmp_name in names(CMPS)) {
    cmp_file <- file.path(SITE_DIR, mod_dir, "site_diff", cmp_name,
                          paste0(cmp_name, ".", mod_file, ".per_site.add_name.xls.gz"))
    if (!file.exists(cmp_file)) { cat(sprintf("  %s %s: MISSING\n", mod_label, cmp_name)); next }
    dt <- fread(cmp_file, select = c("trans","pos","a_pct_modified","b_pct_modified",
                                     "effect_size","map_pvalue","diffType"))
    dt <- dt[pmax(a_pct_modified, b_pct_modified) > 0.7 &
             abs(effect_size) >= 0.1 & map_pvalue < 0.05]
    if (nrow(dt) == 0) next
    dt2 <- sm[dt, on = .(Trans = trans, Pos = pos), mult = "first", nomatch = 0]
    dt2 <- dt2[!is.na(Gene_Struc) & Gene_Struc %in% c("five_prime_UTR","CDS","three_prime_UTR")]
    dt2 <- dt2[!is.na(rel_location) & rel_location >= 0 & rel_location <= 3]
    grp_a <- CMPS[[cmp_name]][1]; grp_b <- CMPS[[cmp_name]][2]   # NC, DayX
    # NC side (collect across comparisons for union dedup)
    nc <- dt2[a_pct_modified > 0 & a_pct_modified <= 1]
    nc_rows <- nc[, .(Trans, Pos, Gene_Struc, rel_location, Fraction = a_pct_modified,
                      effect_size, diffType)]
    out$NC <- rbind(out$NC, nc_rows)
    frac_check[[cmp_name]] <- nc_rows[, .(Trans, Pos, Fraction)]
    # Day side (single comparison per Day group — no dedup needed)
    dd <- dt2[b_pct_modified > 0 & b_pct_modified <= 1]
    out[[grp_b]] <- dd[, .(Trans, Pos, Gene_Struc, rel_location, Fraction = b_pct_modified,
                           effect_size, diffType)]
    cat(sprintf("  %s %s: NC+%d  %s=%d\n", mod_label, cmp_name, nrow(nc_rows), grp_b, nrow(out[[grp_b]])))
    rm(dt, dt2); gc(verbose = FALSE)
  }
  rm(sm); gc(verbose = FALSE)
  # NC union dedup by (Trans,Pos), keep first; report a_pct consistency
  if (!is.null(out$NC)) {
    dup <- duplicated(out$NC, by = c("Trans","Pos"))
    n_dup <- sum(dup)
    # consistency of a_pct across comparisons for duplicated sites
    if (n_dup > 0) {
      frac_check <- rbindlist(frac_check)
      v <- frac_check[, .(spread = max(Fraction) - min(Fraction)), by = .(Trans,Pos)][spread > 1e-9]
      cat(sprintf("  %s NC union: %d -> %d (removed %d cross-cmp dups; a_pct spread>0 sites: %d)\n",
                  mod_label, nrow(out$NC), nrow(out$NC) - n_dup, n_dup, nrow(v)))
    }
    out$NC <- out$NC[!dup]
  }
  out
}

# ── run ──
all_rows <- list()
for (m in names(BASE_MODS)) {
  cat(sprintf("\n=== %s ===\n", m))
  r <- process_mod(m, BASE_MODS[[m]][1], BASE_MODS[[m]][2])
  for (grp in c("NC","Day1","Day3","Day7")) {
    if (!is.null(r[[grp]]) && nrow(r[[grp]]) > 0)
      all_rows[[paste(m, grp)]] <- cbind(Modification = m, Group = grp, r[[grp]])
  }
  rm(r); gc(verbose = FALSE)
}
# Nm = pooled submods, NC deduped across the pool
nm_rows <- list()
for (s in names(NM_SUBS)) {
  cat(sprintf("\n=== Nm sub %s ===\n", s))
  r <- process_mod("Nm", NM_SUBS[[s]], NM_SUBS[[s]])
  for (grp in c("NC","Day1","Day3","Day7")) {
    if (!is.null(r[[grp]]) && nrow(r[[grp]]) > 0)
      nm_rows[[grp]] <- rbind(nm_rows[[grp]], cbind(submod = s, r[[grp]]))
  }
  rm(r); gc(verbose = FALSE)
}
# gate helper: Um-only subset (cache v2's silent definition of "Nm")
um_only <- list()
for (grp in c("NC","Day1","Day3","Day7")) {
  d <- nm_rows[[grp]]
  if (is.null(d)) next
  um_only[[grp]] <- d[submod == "Um"]
  d <- d[!duplicated(d, by = c("Trans","Pos"))]
  cat(sprintf("Nm pooled %s (post cross-submod dedup): %d\n", grp, nrow(d)))
  all_rows[[paste("Nm", grp)]] <- cbind(Modification = "Nm", Group = grp,
                                        d[, setdiff(names(d), "submod"), with = FALSE])
}

res <- rbindlist(all_rows)
res[, Modification := factor(Modification, levels = c("M6A","M5C","PSU","INOSINE","Nm"))]
res[, Group := factor(Group, levels = c("NC","Day1","Day3","Day7"))]

# ── FIDELITY GATE vs cache v2 (Day groups must match exactly) ──
cat("\n=== FIDELITY GATE vs frac_5mod_cache_v2.rds ===\n")
cache <- readRDS(file.path(OUT_DIR, "frac_5mod_cache_v2.rds")); setDT(cache)
cmp_new <- res[, .N, by = .(Modification, Group)]
cmp_old <- cache[, .N, by = .(Modification, Group)]
gate <- merge(cmp_new, cmp_old, by = c("Modification","Group"), suffixes = c("_new","_cache"))
gate[, status := fifelse(Group == "NC",
                         fifelse(N_new <= N_cache, "OK(dedup)", "FAIL(NC grew!)"),
                         fifelse(N_new == N_cache, "OK", "FAIL(Day mismatch!)"))]
# Nm exception: cache v2 Nm = Um-only (overwrite bug) — validate Um subset instead
um_counts <- rbindlist(lapply(names(um_only), function(g)
  data.frame(Group = g, N_um = nrow(um_only[[g]]))))
nm_gate <- merge(um_counts, cmp_old[Modification == "Nm", .(Group, N_cache = N)], by = "Group")
nm_gate[, um_status := fifelse(Group == "NC",
                               fifelse(N_um <= N_cache, "OK(cache==Um-only confirmed)",
                                       "FAIL(Um!=cache NC)"),
                               fifelse(N_um == N_cache, "OK(cache==Um-only confirmed)",
                                       "FAIL(Um!=cache Day)"))]
gate <- gate[Modification != "Nm"]
gate <- rbind(gate,
              nm_gate[, .(Modification = factor("Nm", levels = levels(gate$Modification)),
                          Group = factor(Group, levels = levels(gate$Group)),
                          N_new = N_um, N_cache, status = um_status)])
print(gate)
cat("\nNm pooled (correct Am+Cm+Gm+Um) counts replacing Um-only cache values:\n")
print(cmp_new[Modification == "Nm"])
if (any(grepl("FAIL", gate$status))) stop("FIDELITY GATE FAILED — do not use output")

saveRDS(res, file.path(OUT_DIR, "fig3b_5mod_dedup_rows.rds"))

# ── summary for Fig3B header/caption ──
cat("\n=== DEDUP count matrix (B1 bars) ===\n")
print(dcast(cmp_new, Modification ~ Group, value.var = "N"))
cat("\n=== 3'UTR counts (B2 n= labels) ===\n")
u3 <- res[Gene_Struc == "three_prime_UTR", .N, by = .(Modification, Group)]
print(dcast(u3, Modification ~ Group, value.var = "N"))
cat("\n=== region % by mod ===\n")
for (m in levels(res$Modification)) {
  tb <- round(prop.table(table(res[Modification == m]$Gene_Struc)) * 100, 1)
  cat(sprintf("  %-8s 5'UTR=%.1f%% CDS=%.1f%% 3'UTR=%.1f%%\n",
              m, tb["five_prime_UTR"], tb["CDS"], tb["three_prime_UTR"]))
}
cat(sprintf("\nTOTAL unique differential sites (caption n): %s\n", format(nrow(res), big.mark = ",")))
cat("DONE\n")
