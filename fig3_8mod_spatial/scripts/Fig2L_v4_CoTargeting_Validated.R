# =============================================================================
# Fig2 Panel L v4 — Multi-mod co-targeting on VALIDATED direction sets
# 2026-07-22 — replaces v3 (built on unfiltered trans_diff noise; the m6A
# DOWN transcripts live at 2.4% median stoichiometry). Now: pairwise Fisher
# overlap among validated HC direction sets (5-mod x direction), per-gene
# mod-count distribution, and spotlight of multi-mod hotspot genes.
# Background: union of genes with any differential site call in any caller
# (testable universe from the 8 <mod>_diff_site_context_NEWPATH.tsv.gz).
# Mock-testable: FIG2L_DATA_DIR with upload_5mod_validated/ + bg_genes.txt
# =============================================================================
Sys.setlocale("LC_ALL", "C.utf8")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

dd <- Sys.getenv("FIG2L_DATA_DIR", "../data_placeholder")          # unused default
set_dir <- file.path(dd, "upload_5mod_validated")
ctx_dir <- Sys.getenv("FIG2L_CTX_DIR", file.path(dd))
out_dir <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

show_sets <- c("M6A_UP", "M6A_DOWN", "M5C_UP", "M5C_DOWN",
               "PSU_UP", "INOSINE_UP", "NM_UP", "NM_DOWN")
set_labels <- c(M6A_UP = "m6A↑", M6A_DOWN = "m6A↓(resid.)", M5C_UP = "m5C↑",
                M5C_DOWN = "m5C↓", PSU_UP = "PSU↑", INOSINE_UP = "inosine↑",
                NM_UP = "Nm↑", NM_DOWN = "Nm↓")
# short labels for crowded spotlight tags
set_short <- c(M6A_UP = "m6A↑", M6A_DOWN = "m6A↓r", M5C_UP = "m5C↑",
               M5C_DOWN = "m5C↓", PSU_UP = "PSU↑", INOSINE_UP = "ino↑",
               NM_UP = "Nm↑", NM_DOWN = "Nm↓")

rd <- function(f) { g <- readLines(f); unique(g[!g %in% c("-", "", "NA")]) }
sets <- lapply(show_sets, function(s) rd(file.path(set_dir, paste0(s, ".txt"))))
names(sets) <- show_sets

# ---- background: testable universe ------------------------------------------
bg_file <- file.path(ctx_dir, "bg_genes.txt")
if (file.exists(bg_file)) {
  bg <- rd(bg_file)
} else {
  mods <- c("m6A", "m5c", "psu", "inosine", "um", "am", "cm", "gm")
  bg <- unique(unlist(lapply(mods, function(m) {
    f <- file.path(ctx_dir, sprintf("%s_diff_site_context_NEWPATH.tsv.gz", m))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, select = "name")
    dt[name != "-"]$name
  })))
  writeLines(bg, bg_file)
}
cat(sprintf("background (testable universe): %s genes\n", format(length(bg), big.mark = ",")))
sets <- lapply(sets, intersect, bg)
print(sapply(sets, length))

# ---- pairwise Fisher overlap -------------------------------------------------
K <- length(show_sets)
or_mat <- matrix(NA, K, K, dimnames = list(show_sets, show_sets))
pv_mat <- or_mat
ov_mat <- matrix(0L, K, K, dimnames = list(show_sets, show_sets))
for (i in seq_len(K)) for (j in seq_len(K)) {
  if (i == j) { or_mat[i, j] <- NA; next }
  a <- length(intersect(sets[[i]], sets[[j]]))
  b <- length(setdiff(sets[[i]], sets[[j]]))
  cc <- length(setdiff(sets[[j]], sets[[i]]))
  d <- length(bg) - a - b - cc
  ft <- fisher.test(matrix(c(a, b, cc, d), 2, byrow = TRUE))
  or_mat[i, j] <- unname(ft$estimate)
  pv_mat[i, j] <- ft$p.value
  ov_mat[i, j] <- a
}

or_dt <- as.data.table(or_mat, keep.rownames = "set1")
or_l <- melt(or_dt, id.vars = "set1", variable.name = "set2", value.name = "OR")
ov_dt <- as.data.table(ov_mat, keep.rownames = "set1")
ov_l <- melt(ov_dt, id.vars = "set1", variable.name = "set2", value.name = "overlap")
pv_dt <- as.data.table(pv_mat, keep.rownames = "set1")
pv_l <- melt(pv_dt, id.vars = "set1", variable.name = "set2", value.name = "p")
hm <- merge(merge(or_l, ov_l, by = c("set1", "set2")), pv_l, by = c("set1", "set2"))
hm[, log2OR := log2(pmax(OR, 1e-3))]
hm[, sig := ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
             ifelse(p < 0.05, "*", "ns"))))]
hm[, set1 := factor(set1, levels = show_sets, labels = set_labels[show_sets])]
hm[, set2 := factor(set2, levels = rev(show_sets), labels = set_labels[rev(show_sets)])]

p_hm <- ggplot(hm, aes(x = set1, y = set2, fill = log2OR)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(OR), "", paste0(sprintf("%.1f", log2OR), "\n",
                sig, " (", overlap, ")"))), size = 2.4, lineheight = 0.9) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#B22222",
                       midpoint = 0, na.value = "grey92",
                       guide = guide_colorbar(title = "log2(OR)",
                                              title.position = "top", barwidth = 0.6,
                                              barheight = 5)) +
  labs(x = NULL, y = NULL,
       title = "A  Cross-mod overlap of validated direction sets",
       subtitle = sprintf("Fisher exact, background = %s genes with any diff site call; log2OR / significance / (overlap n)",
                          format(length(bg), big.mark = ","))) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7, color = "grey30"),
        legend.position = "right")

# ---- per-gene validated-mod count -------------------------------------------
up_sets <- c("M6A_UP", "M5C_UP", "PSU_UP", "INOSINE_UP", "NM_UP")
dn_sets <- c("M6A_DOWN", "M5C_DOWN", "NM_DOWN")
mem <- rbindlist(lapply(names(sets), function(s) {
  if (!length(sets[[s]])) return(NULL)
  data.table(gene = sets[[s]], set = s,
             dir = ifelse(s %in% up_sets, "gain-set", "loss-set"))
}))
gene_n <- mem[, .(n_sets = uniqueN(set),
                  dirs = paste(sort(unique(dir)), collapse = "+")),
              by = gene]
cnt <- gene_n[, .N, by = n_sets][order(n_sets)]
cnt[, pct := round(100 * N / sum(N), 1)]

p_cnt <- ggplot(cnt, aes(x = factor(n_sets), y = N)) +
  geom_col(fill = "#6A3D9A", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = paste0(N, "\n(", pct, "%)")), vjust = -0.3, size = 2.7, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Validated direction sets hit per gene",
       y = "Genes",
       title = "B  Single-mod dominates",
       subtitle = sprintf("%s genes in >=1 validated set;\nmulti-set = co-targeting candidates",
                          format(nrow(gene_n), big.mark = ","))) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 8.5),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7, color = "grey30"))

# ---- spotlight: genes in >=3 sets --------------------------------------------
spot <- gene_n[n_sets >= 3][order(-n_sets, gene)]
spot_sets <- mem[gene %in% spot$gene,
                 .(sets = paste(sort(set_short[set]), collapse = " ")), by = gene]
spot <- merge(spot, spot_sets, by = "gene")
fwrite(spot, file.path(out_dir, "Fig2L_v4_multimod_spotlight.tsv"), sep = "\t")
fwrite(hm, file.path(out_dir, "Fig2L_v4_pairwise_overlap.tsv"), sep = "\t")
cat(sprintf("spotlight genes (>=3 sets): %d\n", nrow(spot)))
print(head(spot, 15))

top_lab <- head(spot[order(-n_sets, gene)], 12)
p_spot <- ggplot(top_lab, aes(x = reorder(gene, n_sets), y = n_sets)) +
  geom_col(fill = "#B22222", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sets), hjust = -0.05, size = 2.3, color = "grey25") +
  coord_flip(clip = "off") +
  scale_y_continuous(breaks = 0:8, expand = expansion(mult = c(0, 0.9))) +
  labs(x = NULL, y = "Validated sets hit",
       title = "C  Multi-mod hotspot genes",
       subtitle = "top genes by # validated direction sets (labels = sets)") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 7.5),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7, color = "grey30"),
        plot.margin = margin(5.5, 60, 5.5, 5.5))

p_out <- (p_hm | p_cnt | p_spot) + plot_layout(widths = c(1.7, 0.85, 1.15))
ggsave(file.path(out_dir, "Fig2L_v4_CoTargeting_Validated.pdf"), p_out,
       width = 14.5, height = 4.8, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig2L_v4_CoTargeting_Validated.png"), p_out,
       width = 14.5, height = 4.8, dpi = 600)
cat("saved Fig2L_v4_CoTargeting_Validated.{pdf,png}\n")
