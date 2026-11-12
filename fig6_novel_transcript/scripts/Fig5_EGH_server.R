# Fig5 Panels E, G, H — Server-side computation
# =============================================================================
# Panel E: Real m6A/m5C modification density per class code
# Panel G: Class code × cell type enrichment (scRNA-seq)
# Panel H: Temporal specificity × cell type enrichment
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(dplyr)
})

setwd("~/fig5_analysis")
dir.create("results", FALSE)

# =============================================================================
# ── Panel E: Modification Density ───────────────────────────────────────────
# =============================================================================
cat("\n=== Panel E: Modification Density ===\n")

# Load m6A gene summary (3'UTR data)
m6a_3utr <- fread("zcat < ~/m6a_gene_3utr_summary.csv.gz")
setnames(m6a_3utr, c("Gene","n_sites_Day1","n_sites_NC","mean_frac_Day1","mean_frac_NC","delta"))

# Also check CDS and 5'UTR
m6a_cds <- fread("zcat < ~/m6a_gene_cds_summary.csv.gz")
setnames(m6a_cds, c("Gene","n_sites_Day1","n_sites_NC","mean_frac_Day1","mean_frac_NC","delta"))

m6a_5utr <- fread("zcat < ~/m6a_gene_5utr_summary.csv.gz")
setnames(m6a_5utr, c("Gene","n_sites_Day1","n_sites_NC","mean_frac_Day1","mean_frac_NC","delta"))

cat(sprintf("m6A 3'UTR: %d genes, CDS: %d, 5'UTR: %d\n",
            nrow(m6a_3utr), nrow(m6a_cds), nrow(m6a_5utr)))

# Build ENSEMBL → gene symbol mapping
cat("Building gene symbol mapping...\n")
map_data <- read.table("~/scRNA_analysis/data/ShamA/features.tsv.gz",
                        header=FALSE, sep="\t", col.names=c("ensembl","symbol","type"))
ens2sym <- setNames(map_data$symbol, map_data$ensembl)

# Map m6A ENSEMBL genes to symbols
m6a_3utr[, Symbol := ens2sym[Gene]]
m6a_cds[, Symbol := ens2sym[Gene]]
m6a_5utr[, Symbol := ens2sym[Gene]]

# Load per-class-code gene symbol sets
class_sets <- list()
for (cc in c("k","j","m","u","x","i","n")) {
  f <- sprintf("data/geneset_class_%s_symbols.txt", cc)
  if (file.exists(f)) {
    class_sets[[cc]] <- scan(f, what="")
    cat(sprintf("  Class %s: %d symbols\n", cc, length(class_sets[[cc]])))
  } else {
    class_sets[[cc]] <- character(0)
    cat(sprintf("  Class %s: file not found, skipping\n", cc))
  }
}
all_novel_genes <- scan("data/geneset_Novel_HostGenes_symbols.txt", what="")
cat(sprintf("  All Novel Host Genes: %d symbols\n", length(all_novel_genes)))

# Compute m6A sites per gene for known vs novel categories
# Known: ENSEMBL genes NOT in novel host gene set
# Novel: genes in novel host gene set (by class code)
compute_density <- function(gene_set, label, m6a_data) {
  genes_in_data <- intersect(gene_set, m6a_data$Symbol)
  total_sites <- sum(m6a_data[Symbol %in% genes_in_data, n_sites_NC], na.rm=TRUE)
  # Approximate total gene length from mean_frac (fraction of reads modified)
  # Use n_sites / mean_frac as proxy for gene length
  subset_data <- m6a_data[Symbol %in% genes_in_data]
  subset_data <- subset_data[mean_frac_NC > 0 & n_sites_NC > 0]
  if (nrow(subset_data) == 0) return(data.table(Category=label, N_Genes=0,
         Total_Sites=0, Sites_per_Gene=0, Mean_Frac_NC=0, Mean_Frac_Day1=0))

  data.table(
    Category = label,
    N_Genes = length(genes_in_data),
    Total_Sites = sum(subset_data$n_sites_NC),
    Sites_per_Gene = mean(subset_data$n_sites_NC),
    Mean_Frac_NC = mean(subset_data$mean_frac_NC),
    Mean_Frac_Day1 = mean(subset_data$mean_frac_Day1)
  )
}

# Known mRNA genes (all ENSEMBL genes minus novel host genes)
all_ensembl_genes <- unique(c(m6a_3utr$Symbol, m6a_cds$Symbol, m6a_5utr$Symbol))
all_ensembl_genes <- na.omit(all_ensembl_genes)
known_genes <- setdiff(all_ensembl_genes, all_novel_genes)
cat(sprintf("\nKnown mRNA genes (ENSEMBL - novel): %d\n", length(known_genes)))

# Compute density for each category (using 3'UTR m6A as primary metric)
e_results <- rbindlist(list(
  compute_density(known_genes, "Known mRNA", m6a_3utr),
  compute_density(class_sets[["k"]], "Novel (k)", m6a_3utr),
  compute_density(class_sets[["j"]], "Novel (j)", m6a_3utr),
  compute_density(class_sets[["m"]], "Novel (m)", m6a_3utr),
  compute_density(class_sets[["u"]], "Novel (u)", m6a_3utr),
  compute_density(class_sets[["x"]], "Novel (x)", m6a_3utr),
  compute_density(class_sets[["i"]], "Novel (i)", m6a_3utr)
))

# Also compute using CDS m6A
e_cds <- rbindlist(list(
  compute_density(known_genes, "Known mRNA", m6a_cds),
  compute_density(class_sets[["k"]], "Novel (k)", m6a_cds),
  compute_density(class_sets[["j"]], "Novel (j)", m6a_cds),
  compute_density(class_sets[["m"]], "Novel (m)", m6a_cds)
))

cat("\nPanel E — m6A Modification Density (3'UTR):\n")
print(e_results)

fwrite(e_results, "results/PanelE_m6A_density_3UTR.csv")
fwrite(e_cds, "results/PanelE_m6A_density_CDS.csv")

# =============================================================================
# ── Panel G: Class Code × Cell Type Enrichment ──────────────────────────────
# =============================================================================
cat("\n=== Panel G: Class Code × Cell Type ===\n")

cat("Loading Seurat object...\n")
so <- readRDS("~/scRNA_analysis/results/seurat_merged_annotated.rds")
DefaultAssay(so) <- "RNA"
Idents(so) <- so$celltype
ct_list <- levels(so$celltype)

# For each class code with sufficient genes, compute mean expression per cell type
class_for_enrich <- c("k","j","m")
g_results <- data.table(CellType=ct_list)

for (cc in class_for_enrich) {
  genes <- intersect(class_sets[[cc]], rownames(so))
  cat(sprintf("  Class %s: %d genes in scRNA-seq\n", cc, length(genes)))

  if (length(genes) >= 20) {
    avg_expr <- AverageExpression(so, assays="RNA", slot="data",
                                   features=genes, group.by="celltype")$RNA
    mean_expr <- colMeans(as.matrix(avg_expr))
    zscore <- as.numeric(scale(mean_expr))
    g_results[[paste0("Class_", cc, "_Z")]] <- zscore
  }
}

# Also compute for All Novel
novel_in_sc <- intersect(all_novel_genes, rownames(so))
cat(sprintf("  All Novel: %d genes in scRNA-seq\n", length(novel_in_sc)))
if (length(novel_in_sc) >= 20) {
  avg_expr <- AverageExpression(so, assays="RNA", slot="data",
                                 features=novel_in_sc, group.by="celltype")$RNA
  mean_expr <- colMeans(as.matrix(avg_expr))
  g_results[["All_Novel_Z"]] <- as.numeric(scale(mean_expr))
}

cat("\nPanel G — Class Code × Cell Type Z-scores:\n")
print(g_results)
fwrite(g_results, "results/PanelG_ClassCode_CellType.csv")

# =============================================================================
# ── Panel H: Temporal Specificity × Cell Type ───────────────────────────────
# =============================================================================
cat("\n=== Panel H: Temporal × Cell Type ===\n")

# Use the novel transcript expression data to classify D1-responsive genes
expr_data <- fread("~/fig5_analysis/data/novel_transcript_expression.csv")
cat(sprintf("Expression data: %d transcripts\n", nrow(expr_data)))

# Classify by fold change
expr_data[, FC := Day1_mean / (NC_mean + 0.01)]
expr_data[, Category := "Stable"]
expr_data[FC > 2 & Day1_mean > 5, Category := "D1_Up"]
expr_data[FC < 0.5 & NC_mean > 5, Category := "D1_Down"]

cat("Temporal categories:\n")
print(expr_data[, .N, by=Category])

# Map transcripts to gene symbols (using the expression data ID → class code → symbol)
# For each temporal category, get gene symbols
get_symbols_for_category <- function(cat, expr_dt) {
  tx_ids <- expr_dt[Category == cat, ID]
  # These are transcript IDs (novelXXXX.tX or MGP_...)
  # Extract base gene IDs
  gene_ids <- unique(gsub("\\..*", "", tx_ids))
  # Map to symbols
  symbols <- intersect(gene_ids, all_novel_genes)
  # If gene_ids are MGP format, try direct symbol lookup
  if (length(symbols) < 10) {
    symbols <- intersect(gene_ids, c(class_sets[["k"]], class_sets[["j"]],
                                      class_sets[["m"]], class_sets[["u"]],
                                      class_sets[["x"]], class_sets[["i"]],
                                      class_sets[["n"]]))
  }
  symbols
}

d1_up_genes <- get_symbols_for_category("D1_Up", expr_data)
d1_down_genes <- get_symbols_for_category("D1_Down", expr_data)
stable_genes <- get_symbols_for_category("Stable", expr_data)

cat(sprintf("D1_Up genes: %d, D1_Down: %d, Stable: %d\n",
            length(d1_up_genes), length(d1_down_genes), length(stable_genes)))

# Cell type enrichment for temporal categories
h_results <- data.table(CellType=ct_list)

for (cat_name in c("D1_Up","D1_Down","Stable")) {
  genes <- switch(cat_name,
    "D1_Up" = intersect(d1_up_genes, rownames(so)),
    "D1_Down" = intersect(d1_down_genes, rownames(so)),
    "Stable" = intersect(stable_genes, rownames(so))
  )
  cat(sprintf("  %s: %d genes in scRNA-seq\n", cat_name, length(genes)))

  if (length(genes) >= 10) {
    avg_expr <- AverageExpression(so, assays="RNA", slot="data",
                                   features=genes, group.by="celltype")$RNA
    mean_expr <- colMeans(as.matrix(avg_expr))
    h_results[[paste0(cat_name, "_Z")]] <- as.numeric(scale(mean_expr))
  }
}

cat("\nPanel H — Temporal × Cell Type Z-scores:\n")
print(h_results)
fwrite(h_results, "results/PanelH_Temporal_CellType.csv")

# Also save temporal category counts
temporal_summary <- expr_data[, .(
  N_Transcripts = .N,
  Mean_FC = mean(FC),
  N_Up = sum(Category == "D1_Up"),
  N_Down = sum(Category == "D1_Down"),
  N_Stable = sum(Category == "Stable")
)]
fwrite(temporal_summary, "results/PanelH_Temporal_Summary.csv")

# =============================================================================
cat("\n=== All panels complete ===\n")
cat("Results saved to ~/fig5_analysis/results/\n")
cat("Files: PanelE_m6A_density_*.csv, PanelG_ClassCode_CellType.csv, PanelH_*.csv\n")
