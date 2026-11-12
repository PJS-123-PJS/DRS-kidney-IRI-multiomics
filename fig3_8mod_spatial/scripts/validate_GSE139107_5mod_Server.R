# ============================================================
# validate_GSE139107_5mod.R — PT-segment validation, VALIDATED sets
# 2026-07-22 — mirrors validate_GSE139107_v2.R scoring EXACTLY
# (mean(colMeans) per cell type -> scale per set), but gene sets are
# the validated HC direction sets (no wave framing):
#   M6A_CDS_UP / M6A_CDS_DOWN(residual) / M5C_CDS_DOWN / M5C_CDS_UP
#   INOSINE_UP / PSU_UP / NM_UP / NM_DOWN
# Run: nohup Rscript validate_GSE139107_5mod.R > val5mod.log 2>&1 &
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
})
setwd("~/scRNA_analysis/validation")
res_dir <- "../results"

so_file <- file.path(res_dir, "seurat_merged_validation.rds")
if (file.exists(so_file)) {
  cat("Loading cached Seurat object...\n")
  so <- readRDS(so_file)
} else {
  # mirror v2: rebuild from DGE matrices (NormalizeData only — scoring uses RNA/data)
  cat("Rebuilding Seurat from DGE matrices...\n")
  meta_all <- read.table("GSE139107_MouseIRI.metadata.txt.gz", header = TRUE,
                         sep = "\t", stringsAsFactors = FALSE, row.names = 1)
  timepoints <- c("control"="Control", "4hours"="4h", "12hours"="12h",
                  "2days"="2d", "14days"="14d", "6weeks"="6w")
  counts_all <- list(); cell_meta <- list()
  for (tp_name in names(timepoints)) {
    cat(sprintf("Loading %s...\n", tp_name))
    file <- sprintf("GSE139107_MouseIRI_%s.dge.txt.gz", tp_name)
    hdr <- strsplit(readLines(gzfile(file, "r"), n = 1), "\t")[[1]]; closeAllConnections()
    # fread: ~10-50x faster than read.table on the 660M-field DGE matrices
    raw <- fread(file, header = FALSE, skip = 1, sep = "\t")
    gene_ids <- raw[[1]]
    mat <- as.matrix(raw[, -1, with = FALSE])
    rownames(mat) <- gene_ids
    rm(raw); gc(verbose = FALSE)
    # header has no gene-column label in this GEO release; tolerate either layout
    if (length(hdr) == ncol(mat) + 1) hdr <- hdr[-1]
    stopifnot(length(hdr) == ncol(mat), !anyDuplicated(gene_ids))
    colnames(mat) <- hdr; mat <- Matrix(mat, sparse = TRUE)
    cells_tp <- intersect(colnames(mat), rownames(meta_all))
    cat(sprintf("  %s: %s cells, %s genes, %s in metadata\n", tp_name,
                format(ncol(mat), big.mark=","), format(nrow(mat), big.mark=","),
                format(length(cells_tp), big.mark=",")))
    counts_all[[tp_name]] <- mat[, cells_tp, drop = FALSE]
    cell_meta[[tp_name]] <- meta_all[cells_tp, ]
  }
  all_genes <- Reduce(intersect, lapply(counts_all, rownames))
  # do.call(cbind, named_list) prefixes colnames ("control.<barcode>") AND
  # unlist(named_list) prefixes expected_cells identically — masking the bug.
  # unname() BOTH: matrices before cbind, and use.names=FALSE for expected.
  sub_mats <- lapply(counts_all, function(m) m[all_genes, ])
  combined_counts <- do.call(cbind, unname(sub_mats))
  expected_cells <- unlist(lapply(counts_all, colnames), use.names = FALSE)
  stopifnot(length(expected_cells) == ncol(combined_counts))
  if (!identical(colnames(combined_counts), expected_cells))
    colnames(combined_counts) <- expected_cells
  combined_meta <- do.call(rbind, unname(cell_meta))   # unname: rbind prefixes rownames too
  # hard assertion: positional identity expected; repair by name or die loudly
  if (!identical(rownames(combined_meta), colnames(combined_counts))) {
    if (!all(colnames(combined_counts) %in% rownames(combined_meta)))
      stop("cell-name mismatch between counts and metadata (rbind/cbind prefix bug?)")
    combined_meta <- combined_meta[colnames(combined_counts), ]
  }
  # safety net: drop any unassigned-celltype cells so they cannot form a
  # spurious "" level passing the >50-cell filter downstream (none expected
  # in the official Kirita metadata, verified intact 2026-07-23, md5 7937f622)
  ok_ct <- !is.na(combined_meta$celltype) & combined_meta$celltype != ""
  cat(sprintf("Dropping %d unassigned-celltype cells\n", sum(!ok_ct)))
  if (all(!ok_ct)) stop("ALL cells dropped — metadata misalignment (see debug/ diag scripts)")
  combined_counts <- combined_counts[, ok_ct, drop = FALSE]
  combined_meta <- combined_meta[ok_ct, , drop = FALSE]
  so <- CreateSeuratObject(combined_counts, meta.data = combined_meta,
                           min.cells = 3, min.features = 100)
  # min.features drops cells at object creation -> combined_meta has MORE rows
  # than ncol(so); Seurat v5 also requires NAMED vectors for $<- metadata.
  # Subset+align by surviving cell names, name every vector, factor celltype.
  keep <- colnames(so)
  so$timepoint <- setNames(combined_meta[keep, "Group"], keep)
  so$celltype <- factor(setNames(combined_meta[keep, "celltype"], keep))
  so <- NormalizeData(so)
  saveRDS(so, so_file)
}
cat(sprintf("Seurat: %d cells x %d genes\n", ncol(so), nrow(so)))
cat("timepoint values:", paste(unique(so$timepoint), collapse = ", "), "\n")

map_data <- read.table("../data/ShamA/features.tsv.gz", header = FALSE, sep = "\t",
                       col.names = c("ensembl", "symbol", "type"))
ens2sym <- setNames(map_data$symbol, map_data$ensembl)

load_geneset <- function(file_path, desc) {
  if (!file.exists(file_path)) { cat(sprintf("  SKIP %s\n", desc)); return(character(0)) }
  genes_raw <- readLines(file_path)
  genes_raw <- genes_raw[!grepl("^#|^$", genes_raw)]
  if (length(genes_raw) == 0) return(character(0))
  genes_sym <- if (grepl("^ENSMUSG", genes_raw[1])) na.omit(ens2sym[genes_raw]) else genes_raw
  gs <- intersect(genes_sym, rownames(so))
  cat(sprintf("  %-25s: %d genes\n", desc, length(gs)))
  gs
}

cat("\n=== Loading validated gene sets ===\n")
gd <- "../gene_sets_5mod"
gene_sets <- list(
  "M6A_CDS_UP"   = load_geneset(file.path(gd, "M6A_CDS_UP.txt"),   "m6A CDS Up (validated)"),
  "M6A_CDS_DOWN" = load_geneset(file.path(gd, "M6A_CDS_DOWN.txt"), "m6A CDS Down (residual)"),
  "M5C_CDS_DOWN" = load_geneset(file.path(gd, "M5C_CDS_DOWN.txt"), "m5C CDS Down (validated)"),
  "M5C_CDS_UP"   = load_geneset(file.path(gd, "M5C_CDS_UP.txt"),   "m5C CDS Up"),
  "INOSINE_UP"   = load_geneset(file.path(gd, "INOSINE_UP.txt"),   "inosine Up"),
  "PSU_UP"       = load_geneset(file.path(gd, "PSU_UP.txt"),       "PSU Up"),
  "NM_UP"        = load_geneset(file.path(gd, "NM_UP.txt"),        "Nm Up (union)"),
  "NM_DOWN"      = load_geneset(file.path(gd, "NM_DOWN.txt"),      "Nm Down (union)")
)

cat("\n=== Cell-type enrichment (GSE139107) ===\n")
Idents(so) <- so$celltype
ct_list <- levels(so$celltype)
ct_list <- ct_list[table(so$celltype)[ct_list] > 50]
cat(sprintf("Cell types with >50 cells: %d\n", length(ct_list)))

mat <- GetAssayData(so, assay = "RNA", layer = "data")
ct_idx <- lapply(ct_list, function(ct) WhichCells(so, idents = ct))
names(ct_idx) <- ct_list

enrich_raw <- sapply(gene_sets, function(gs) {
  if (length(gs) == 0) return(rep(NA_real_, length(ct_list)))
  vapply(ct_idx, function(cells) mean(colMeans(mat[gs, cells, drop = FALSE])), numeric(1))
})
rownames(enrich_raw) <- ct_list
enrich_z <- scale(enrich_raw)
cat("\n=== GSE139107 validated-set Z-scores ===\n")
print(round(enrich_z, 2))
write.csv(round(enrich_z, 4), file.path(res_dir, "validation_GSE139107_5mod_zscore.csv"))
write.csv(round(enrich_raw, 6), file.path(res_dir, "validation_GSE139107_5mod_meanexpr.csv"))

# per-timepoint breakdown for PT subtypes (control vs 12hours vs 2days = our Sham/Day1-ish/Day3-ish)
# Group values in metadata: Control / 4hours / 12hours / 2days / 14days / 6weeks
cat("\n=== Per-timepoint PT-subtype scoring ===\n")
pt_cts <- intersect(c("PTS1", "PTS2", "PTS3", "InjuredPT", "NewPT1", "NewPT2"), ct_list)
cat("PT subtypes present:", paste(pt_cts, collapse = ", "), "\n")
tp_list <- c("Control", "12hours", "2days", "14days")
tp_z_all <- list()
for (tp in tp_list) {
  cells_tp <- WhichCells(so, cells = colnames(so)[so$timepoint == tp])
  if (length(cells_tp) < 100) next
  raw_tp <- sapply(gene_sets, function(gs) {
    if (length(gs) == 0) return(rep(NA_real_, length(pt_cts)))
    vapply(pt_cts, function(ct) {
      cells <- intersect(cells_tp, ct_idx[[ct]])
      if (length(cells) < 20) return(NA_real_)
      mean(colMeans(mat[gs, cells, drop = FALSE]))
    }, numeric(1))
  })
  rownames(raw_tp) <- pt_cts
  tp_z_all[[tp]] <- scale(raw_tp)
}
for (tp in names(tp_z_all)) {
  out <- as.data.table(tp_z_all[[tp]], keep.rownames = "CellType")
  out[, Timepoint := tp]
  fwrite(out, file.path(res_dir, paste0("validation_GSE139107_5mod_zscore_", tp, ".csv")))
}
cat("== GSE139107 5-MOD VALIDATION DONE ==\n")
