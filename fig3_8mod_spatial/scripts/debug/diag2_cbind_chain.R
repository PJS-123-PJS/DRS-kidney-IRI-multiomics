# Diagnostic 2: full chain — row-subset + cbind of 2 timepoints
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(Matrix) })
setwd("~/scRNA_analysis/validation")
meta_all <- read.table("GSE139107_MouseIRI.metadata.txt.gz", header = TRUE,
                       sep = "\t", stringsAsFactors = FALSE, row.names = 1)

load_tp <- function(tp_name) {
  file <- sprintf("GSE139107_MouseIRI_%s.dge.txt.gz", tp_name)
  hdr <- strsplit(readLines(gzfile(file, "r"), n = 1), "\t")[[1]]; closeAllConnections()
  raw <- fread(file, header = FALSE, skip = 1, sep = "\t")
  gene_ids <- raw[[1]]
  mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(mat) <- gene_ids; rm(raw)
  if (length(hdr) == ncol(mat) + 1) hdr <- hdr[-1]
  colnames(mat) <- hdr
  mat <- Matrix(mat, sparse = TRUE)
  cells_tp <- intersect(colnames(mat), rownames(meta_all))
  cat(sprintf("%s: ncol=%d cells_tp=%d mat_class=%s\n", tp_name, ncol(mat), length(cells_tp), class(mat)[1]))
  mat[, cells_tp, drop = FALSE]
}

counts_all <- list(control = load_tp("control"), `4hours` = load_tp("4hours"))
cat("per-tp colnames null?", sapply(counts_all, function(m) is.null(colnames(m))), "\n")
cat("per-tp colnames len:", sapply(counts_all, function(m) length(colnames(m))), "\n")
cat("gene-id overlap head:", paste(head(rownames(counts_all[[1]]), 3), collapse=" | "), "\n")

all_genes <- Reduce(intersect, lapply(counts_all, rownames))
cat("all_genes:", length(all_genes), "\n")

sub1 <- counts_all[[1]][all_genes, ]
cat("after row-subset: colnames null?", is.null(colnames(sub1)), " len:", length(colnames(sub1)), "\n")

combined_counts <- do.call(cbind, lapply(counts_all, function(m) m[all_genes, ]))
cat("after cbind: dim", nrow(combined_counts), "x", ncol(combined_counts),
    " colnames null?", is.null(colnames(combined_counts)),
    " len:", length(colnames(combined_counts)), "\n")

expected_cells <- unlist(lapply(counts_all, colnames))
cat("expected_cells:", length(expected_cells), " null?", is.null(expected_cells), "\n")
if (!identical(colnames(combined_counts), expected_cells)) {
  colnames(combined_counts) <- expected_cells
  cat("reassigned colnames; null now?", is.null(colnames(combined_counts)), "\n")
}

combined_meta <- do.call(rbind, list(meta_all[colnames(counts_all[[1]]), ],
                                     meta_all[colnames(counts_all[[2]]), ]))
combined_meta <- combined_meta[colnames(combined_counts), ]
cat("meta aligned:", nrow(combined_meta), "\n")
so <- CreateSeuratObject(combined_counts, meta.data = combined_meta, min.cells = 3, min.features = 100)
cat("SEURAT OK:", ncol(so), "cells x", nrow(so), "genes\n")
