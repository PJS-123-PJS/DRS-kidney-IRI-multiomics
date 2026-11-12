# Minimal diagnostic: where do colnames vanish?
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
setwd("~/scRNA_analysis/validation")
meta_all <- read.table("GSE139107_MouseIRI.metadata.txt.gz", header = TRUE,
                       sep = "\t", stringsAsFactors = FALSE, row.names = 1)
cat("meta dim:", dim(meta_all), " cols:", paste(names(meta_all), collapse=","), "\n")
cat("meta rownames head:", paste(head(rownames(meta_all),2), collapse=" | "), "\n")
cat("celltype class:", class(meta_all$celltype), " empty:", sum(meta_all$celltype=="", na.rm=TRUE), " NA:", sum(is.na(meta_all$celltype)), "\n")

file <- "GSE139107_MouseIRI_control.dge.txt.gz"
hdr <- strsplit(readLines(gzfile(file, "r"), n = 1), "\t")[[1]]; closeAllConnections()
mat <- as.matrix(read.table(file, header = FALSE, skip = 1, row.names = 1))
cat("hdr len:", length(hdr), " mat dim:", nrow(mat), "x", ncol(mat), "\n")
cat("hdr head:", paste(head(hdr,2), collapse=" | "), "\n")
if (length(hdr) == ncol(mat) + 1) hdr <- hdr[-1]
colnames(mat) <- hdr; mat <- Matrix(mat, sparse = TRUE)
cells_tp <- intersect(colnames(mat), rownames(meta_all))
cat("cells_tp:", length(cells_tp), "of", ncol(mat), "\n")
m1 <- mat[, cells_tp, drop = FALSE]
cat("m1 colnames null?", is.null(colnames(m1)), " len:", length(colnames(m1)), "\n")

combined_counts <- m1[head(rownames(m1), 2000), ]
combined_meta <- meta_all[colnames(combined_counts), ]
cat("cb meta dim:", dim(combined_meta), " celltype NA:", sum(is.na(combined_meta$celltype)), "\n")
ok_ct <- !is.na(combined_meta$celltype) & combined_meta$celltype != ""
cat("ok_ct TRUE:", sum(ok_ct), " FALSE:", sum(!ok_ct), "\n")
cc <- combined_counts[, ok_ct, drop = FALSE]
cat("cc dim:", nrow(cc), "x", ncol(cc), " colnames null?", is.null(colnames(cc)), " len:", length(colnames(cc)), "\n")
so <- CreateSeuratObject(cc, meta.data = combined_meta[ok_ct, , drop=FALSE], min.cells = 3, min.features = 100)
cat("SEURAT OK:", ncol(so), "cells x", nrow(so), "genes\n")
keep <- colnames(so)
so$timepoint <- setNames(combined_meta[keep, "Group"], keep)
so$celltype <- factor(setNames(combined_meta[keep, "celltype"], keep))
cat("metadata assigned OK; timepoints:", paste(unique(so$timepoint), collapse=","), "\n")
