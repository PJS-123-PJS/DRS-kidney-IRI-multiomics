#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
BASE <- "/mnt/f/BN2505012AH02S66N1"
drs <- fread(file.path(BASE,"DirectRNA/result/03.expression/gene/gene.TPM.matrix.xls"))
rna <- fread(file.path(BASE,"RNAseq/Result/03.expression/gene/gene.FPKM.matrix.xls"))
drs <- as.data.frame(drs); rownames(drs) <- drs$ID; drs$ID <- NULL
rna <- as.data.frame(rna); rownames(rna) <- rna$ID; rna$ID <- NULL
common <- intersect(rownames(drs), rownames(rna))
drs <- as.matrix(drs[common,]); rna <- as.matrix(rna[common,])
cat("column order:", paste(colnames(drs), collapse=" "), "\n")
ld <- log2(drs+0.01); lr <- log2(rna+0.01)
ps <- sapply(seq_len(15), function(i) cor(lr[,i], ld[,i]))
for (i in 1:15) cat(sprintf("  %-8s r = %.5f\n", colnames(drs)[i], ps[i]))
cat(sprintf("pooled log2(x+0.01): %.5f\n", cor(as.vector(lr), as.vector(ld))))
cat(sprintf("pooled log2(x+1):    %.5f\n", cor(as.vector(log2(rna+1)), as.vector(log2(drs+1)))))
cat(sprintf("pooled raw:          %.5f\n", cor(as.vector(rna), as.vector(drs))))
# single-sample variants for the sample nearest 0.799
j <- which.min(abs(ps-0.799))
cat(sprintf("nearest sample: %s (%.5f); raw %.5f; log2(x+1) %.5f\n",
    colnames(drs)[j], ps[j], cor(rna[,j], drs[,j]), cor(log2(rna[,j]+1), log2(drs[,j]+1))))
