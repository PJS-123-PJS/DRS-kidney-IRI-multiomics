# Fig7 Panel G (new): Metabolome × Epitranscriptome — DYNAMIC from Excel, NO hardcoding
# ============================================================================
# v2 fix: All values dynamically extracted from real metabolomics Excel file.
# Lactate and other key metabolites found by name pattern matching — no hardcoded values.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(readxl) })

source("../../shared/palette/color_palette.R")
source("Fig6_config.R")

metab_file <- "/mnt/e/BN2505012AH02S66N1/Metabolomics/result/02.Statistical Analysis/Differentially Expressed Metabolites.xlsx"

cat("Loading real metabolomics data...\n")
nc_d1 <- as.data.table(read_excel(metab_file, sheet="NCvsDay1"))
nc_d3 <- as.data.table(read_excel(metab_file, sheet="NCvsDay3"))
nc_d7 <- as.data.table(read_excel(metab_file, sheet="NCvsDay7"))

theme_nature <- theme_classic(base_size=6.5,base_family="Arial")+
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
    panel.grid=element_blank(),legend.title=element_text(size=6),legend.text=element_text(size=5.5))

# =============================================================================
# EXACT metabolite lookup (v3 fix, 2026-07-14)
# -----------------------------------------------------------------------------
# BUG FIXED: v2 used grepl() SUBSTRING matching on MS2 name, which silently
# grabbed WRONG compounds under clean biological labels, e.g.:
#   "SAM"        -> matched "N-Acetylgalacto[SAM]ine"  (a sugar, NOT S-Adenosylmethionine)
#   "Succinate"  -> matched a synthetic indazole-...-dimethyl[succinic acid]
#   "Citrate"    -> matched "Homo[citric] acid"
#   "Glutathione"-> matched "S-Formyl[glutathione]"
# S-Adenosylmethionine (SAM) is NOT present in this dataset at all — the panel
# must NOT claim to show it. Only compounds whose MS2 name matches EXACTLY are used.
# =============================================================================
find_metabolite <- function(dt, exact_names, label) {
  for(nm in exact_names) {
    hits <- dt[tolower(trimws(dt[["MS2 name"]])) == tolower(nm)]
    if(nrow(hits) > 0) {
      hit <- hits[which.max(abs(hits[["Log_Foldchange"]]))]  # strongest if duplicates
      return(data.table(
        Metabolite=label, log2FC=hit[["Log_Foldchange"]],
        pvalue=hit[["P-Value"]], SuperClass=hit[["Super.Class"]],
        VIP=if("VIP" %in% colnames(hit)) hit[["VIP"]] else NA_real_,
        MS2_name=hit[["MS2 name"]]
      ))
    }
  }
  return(data.table(Metabolite=label, log2FC=NA, pvalue=NA, SuperClass=NA, VIP=NA, MS2_name=NA))
}

# Curated to compounds VERIFIED present by exact MS2-name match (checked against
# real Excel 2026-07-14). SAM itself is absent — the methyl-cycle depletion is
# evidenced by its precursor (Methionine) and product (SAH), both significantly down.
metab_patterns <- list(
  list(patterns=c("Methionine"),              label="Methionine\n(SAM precursor)"),
  list(patterns=c("S-Adenosylhomocysteine"),  label="SAH\n(SAM product)"),
  list(patterns=c("Pyruvate"),                label="Pyruvate\n(glycolysis)"),
  list(patterns=c("Lactate"),                 label="Lactate\n(glycolysis)"),
  list(patterns=c("S-Formylglutathione"),     label="S-Formylglutathione\n(redox)"),
  list(patterns=c("Homocitric acid"),         label="Homocitric acid\n(TCA-related)")
)

cat("\n=== Dynamic metabolite extraction (NC vs Day1) ===\n")
metab_dynamic <- rbindlist(lapply(metab_patterns, function(mp) {
  res <- find_metabolite(nc_d1, mp$patterns, mp$label)
  if(!is.na(res$log2FC)) cat(sprintf("  %-30s log2FC=%+.2f  p=%.4f\n", mp$label, res$log2FC, res$pvalue))
  else cat(sprintf("  %-30s NOT FOUND\n", mp$label))
  res
}), fill=TRUE)

# Remove NA rows
metab_dynamic <- metab_dynamic[!is.na(log2FC)]
metab_dynamic[, significance := fifelse(pvalue < 0.01, "p<0.01", fifelse(pvalue < 0.05, "p<0.05", "n.s."))]
metab_dynamic[, Category := fifelse(grepl("SAM|Meth|SAH", Metabolite), "Methyl Cycle",
                             fifelse(grepl("Glut|Formyl", Metabolite), "Redox", "Energy"))]
metab_dynamic[, Metabolite := factor(Metabolite, levels=rev(Metabolite))]

# QC flag: lactate decrease contradicts IRI biology
lactate_row <- metab_dynamic[grepl("Lactate", Metabolite)]
lactate_flag <- if(nrow(lactate_row)>0 && !is.na(lactate_row$log2FC) && lactate_row$log2FC < 0) {
  sprintf("QC flag: Lactate decrease (log2FC=%+.2f) contradicts IRI glycolysis. Verify data integrity.", lactate_row$log2FC)
} else ""

# =============================================================================
# Panel G1: Key Metabolites Bar Chart (DYNAMIC)
# =============================================================================
cat("\nPanel G1: Key metabolites...\n")
p1 <- ggplot(metab_dynamic, aes(x=log2FC, y=Metabolite, fill=Category))+
  geom_col(width=0.65, alpha=0.88)+
  geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.8)+
  geom_text(aes(label=sprintf("%+.2f (%s)", log2FC, significance),
                hjust=ifelse(log2FC>0, -0.08, 1.08)), size=2.8, family="Arial")+
  scale_fill_manual(values=c("Methyl Cycle"="#D6404C","Redox"="#3D7EA6","Energy"="#E8753A"), name="Pathway")+
  scale_x_continuous(expand=expansion(mult=c(0.42,0.42)))+  # wide L/R margin so value labels on longest bars are not clipped
  labs(x="log2 FC (Day1 vs NC)", y=NULL,
       title="Key Metabolites: Methionine Cycle & Energy (Day1)",
       subtitle=sprintf("Exact MS2-name match from %s. %s", basename(metab_file), lactate_flag))+
  theme_nature+theme(legend.position="bottom")

# =============================================================================
# Panel G2: Super-Class Distribution
# =============================================================================
cat("Panel G2: Super class distribution...\n")
sig_d1 <- nc_d1[`P-Value` < 0.05 & abs(Log_Foldchange) > 1]
sig_d1 <- sig_d1[!is.na(Super.Class) & Super.Class != ""]

sc_summary <- sig_d1[, .(N=.N, Mean_abs_log2FC=mean(abs(Log_Foldchange)),
  Pct_Up=sum(Log_Foldchange>0)/.N*100), by=Super.Class][order(-N)]
top_sc <- sc_summary[1:min(8, nrow(sc_summary))]
top_sc[, Super.Class := factor(Super.Class, levels=rev(Super.Class))]

p2 <- ggplot(top_sc, aes(x=N, y=Super.Class, fill=Pct_Up))+
  geom_col(alpha=0.88, width=0.7)+
  geom_text(aes(label=paste0("N=",N,"  ",round(Pct_Up),"% Up")), hjust=-0.1, size=2.8, family="Arial")+
  scale_fill_gradient2(low="#3D7EA6", mid="grey80", high="#D6404C", midpoint=50, name="% Up")+
  scale_x_continuous(expand=expansion(mult=c(0,0.25)))+
  labs(x="Significant Metabolites (|log2FC|>1, p<0.05)", y=NULL,
       title="Metabolite Super-Class Distribution (NC vs Day1)",
       subtitle=sprintf("%d annotated metabolites significantly changed", nrow(sig_d1)))+
  theme_nature+theme(legend.position="bottom")

# =============================================================================
# Panel G3 REMOVED (2026-07-14): methyl-cycle metabolites (Methionine, SAH) are
# significant ONLY at Day1 (absent from NCvsDay3 / NCvsDay7 differential sheets).
# A "time-course" line with a single real timepoint would be misleading — the
# Day1-acute-only nature IS the finding, already conveyed by the G1 bar chart.
# =============================================================================
n_acute_only <- nrow(nc_d3[tolower(trimws(`MS2 name`)) %in% c("methionine","s-adenosylhomocysteine")]) == 0
cat(sprintf("Panel G3: skipped — methyl-cycle metabolites Day1-acute only (present in Day3 sheet: %s)\n",
            !n_acute_only))

# =============================================================================
# Export — clean 2-subpanel (G1 key metabolites | G2 super-class); no time-course
# =============================================================================
p_combined <- (p1 | p2) +
  plot_annotation(
    title="Metabolome × Epitranscriptome: Methionine-Cycle Depletion (Methionine + SAH down, Day1-acute) Parallels m6A Decrease",
    subtitle=paste0("Source: Differentially Expressed Metabolites.xlsx (n=",format(nrow(nc_d1),big.mark=",")," metabolites). Exact MS2-name match; SAM itself not detected in dataset."),
    theme=theme(plot.title=element_text(size=11,face="bold",family="Arial"),
                plot.subtitle=element_text(size=7.5,color="grey40",family="Arial")))

ggsave(file.path(OUTDIR,"Fig7C_Metabolome_Dynamic.pdf"),p_combined,width=10,height=4.5,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig7C_Metabolome_Dynamic.png"),p_combined,width=10,height=4.5,dpi=400)

# Save dynamic extraction table
fwrite(metab_dynamic, file.path(OUTDIR,"Fig7C_Metabolite_Dynamic.csv"))
cat(sprintf("\nPanel G complete — %d metabolites dynamically extracted, 0 hardcoded\n", nrow(metab_dynamic)))
if(nchar(lactate_flag) > 0) cat(sprintf("WARNING: %s\n", lactate_flag))
