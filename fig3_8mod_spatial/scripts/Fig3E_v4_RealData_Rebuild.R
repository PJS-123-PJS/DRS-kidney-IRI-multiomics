# Fig3 Panel E v5 — integrated summary, DAY1-FOCUSED (user directive 2026-07-29)
# ============================================================================
# v5: summary basis = IRI Day1 response (NC_vs_Day1 differential sites), NOT
#   the 4-group pool. v4 pooled NC+Day1+Day3+Day7 in p_sites/p_region/p_table,
#   making the summary Sham-diluted (and under v2 cache, NC was even x3
#   duplicated). User: "这种summary不要以NC为主" — the paper showcases the
#   Day1 IRI model. p_ct stays as-is: it already uses validated Day1 direction
#   sets (M6A_UP/M5C_DOWN/PSU_UP/INOSINE_UP/NM_UP).
#   Also: reads frac_5mod_cache_v3.rds (true pooled Nm; v2 was Um-only);
#   palette harmonized to the Fig3 canonical MOD_COLORS (Fig3A/B);
#   3'UTR flag: "* = low n (<20)" — no mod has 0 Day1 3'UTR sites (v1-cache
#   era "** no 3UTR sites" text is retired; user already removed it in layout).
# v4 lineage: real-data-only rebuild; HARD FAIL if enrichment CSV missing.
# Env override: FIG3E_ENRICH_CSV
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales) })

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
ENRICH_CSV <- Sys.getenv("FIG3E_ENRICH_CSV",
  "/mnt/d/desktop/DRS/DRS_TechShowcase/fig2_8mod_temporal/results/enrichment_5mod_zscore.csv")
FDR_CSV <- sub("zscore", "permFDR_specific", ENRICH_CSV)
if (!file.exists(ENRICH_CSV)) stop("enrichment_5mod_zscore.csv not found — download from server first (no fallback allowed)")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

theme_nature <- theme_classic(base_size=6.5,base_family="Arial") +
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
    axis.ticks.length=unit(1,"pt"),legend.title=element_text(size=6),
    legend.text=element_text(size=5.5),legend.key.size=unit(8,"pt"),
    strip.text=element_text(size=6.5,face="bold"),strip.background=element_blank(),
    plot.title=element_text(size=7.5,face="bold",hjust=0.5),
    plot.subtitle=element_text(size=6,color="grey30",hjust=0.5),
    panel.grid=element_blank(),panel.spacing=unit(1.5,"mm"))

# Fig3 canonical palette (Fig3A/B)
MOD_COLORS <- c("M6A"="#D6404C","M5C"="#E8753A","PSU"="#3D7EA6","INOSINE"="#7B5294","Nm"="#2E7D6F")
REGION_FILL <- c("5'UTR"="#F2E8D0","CDS"="#B8D8B8","3'UTR"="#6090C0")
MOD_ORDER <- c("M6A","M5C","PSU","INOSINE","Nm")
# validated headline direction set per mod class (Day1 response sets)
SET_MAP <- c("M6A"="M6A_UP","M5C"="M5C_DOWN","PSU"="PSU_UP","INOSINE"="INOSINE_UP","Nm"="NM_UP")

frac_all <- readRDS(file.path(OUT_DIR,"frac_5mod_cache_v3.rds")); setDT(frac_all)
frac_all[, ModLabel:=factor(Modification,levels=MOD_ORDER)]
frac <- frac_all[Group=="Day1"]                       # ← Day1 basis for site metrics
cat(sprintf("Day1 basis: %s sites\n", format(nrow(frac),big.mark=",")))

# ── 1. Day1 site count (with power tier legend) ──
site_tot <- frac[, .(N=.N), by=ModLabel]
site_tot[, Tier:=fifelse(N>=10000,"★★★★★",fifelse(N>=5000,"★★★★",
  fifelse(N>=1000,"★★★",fifelse(N>=500,"★★","★"))))]
site_tot[, TierLabel:=paste0(format(N,big.mark=",")," (",Tier,")")]

p_sites <- ggplot(site_tot, aes(x=ModLabel,y=N,fill=ModLabel)) +
  geom_col(width=0.5) +
  geom_text(aes(label=TierLabel),vjust=-0.3,size=2.2,family="Arial",lineheight=0.9) +
  scale_fill_manual(values=MOD_COLORS,guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.3)),labels=comma) +
  labs(x=NULL,y="Sites",title="Day1 Differential Sites (vs Sham)",
    subtitle="★=statistical power (3≥1K, 2≥500, 1<500)") +
  theme_nature + theme(plot.title=element_text(hjust=0.5))

# ── 2. Day1 regional distribution ──
reg_sum <- frac[, .N, by=.(ModLabel,Gene_Struc)]
reg_sum[, Proportion:=N/sum(N),by=ModLabel]
reg_sum[, Region:=factor(Gene_Struc,levels=c("three_prime_UTR","CDS","five_prime_UTR"),
  labels=c("3'UTR","CDS","5'UTR"))]
utr3_n <- reg_sum[Region=="3'UTR",.(ModLabel,N_3UTR=N)]
reg_sum <- merge(reg_sum,utr3_n,by="ModLabel",all.x=TRUE)
reg_sum[is.na(N_3UTR),N_3UTR:=0]
reg_sum[, UTR3_flag:=fifelse(N_3UTR>=20,"","*")]
reg_sum[, UTR3_anno:=paste0("3UTR:",N_3UTR,UTR3_flag)]

p_region <- ggplot(reg_sum,aes(x=ModLabel,y=Proportion,fill=Region)) +
  geom_col(width=0.5,color="white",linewidth=0.2) +
  geom_text(aes(label=ifelse(Proportion>0.04,percent(Proportion,accuracy=1),""),color=Region),
    position=position_stack(vjust=0.5),size=2.8,family="Arial") +
  geom_text(aes(y=1.03,label=UTR3_anno),size=2.1,color="grey40",family="Arial") +
  scale_fill_manual(values=REGION_FILL,name=NULL) +
  scale_color_manual(values=c("5'UTR"="grey30","CDS"="grey30","3'UTR"="white"),guide="none") +
  scale_y_continuous(labels=percent_format(),expand=expansion(mult=c(0,0.12))) +
  labs(x=NULL,y="Proportion",title="Regional Distribution at Day1",
    subtitle="* = low 3'UTR n (<20 sites)") +
  theme_nature + theme(plot.title=element_text(hjust=0.5),legend.position="bottom",
    legend.margin=margin(t=-3))

# ── 3. Cell-type enrichment heatmap — REAL server z-scores (Day1 validated sets) ──
ez <- fread(ENRICH_CSV); setnames(ez, 1, "CellType")
CT_ORDER <- intersect(c("PT","TAL","DCT","CNT_CD","ICB","ICA","Podo","Endo","Stromal",
                        "Macrophage","T_NK","B_cell","Neutrophil","Proliferating"), ez$CellType)
miss <- setdiff(SET_MAP, names(ez))
if (length(miss)) stop(paste("sets missing in enrichment csv:", paste(miss, collapse=", ")))
ct_long <- melt(ez[CellType %in% CT_ORDER, c("CellType", SET_MAP), with=FALSE],
                id.vars="CellType", variable.name="Set", value.name="Z")
ct_long[, ModLabel := factor(names(SET_MAP)[match(Set, SET_MAP)], levels=MOD_ORDER)]
ct_long[, CellType := factor(CellType, levels=rev(CT_ORDER))]

if (file.exists(FDR_CSV)) {
  ef <- fread(FDR_CSV); setnames(ef, 1, "CellType")
  fdr_long <- melt(ef[CellType %in% CT_ORDER, c("CellType", SET_MAP), with=FALSE],
                   id.vars="CellType", variable.name="Set", value.name="FDR")
  ct_long <- merge(ct_long, fdr_long, by=c("CellType","Set"), all.x=TRUE)
} else {
  ct_long[, FDR := NA_real_]
  warning("permFDR_specific csv not found — no significance stars")
}
ct_long[, star := fifelse(is.na(FDR), "", fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", "")))]

p_ct <- ggplot(ct_long, aes(x=ModLabel, y=CellType, fill=Z)) +
  geom_tile(color="white",linewidth=0.5) +
  geom_text(aes(label=star),size=3.2,color="grey15",family="Arial") +
  scale_fill_gradient2(low="#3B4F8A",mid="white",high="#D6404C",midpoint=0,
    name="Z",guide=guide_colorbar(barwidth=unit(2,"mm"),barheight=unit(12,"mm"))) +
  labs(x=NULL,y=NULL,title="Cell-Type Enrichment (Day1 validated sets)",
    subtitle="Z = GSE274819 scRNA enrichment; * FDR<0.05 ** <0.01 (1,000 perms)") +
  theme_nature + theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(angle=0,size=6),
    legend.position="right")

# ── 4. Summary table — Day1 metrics (Peak keeps full trajectory, excl. Sham) ──
peak_dt <- frac_all[Group != "NC", .N, by=.(ModLabel,Group)][, .SD[which.max(N)], by=ModLabel]
cds_dt <- frac[, .(CDS_pct = mean(Gene_Struc=="CDS")), by=ModLabel]
ct_top <- ct_long[, .SD[which.max(Z)], by=ModLabel]
summary_dt <- data.table(Mod=MOD_ORDER)
summary_dt[, Sites := site_tot[match(Mod,ModLabel), format(N,big.mark=",")]]
summary_dt[, CDS := percent(cds_dt[match(Mod,ModLabel), CDS_pct], accuracy=1)]
u3sub <- reg_sum[Region=="3'UTR"]
u3idx <- match(MOD_ORDER, u3sub$ModLabel)
u3n <- u3sub$N_3UTR[u3idx]; u3f <- u3sub$UTR3_flag[u3idx]
u3n[is.na(u3n)] <- 0; u3f[is.na(u3f)] <- "*"
summary_dt[, UTR3 := paste0(u3n, u3f)]
summary_dt[, CT := ct_top[match(Mod,ModLabel), as.character(CellType)]]
summary_dt[, Peak := peak_dt[match(Mod,ModLabel), Group]]
summary_dt[, ypos:=rev(seq_len(.N))]
summary_dt[, label:=sprintf("%-7s %6s %4s %5s %-6s %-4s",
  Mod,Sites,CDS,UTR3,CT,Peak)]

p_table <- ggplot(summary_dt,aes(x=0,y=ypos)) +
  geom_text(aes(label=label),hjust=0,size=2.4,family="mono") +
  annotate("text",x=0,y=nrow(summary_dt)+0.8,
    label=sprintf("%-7s %6s %4s %5s %-6s %-4s","Mod","Sites","CDS","3'UTR","CT","Peak"),
    hjust=0,size=2.4,color="grey50",family="mono") +
  xlim(-0.1,5.5) + ylim(0,nrow(summary_dt)+1.2) + theme_void() +
  labs(title="Summary (Day1 vs Sham)",
    subtitle="Sites/CDS/3'UTR = Day1 diff sites; CT = top-enriched cell type\nPeak = max diff-site timepoint (excl. Sham)") +
  theme(plot.title=element_text(size=7.5,face="bold",hjust=0,family="Arial"),
    plot.subtitle=element_text(size=5.5,color="grey40",hjust=0,family="Arial"))

# ── Assemble ──
top_row <- wrap_plots(p_sites,p_region,ncol=2,widths=c(0.42,0.58))
bot_row <- wrap_plots(p_ct,p_table,ncol=2,widths=c(0.52,0.48))

combined <- wrap_plots(top_row,bot_row,ncol=1,heights=c(0.45,0.55)) +
  plot_annotation(
    title="Fig3E: Integrated Summary — 5 RNA Modifications at IRI Day1",
    subtitle=paste0(
      "Day1 differential sites (vs Sham) → Regional distribution → Cell-type enrichment → Summary. ",
      "Differential sites: fraction>0.7, |ES|>=0.1, p<0.05."
    ),
    caption=paste0("n=",format(nrow(frac),big.mark=",")," Day1 differential sites (NC_vs_Day1). ",
      "Nm=Am+Cm+Gm+Um (NM_UP validated union). scRNA enrichment GSE274819, validated Day1 direction sets."),
    theme=theme(plot.title=element_text(size=8.5,face="bold",family="Arial"),
      plot.subtitle=element_text(size=6.5,color="grey30",family="Arial"),
      plot.caption=element_text(size=5.5,color="grey50",hjust=0,family="Arial")))

W_MM<-183; H_MM<-120
ggsave(file.path(OUT_DIR,"Fig3E_Integrated_Summary.pdf"),combined,width=W_MM/25.4,height=H_MM/25.4,device=cairo_pdf)
ggsave(file.path(OUT_DIR,"Fig3E_Integrated_Summary.png"),combined,width=W_MM/25.4,height=H_MM/25.4,dpi=600)
cat(sprintf("Fig3E v5 (Day1-focused): %.0f x %.0f mm\n",W_MM,H_MM))
print(summary_dt[, .(Mod,Sites,CDS,UTR3,CT,Peak)])
