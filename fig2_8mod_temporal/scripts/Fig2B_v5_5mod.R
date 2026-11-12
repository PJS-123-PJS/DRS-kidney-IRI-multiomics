# Fig2B v5 — 5-mod Ridge plot (Nature-optimized, Nm merged)
# Fixes: 2'-O-Me→Nm | Nature theme | unified palette | proper legend
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggridges) })

BASE <- "/mnt/e/BN2505012AH02S66N1/DirectRNA/05.methy"
OUTDIR <- "../outputs"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

BASE_MODS <- c("m6a"="m6A","m5c"="m5C","psu"="PSU","inosine"="Inosine")
NM_MODS <- c("Am","Cm","Gm","Um")
COMPS <- c("NC_vs_Day1","NC_vs_Day3","NC_vs_Day7")
COMP_LABELS <- c("NC→Day1","NC→Day3","NC→Day7")

# Read per-site diff files
all_data <- rbindlist(lapply(c(names(BASE_MODS),NM_MODS), function(mod) {
  rbindlist(lapply(seq_along(COMPS), function(ci) {
    f <- file.path(BASE,mod,"site_diff",COMPS[ci],
                   paste0(COMPS[ci],".",mod,".per_site.diff.add_name.xls"))
    if(!file.exists(f)) return(NULL)
    dt <- fread(f,select=c("effect_size","diffType"))
    dt <- dt[diffType!="No Significant Difference"]
    dt[,`:=`(Mod=mod,Comparison=COMP_LABELS[ci],abs_es=abs(effect_size))]
  }))
}))

# Map to 5-mod labels
all_data[Mod %in% NM_MODS, Mod:="Nm"]
for(bm in names(BASE_MODS)) all_data[Mod==bm, Mod:=BASE_MODS[bm]]
all_data[,Comparison:=factor(Comparison,levels=COMP_LABELS)]
all_data[,Mod:=factor(Mod,levels=c("m6A","m5C","PSU","Inosine","Nm"))]

# Stats
stats_dt <- all_data[,.N,by=.(Mod,Comparison)][order(Comparison,Mod)]
cat("=== Per-mod per-comparison N ===\n"); print(stats_dt)

# Nature palette + theme
MOD_COLORS <- c("m6A"="#D6404C","m5C"="#E8753A","PSU"="#3D7EA6","Inosine"="#7B5294","Nm"="#2E7D6F")
theme_nature <- theme_classic(base_size=7,base_family="Arial")+
  theme(axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
    strip.text=element_text(size=8,face="bold"),strip.background=element_blank(),
    plot.title=element_text(size=9,face="bold",hjust=0.5),
    plot.subtitle=element_text(size=7,color="grey30",hjust=0.5),
    panel.grid.major.x=element_line(color="grey92",linewidth=0.2))

p <- ggplot(all_data,aes(x=abs_es,y=Mod,fill=Mod))+
  geom_density_ridges(scale=1.2,alpha=0.75,linewidth=0.25,
    quantile_lines=TRUE,quantiles=2,color="grey40")+
  facet_wrap(~Comparison,ncol=3)+
  scale_fill_manual(values=MOD_COLORS,guide="none")+
  scale_x_continuous(limits=c(0,0.8),breaks=seq(0,0.8,0.2),expand=c(0,0))+
  labs(x=expression("|"*Delta*"| (per-site fraction)"),y=NULL,
    title="Continuous remodeling spectrum across IRI recovery (5-mod architecture)",
    subtitle="Nm=Am+Cm+Gm+Um merged. Differential sites: |ES|>=0.1, p<0.05.")+
  theme_nature

W_MM<-183; H_MM<-90
ggsave(file.path(OUTDIR,"Fig2B_Ridge_v5_5mod.pdf"),p,width=W_MM/25.4,height=H_MM/25.4,device=cairo_pdf)
ggsave(file.path(OUTDIR,"Fig2B_Ridge_v5_5mod.png"),p,width=W_MM/25.4,height=H_MM/25.4,dpi=600)
cat(sprintf("Fig2B v5: %.0f x %.0f mm\n",W_MM,H_MM))
