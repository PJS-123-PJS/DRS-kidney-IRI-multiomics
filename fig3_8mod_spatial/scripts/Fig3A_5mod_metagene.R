# Fig3 Panel A v4 — canonical structure mapping + Sham labels
# ============================================================================
# v4 (2026-07-23): reads frac_5mod_cache_v2.rds (canonical all_*.rds join —
# v1 cache had a GTF-transcript-model mismatch bug collapsing 3'UTR to ~1%);
# NC relabeled Sham per project convention.
# v5 (2026-07-29): REMOVED hardcoded "0 sites in 3'UTR" INOSINE annotation —
# it was a v1-cache-era leftover; v2 data has 106 INOSINE 3'UTR sites and the
# real 3'UTR boxplot is drawn. Audit: Fig3B checkpoint v8.
# v6 (2026-07-29): reads frac_5mod_cache_v3.rds — true pooled Nm (v2's "Nm"
# was silently Um-only) + NC union-deduped across the 3 comparisons.
# v3 fixes: opacity legend | transcript axis labels | color contrast | criteria in title
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

OUT_DIR <- "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_nature <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line=element_line(linewidth=0.3), axis.ticks=element_line(linewidth=0.3),
        axis.ticks.length=unit(1,"pt"), legend.title=element_text(size=6),
        legend.text=element_text(size=5.5), legend.key.size=unit(8,"pt"),
        strip.text=element_text(size=6.5, face="bold"), strip.background=element_blank(),
        plot.title=element_text(size=7.5, face="bold"), plot.subtitle=element_text(size=6,color="grey30"),
        plot.caption=element_text(size=5.5,color="grey50",hjust=0),
        panel.grid=element_blank(), panel.spacing=unit(1.5,"mm"))

# High-contrast palette (colorblind-aware)
MOD_COLORS <- c("M6A"="#D6404C","M5C"="#E8753A","PSU"="#3D7EA6","INOSINE"="#7B5294","Nm"="#2E7D6F")
GROUP_COLORS <- c("Sham"="#6A9FB5","Day1"="#D6404C","Day3"="#E8753A","Day7"="#2E7D6F")
REGION_FILL <- c("5'UTR"="#F2E8D0","CDS"="#B8D8B8","3'UTR"="#6090C0")
MOD_ORDER <- c("M6A","M5C","PSU","INOSINE","Nm")
TIME_ORDER <- c("Sham","Day1","Day3","Day7")

region_bg <- list(
  annotate("rect",xmin=0,xmax=10,ymin=-Inf,ymax=Inf,fill="#F2E8D0",alpha=0.45),
  annotate("rect",xmin=70,xmax=100,ymin=-Inf,ymax=Inf,fill="#6090C0",alpha=0.25),
  geom_vline(xintercept=c(10,70),linetype="dotted",color="grey60",linewidth=0.2),
  annotate("text",x=5,y=Inf,label="5'UTR",vjust=2,size=2.2,color="grey50",family="Arial"),
  annotate("text",x=40,y=Inf,label="CDS",vjust=2,size=2.2,color="grey50",family="Arial"),
  annotate("text",x=85,y=Inf,label="3'UTR",vjust=2,size=2.2,color="grey50",family="Arial")
)

frac <- readRDS(file.path(OUT_DIR,"frac_5mod_cache_v3.rds")); setDT(frac)
frac[, Group:=factor(Group, levels=c("NC","Day1","Day3","Day7"),
                     labels=c("Sham","Day1","Day3","Day7"))]
frac[, ModLabel:=factor(Modification,levels=MOD_ORDER)]
frac[, Group:=factor(Group,levels=TIME_ORDER)]
frac[, norm_pos:=fcase(
  Gene_Struc=="five_prime_UTR", rel_location*10,
  Gene_Struc=="CDS", 10+(rel_location-1)*60,
  Gene_Struc=="three_prime_UTR", 70+(rel_location-2)*30
)]

# ── UPPER: Density with proper legend ──
compute_density <- function(dt, nbins=200) {
  if(nrow(dt)<30) return(NULL)
  dens <- density(dt$norm_pos,from=0,to=100,n=nbins)
  data.table(Position=dens$x, Density=dens$y, N=nrow(dt))
}
dens_list <- frac[, compute_density(.SD), by=.(ModLabel,Group)]
dens_list <- dens_list[!is.na(Density)]
dens_list[, Density_scaled:=Density/max(Density), by=ModLabel]
dens_list[, Tier:=fifelse(N>=200,"N>=200",fifelse(N>=50,"N>=50","N<50"))]
dens_list[, Tier:=factor(Tier,levels=c("N>=200","N>=50","N<50"))]

# Dummy layer for opacity legend
p_dummy <- ggplot(dens_list, aes(x=Position,y=Density_scaled,alpha=Tier,linewidth=Tier)) +
  geom_line(aes(color="dummy")) + scale_color_manual(values=c(dummy="grey50")) +
  scale_alpha_manual(values=c("N>=200"=0.9,"N>=50"=0.5,"N<50"=0.2), name="Sample size") +
  scale_linewidth_manual(values=c("N>=200"=0.5,"N>=50"=0.3,"N<50"=0.18), name="Sample size") +
  guides(color="none")

p_density <- ggplot(dens_list, aes(x=Position,y=Density_scaled,
  color=Group,group=interaction(Group,ModLabel),alpha=Tier,linewidth=Tier)) +
  region_bg + geom_line() +
  facet_wrap(~ModLabel,ncol=5,scales="free_y") +
  scale_color_manual(values=GROUP_COLORS,name=NULL,guide=guide_legend(override.aes=list(linewidth=0.5))) +
  scale_alpha_manual(values=c("N>=200"=0.9,"N>=50"=0.5,"N<50"=0.2),
    name="Sample size", guide=guide_legend(override.aes=list(linewidth=0.5,color="grey50"))) +
  scale_linewidth_manual(values=c("N>=200"=0.45,"N>=50"=0.3,"N<50"=0.18), guide="none") +
  scale_x_continuous(breaks=c(0,25,50,75,100),expand=c(0,0),
    labels=c("0%","25%","50%","75%","100%")) +
  labs(x="Transcript position (% length)", y="Site density (scaled)",
    title="Site Density Along Transcript — 5 Modifications × 4 Timepoints") +
  theme_nature + theme(legend.position="bottom",legend.box="vertical",
    legend.margin=margin(t=-2))

# ── LOWER: Boxplot with proper legend ──
frac_box <- frac[!is.na(Group) & ModLabel %in% MOD_ORDER]
frac_box[, Region:=factor(Gene_Struc,levels=c("five_prime_UTR","CDS","three_prime_UTR"),
  labels=c("5'UTR","CDS","3'UTR"))]
frac_box[, Group:=factor(Group,levels=TIME_ORDER)]
frac_box[, N_reg:=.N, by=.(ModLabel,Group,Region)]
frac_box[, Tier:=fifelse(N_reg>=100,"N>=100",fifelse(N_reg>=20,"N>=20","N<20"))]
frac_box[, Tier:=factor(Tier,levels=c("N>=100","N>=20","N<20"))]

p_box <- ggplot(frac_box, aes(x=Group,y=Fraction,fill=Region,alpha=Tier)) +
  geom_boxplot(outlier.size=0.2,outlier.alpha=0.25,linewidth=0.25,
    position=position_dodge(width=0.8),width=0.6,median.linewidth=0.4) +
  facet_wrap(~ModLabel,ncol=5,scales="free_y") +
  scale_fill_manual(values=REGION_FILL,name=NULL,
    guide=guide_legend(override.aes=list(alpha=1))) +
  scale_alpha_manual(values=c("N>=100"=1,"N>=20"=0.7,"N<20"=0.35),
    name="Sample size per region", guide=guide_legend(override.aes=list(fill="grey60"))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x=NULL,y="Modification fraction",
    title="Per-Site Modification Fraction — Distribution by Region × Timepoint") +
  theme_nature + theme(legend.position="bottom",legend.box="vertical",
    legend.margin=margin(t=-2),axis.text.x=element_text(size=6))

# (v5: stale hardcoded "0 sites in 3'UTR" INOSINE note removed — v2 data shows
#  real INOSINE 3'UTR sites; the boxplot layer above renders them directly.)

# ── Assemble ──
FILTER_CRITERIA <- "Differential sites: fraction>0.7, |effect_size|>=0.1, p<0.05"
combined <- wrap_plots(p_density,p_box,ncol=1,heights=c(0.4,0.6)) +
  plot_annotation(
    title="Fig3A: Spatial Distribution of Differential RNA Modification Sites Along Transcript Body",
    subtitle=paste0("Upper: Site count density. Lower: Per-site modification fraction. ",FILTER_CRITERIA,"."),
    caption=paste0("n=",format(nrow(frac),big.mark=",")," sites. Nm=Am+Cm+Gm+Um (2'-O-Me merged). ",
      "5'UTR[0-10], CDS[10-70], 3'UTR[70-100]. ",
      "PSU/INOSINE: underpowered for regional analysis (n<200/group)."),
    theme=theme(plot.title=element_text(size=8.5,face="bold",family="Arial"),
      plot.subtitle=element_text(size=6.5,color="grey30",family="Arial",margin=margin(t=2,b=4)),
      plot.caption=element_text(size=5.5,color="grey50",hjust=0,family="Arial",margin=margin(t=4))))

W_MM<-183; H_MM<-115
ggsave(file.path(OUT_DIR,"Fig3A_5mod_Metagene.pdf"),combined,width=W_MM/25.4,height=H_MM/25.4,device=cairo_pdf)
ggsave(file.path(OUT_DIR,"Fig3A_5mod_Metagene.png"),combined,width=W_MM/25.4,height=H_MM/25.4,dpi=600)
cat(sprintf("Fig3A v3: %.0f x %.0f mm\n",W_MM,H_MM))
