# Fig7 Composite v8.7 — 7-Panel PT-S3 Multi-Omics
# ============================================================================
# A(Circos) B(Layers) C(PT-S3) D(HubNet) E(CellComp) F(Metabolome) G(DRS Adv)
library(png); library(grid)
source("Fig6_config.R")

sec_col <- c("Landscape"="#D6404C","PT-S3"="#E8753A","Mechanism"="#3D7EA6","Advantage"="#2E7D6F")

panels <- list(
  list(id="a",file="Fig7A_4Omics_Circos_Journal",  section="Landscape", row=1,col=1,cspan=2),
  list(id="b",file="Fig7B_LayerOverlap",            section="Landscape", row=2,col=1,cspan=2),
  list(id="c",file="Fig7D_PTS3_MultiLayer",         section="PT-S3",     row=3,col=1,cspan=2),
  list(id="d",file="Fig7F_HubGene_Network",         section="PT-S3",     row=4,col=1,cspan=1),
  list(id="e",file="Fig7I_CellComposition_LayerLink",section="Mechanism", row=4,col=2,cspan=1),
  list(id="f",file="Fig7G_Metabolome_Dynamic",       section="Mechanism", row=5,col=1,cspan=1),
  list(id="g",file="Fig7J_DRS_Advantage",            section="Advantage", row=5,col=2,cspan=1)
)

total_w<-7.1; total_h<-8.5; lm<-0.55; rm<-0.15; tm<-0.25; bm<-0.15
gap_x<-0.08; gap_y<-0.06; ncols<-2
col_w<-(total_w-lm-rm-gap_x)/ncols
row_h<-c(1.6, 1.4, 1.4, 1.4, 1.4)  # 5 rows

compute_rect<-function(row,col,cspan=1){
  y_above<-sum(row_h[1:(row-1)])+(row-1)*gap_y
  y_top<-tm+y_above
  x_left<-lm+(col-1)*(col_w+gap_x)
  w<-col_w*cspan+gap_x*(cspan-1)
  c(x=x_left,y=total_h-y_top,w=w,h=row_h[row])
}

sections<-list(
  list(label="Multi-Omics\nLandscape",rows=1:2,sec="Landscape"),
  list(label="PT-S3\nConvergence",rows=3:4,sec="PT-S3"),
  list(label="Mechanism +\nDRS Advantage",rows=5:5,sec="Mechanism"))

cat("Loading...\n")
panel_grobs<-list()
for(p in panels){
  fpath<-file.path(OUTDIR,paste0(p$file,".png"))
  if(!file.exists(fpath)) fpath<-file.path("..",fpath)
  if(!file.exists(fpath)){cat(sprintf("  MISSING: %s\n",p$file));next}
  img<-readPNG(fpath)
  panel_grobs[[p$id]]<-img
  cat(sprintf("  %s: %dx%d\n",p$id,dim(img)[2],dim(img)[1]))
}

draw<-function(){
  grid.newpage(); grid.rect(gp=gpar(fill="white",col=NA))
  for(s in sections){
    rs<-compute_rect(s$rows[1],1,2); re<-compute_rect(s$rows[length(s$rows)],1,2)
    ym<-(rs["y"]+re["y"]+re["h"])/2
    pushViewport(viewport(x=unit(0.22,"inches"),y=unit(ym,"inches"),
      width=unit(0.44,"inches"),height=unit(1,"npc"),just="center"))
    grid.text(s$label,rot=90,gp=gpar(fontsize=5.5,fontface="bold",col=sec_col[s$sec],fontfamily="Arial"))
    popViewport()
  }
  for(p in panels){
    if(!p$id %in% names(panel_grobs)) next
    rect<-compute_rect(p$row,p$col,p$cspan); img<-panel_grobs[[p$id]]
    pushViewport(viewport(x=unit(rect["x"],"inches"),y=unit(rect["y"],"inches"),
      width=unit(rect["w"],"inches"),height=unit(rect["h"],"inches"),just=c(0,0)))
    grid.rect(gp=gpar(fill=NA,col="grey85",lwd=0.3))
    grid.raster(img,width=unit(1,"npc"),height=unit(1,"npc"),interpolate=TRUE)
    grid.text(p$id,x=unit(0.02,"npc"),y=unit(0.98,"npc"),just=c(0,1),
      gp=gpar(fontsize=8,fontface="bold",col=sec_col[p$section],fontfamily="Arial"))
    popViewport()
  }
  grid.text("PT-S3 Multi-Omics Convergence: DRS Captures 8 Regulatory Layers in a Single Sequencing Run",
    x=unit(total_w/2,"inches"),y=unit(0.06,"inches"),just=c(0.5,0),
    gp=gpar(fontsize=8,fontface="bold",col="grey20",fontfamily="Arial"))
}

cairo_pdf(file.path(OUTDIR,"Fig7_Composite_7panel.pdf"),w=total_w,h=total_h,family="Arial")
draw(); dev.off(); cat("PDF saved\n")
png(file.path(OUTDIR,"Fig7_Composite_7panel.png"),w=total_w,h=total_h,u="in",res=600)
draw(); dev.off(); cat("PNG 600dpi saved\n")
if(requireNamespace("ragg",quietly=TRUE)) {
  ragg::agg_tiff(file.path(OUTDIR,"Fig7_Composite_7panel.tiff"),w=total_w,h=total_h,u="in",res=600,compression="lzw")
  draw(); dev.off(); cat("TIFF LZW saved\n")
}
cat(sprintf("Fig7 7-panel: %.1f x %.1f in\n",total_w,total_h))
