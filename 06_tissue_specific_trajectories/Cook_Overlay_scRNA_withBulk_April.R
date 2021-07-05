library(hexbin)
library(ggplot2)
library(grid)
library(readxl)
library(plotrix)

setwd("/data/pseudospace/explore_cook_sort_pseudotimes")
load("COOK_branched_pseudotime.seurat.RData")

load("/data/pseudospace/res_multiple_pseudospace/A549_mapped_seurat_correct_all_timecourse.RData")
LUAD_scores_EMT<-df_scores_EMT
LUAD_mapping<-res_mapping

load("/data/pseudospace/res_multiple_pseudospace/MCF7_mapped_seurat_correct_all_timecourse.RData")
BRCA_scores_EMT<-df_scores_EMT
BRCA_mapping<-res_mapping

list_cellline<-c("A549","MCF7")
list_tcga<-c("LUAD_mapping","BRCA_mapping")

for(i in 1:length(list_cellline)){

  if(list_cellline[i]=="A549"){
    
    clusters_df<-data.frame(read_excel("/data/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events/LUAD_groups_for_Lucie.xlsx"))
  
  }else{
  
    clusters_df<-data.frame(read_excel("/data/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events/data_for_lucie_BRCA_PRAD_OV.xlsx","BRCA"))
    
  }
  
  lc<-list_cellline[i]

  results_celline<-grep(ls(),pattern=paste(lc,"deconstruct_slingshot",sep="_"),value=T)
  current_cellline<-get(results_celline)
  
  lineage_celline<-get(grep(ls(),pattern=paste(lc,"lineage",sep="_"),value=T))
  
  pca_res_cl<-cbind(data.frame(current_cellline[[2]]),"scRNAseq")
  colnames(pca_res_cl)<-c("x","y","exp")

  connections_cl<-current_cellline[[5]]
  
  current_tcga_pca<-data.frame(get(list_tcga[i]),exp="TCGA")
  current_tcga_pca<-merge(current_tcga_pca,clusters_df,by.x="Samples",by.y="Samples")
  tcga_pca<-cbind(current_tcga_pca[,which(colnames(current_tcga_pca)%in%c("x","y","clusters"))])
  colnames(tcga_pca)<-c("x","y","exp")
  
  setwd("/data/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")
  
  # all_data<-rbind(pca_res_cl,tcga_pca)
  # all_data$exp<-as.factor(all_data$exp)
  
  # 
  # # p<-ggplot(all_data, aes(x=x, y=y)) + geom_point(aes(color=exp)) +  scale_color_manual(values=c('#999999','#E69F00'))+theme_bw()
  # 
  # plot(x=all_data[,"x"],y=all_data[,"y"],col=c(rep('#999999',nrow(pca_res_cl)),rep('#E69F00',nrow(tcga_pca))))
  # par(new=TRUE)
  # lines(lineage_celline,col="midnightblue")
  # 
  # 
  # 
  # # plot(x=hbin@xcm,y=hbin@ycm,col=colors_to_use,pch=23,)
  # # 
  # # p<-gplot.hexbin(hbin)
  # # pushHexport(p$plot.vp)
  # # grid.points(tcga_pca[,"x"], tcga_pca[,"y"], gp=gpar(col="red"))
  # # lines(lineage_celline,col="midnightblue",ylim=c(-3,3),xlim=c(-3,3.5))
  # # 
  # # plot(x=tcga_pca[,"x"],y=tcga_pca[,"y"],ylim=c(-3,3),xlim=c(-3,3.5),col='#E69F00')
  # # upViewport()
  pdf(paste("TCGA_project_onto",lc,"pdf",sep="."))
  
  hexmap <- function(xcor,ycor,colval){
    plot(c(min(xcor),max(xcor)),c(min(ycor),max(ycor)), type="n", frame.plot=F, xlab="", ylab="")
    data <- data.frame(xcor,ycor,colval)
    apply(data, 1, function(zone) {
      hexagon(as.numeric(zone[1]),as.numeric(zone[2]),col=as.character(zone[3]), unitcell=0.15,border=as.character(zone[3]))})
    # text(xcor+ 0.5,ycor + 0.5, cex=0.7)
  }
  
  hbin <- hexbin(x=pca_res_cl[,"x"],y=pca_res_cl[,"y"], xbins = 40)
  
  rbPal <- colorRampPalette(c('gainsboro','black'))
  
  colors_to_use <- rbPal(10)[as.numeric(cut(hbin@count,breaks = 10))]
  
  hexmap(xcor=hbin@xcm,ycor=hbin@ycm, colval=colors_to_use)
  par(new=TRUE)
  
  groups<-tcga_pca$exp
  
  if(list_cellline[i]=="A549"){
    
    groups[groups==3]<-"midnightblue"
    groups[groups==4]<-"green2"
    groups[groups==5]<-"violet"
    groups[groups==1]<-"red2"
    groups[groups==2]<-"orange2"
    
  }else{
    
    groups[groups==4]<-"midnightblue"
    groups[groups==2]<-"green2"
    groups[groups==1]<-"violet"
    groups[groups==5]<-"orange2"
    groups[groups==3]<-"red2"
    
  }
  
  plot(x=tcga_pca[,"x"],y=tcga_pca[,"y"],ylim=c(min(pca_res_cl[,"y"]),max(pca_res_cl[,"y"])),xlim=c(min(pca_res_cl[,"x"]),max(pca_res_cl[,"x"])),col=groups,pch=19,xaxt="n",yaxt="n")
  par(new=TRUE)
  lines(lineage_celline,col="red",ylim=c(min(pca_res_cl[,"y"]),max(pca_res_cl[,"y"])),xlim=c(min(pca_res_cl[,"x"]),max(pca_res_cl[,"x"])))
  
  x<-hbin@xcm
  y<-hbin@ycm
  count<-hbin@count
  dfRes<-data.frame(x,y,count)
  
  colors_to_use <- unique(rbPal(10)[as.numeric(cut(hbin@count,breaks = 10))])
  
  plot(x=rep(5,10),y=seq(min(hbin@count),max(hbin@count),by=10)[1:10],col=colors_to_use,pch=19,cex=3)
  
  dev.off()
  
}

