library(hexbin)
library(ggplot2)
library(grid)
library(readxl)
library(plotrix)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

# setwd(input_dir)
# load("COOK_branched_pseudotime.seurat.RData")

load("A549_mapped_seurat_correct_all_timecourse.RData")
LUAD_scores_EMT<-df_scores_EMT
LUAD_mapping<-res_mapping

load("MCF7_mapped_seurat_correct_all_timecourse.RData")
BRCA_scores_EMT<-df_scores_EMT
BRCA_mapping<-res_mapping

list_cellline<-c("A549","MCF7")
list_tcga<-c("LUAD_mapping","BRCA_mapping")

for(i in 1:length(list_cellline)){
  
  setwd(input_dir)
  
  if(list_cellline[i]=="A549"){
    
    clusters_df<-data.frame(read_excel("LUAD_groups_for_MultiTRJ.xlsx"))
  
  }else{
  
    clusters_df<-data.frame(read_excel("BRCA_PRAD_OV_groups_for_MultiTRJ.xlsx","BRCA"))
    
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
  
  setwd(output_dir)
  
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

