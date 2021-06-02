library(data.table)
library('biomaRt')
library(phenopath)
library(openxlsx)
library(ggplot2)
library(FNN)
library(irlba)
library(sva)
library(factoextra)
library(monocle)
library(depmixS4)

setwd("/home/guidantoniomt/pseudospace/gdsc")
source("knnForGDSC_MCF.R")

# 
# Load the expression MCF7 data
#
norm_expr_data<-function(FM,pseudo_expr){
  
  FM2 <- FM/monocle:::estimateSizeFactorsForMatrix(FM)
  if (is.null(pseudo_expr))
    pseudo_expr <- 1
  FM2 <- FM2 + pseudo_expr
  FM2 <- log2(FM2)
  
  return(FM2)
}

setwd('/home/guidantoniomt/pseudospace/input_pseudospace')

cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

mock <- Biobase::exprs(cds.list[["Mock"]])

mock_genes<-cds.list[["Mock"]]@featureData@data[,c(1,2)]

mock2<-data.frame(genes=mock_genes,as.matrix(mock))

mock_pData <- pData(cds.list[["Mock"]])[,c("cell","sample","Pseudotime")]

mock_norm <- norm_expr_data(mock2[,-c(1:2)],1)
mock_norm2 <- data.frame(mock2[,c(1:2)],mock_norm)


# 
# List of metastatic genes
#

setwd("/home/guidantoniomt/pseudospace/HMM")
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
markers_genes<-markers_genes_read[,1]

# 
# Analysis GDSC
#
setwd("/home/guidantoniomt/pseudospace/gdsc")

tab_input<-fread(file="/home/guidantoniomt/pseudospace/gdsc/gdsc_rma_exp.txt",data.table=F)

tab_annotation<-read.xlsx("/home/guidantoniomt/pseudospace/gdsc/Cell_Lines_Details.xlsx",sheet=2)
cell_line_to_use<-tab_annotation[-which(tab_annotation$Site%in%c("autonomic_ganglia","bone","central_nervous_system","haematopoietic_and_lymphoid_tissue","NS")),1]

tab_input2<-tab_input[,colnames(tab_input)%in%c("gene_symbol",cell_line_to_use)]

dim(tab_input2)

#
# Now map the bulk data onto the gdsc data and compute the metastatic potential for each tissue
#

gdsc_for_knn<-tab_input2

combined_MCF_gdsc<-merge(mock_norm2,gdsc_for_knn,by.x="genes.gene_short_name",by.y="gene_symbol",all.y=F,all.x=F)

irlba_res <- prcomp(t(combined_MCF_gdsc[,-c(1:2)]))

#
# PCA pre-combat
#
setwd("/home/guidantoniomt/pseudospace/gdsc/")

pdf('pca_MCF_gdsc_nocombat.pdf')
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-2),rep("gdsc",ncol(gdsc_for_knn)-1)))
p<-fviz_pca_ind(irlba_res,
                col.ind = groups, # color by groups
                palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = FALSE, # Concentration ellipses
                legend.title = "Groups",
                repel=FALSE, 
                geom="point")
print(p)
dev.off()

#
# Combat
#
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-2),rep("gdsc",ncol(gdsc_for_knn)-1)))

dat = combined_MCF_gdsc[,-c(1:2)]

dat <- as.matrix(dat)
batch <- as.factor(groups)

zero.rows.lst <- lapply(levels(batch), function(batch_level) {
  if (sum(batch == batch_level) > 1) {
    return(which(apply(dat[, batch == batch_level], 1,
                       function(x) {
                         var(x) == 0
                       })))
  }else {
    return(which(rep(1, 3) == 2))
  }    
})

zero.rows <- Reduce(union, zero.rows.lst)
keep.rows <- setdiff(1:nrow(dat), zero.rows)
if (length(zero.rows) > 0) {
  dat.orig <- dat
  dat <- dat[keep.rows, ]
}

combined_MCF_gdsc2<-combined_MCF_gdsc[keep.rows,]

combat_MCF_gdsc = ComBat(dat=as.matrix(combined_MCF_gdsc2[,-c(1:2)]), batch=groups)
rownames(combat_MCF_gdsc)<-combined_MCF_gdsc2[,1]

irlba_res_combat <- prcomp(t(combat_MCF_gdsc))

pdf('pca_MCF_gdsc_combat.pdf')
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-2),rep("gdsc",ncol(gdsc_for_knn)-1)))
p<-fviz_pca_ind(irlba_res_combat,
                col.ind = groups, # color by groups
                palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = FALSE, # Concentration ellipses
                legend.title = "Groups",
                repel=FALSE, 
                geom="point")
print(p)
dev.off()

#
# Get only the most variable genes (top 30%)
#
vardata<-apply(combat_MCF_gdsc,1,var)
combat_MCF_gdsc2<-cbind(combat_MCF_gdsc,vardata=vardata)
combat_MCF_gdsc2<-combat_MCF_gdsc2[order(combat_MCF_gdsc2[,"vardata"],decreasing=T),]
combat_MCF_gdsc_filter<- combat_MCF_gdsc2[1:round((nrow(combat_MCF_gdsc2)*30)/100),-ncol(combat_MCF_gdsc2)]

#determine number of components to use

pca_for_components <- prcomp(t(combat_MCF_gdsc_filter))

pdf('MCF_gdsc_PrincipalComponents.pdf',width=15)
fviz_eig(pca_for_components,ncp=50,addlabels=T)
#with 30 PC i explain 33.% of variance 30 is also the point in which the variance become 0.1
variance_explained<-sum(get_eig(pca_for_components)[1:30,2])
dev.off()

FM = pca_for_components
num_dim = 30
knn_k = 10
ann_pseudotime = mock_pData

#
#First output created
#

pseudotimegdscderivedMCF<-knnForGDSC_MCF(FM, knn_k=knn_k,num_dim=num_dim,ann_pseudotime)
colnames(pseudotimegdscderivedMCF)[1]<-"gdsc_ID"

#
# QC: markers and along pseudotime
#
combat_MCF_gdsc2_markers<-gdsc_for_knn[gdsc_for_knn[,1]%in%c("CDH1","CRB3","DSP","CDH2","FN1","VIM"),]
combat_MCF_gdsc2_markers_melt<-melt(combat_MCF_gdsc2_markers)
colnames(combat_MCF_gdsc2_markers_melt)[2]<-"gdsc_ID"
combat_MCF_gdsc2_markers_melt2<-merge(combat_MCF_gdsc2_markers_melt,pseudotimegdscderivedMCF,by.x="gdsc_ID",by.y="gdsc_ID")

combat_MCF_gdsc2_markers_melt2$gene_symbol<-factor(combat_MCF_gdsc2_markers_melt2$gene_symbol,levels=c("CDH1","CRB3","DSP","CDH2","FN1","VIM"))

setwd("/home/guidantoniomt/pseudospace/gdsc/")

pdf("gdsc_EMT_markers_along_MCFtrajectory.pdf")

p<-ggplot(combat_MCF_gdsc2_markers_melt2,
          aes(x=all_pseudotime, y=value))+
  geom_point(shape=16,size=2,color="darkorange3",alpha=0.5)+ylim(0,max(combat_MCF_gdsc2_markers_melt2$value))+scale_x_reverse()+
  facet_wrap(.~gene_symbol,nrow=3,ncol=3)+ 
  geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+ 
  scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p)

dev.off()

#
# Compute the EMT of samples in each gdsc samples 
#

epi_genes<-intersect(gdsc_for_knn[,1], markers_genes_read[markers_genes_read$status=="Epithelial_marker",1])
mes_genes<-intersect(gdsc_for_knn[,1], markers_genes_read[markers_genes_read$status=="Mesenchymal_marker",1])

zscore<-function(x){(x-mean(x))/sd(x)}

gdsc_for_knn_zscore<-data.frame(gene_symbol=gdsc_for_knn[,1],t(apply(gdsc_for_knn[,-1],1,zscore)))

gdsc_for_knn_zscore2<-gdsc_for_knn_zscore[gdsc_for_knn_zscore[,1]%in%c(epi_genes,mes_genes),]
rownames(gdsc_for_knn_zscore2)<-gdsc_for_knn_zscore[gdsc_for_knn_zscore[,1]%in%c(epi_genes,mes_genes),1]

all_scores_emt<-NULL

for(iemt in 2:ncol(gdsc_for_knn_zscore2)){
  
  mean_epi<-mean(gdsc_for_knn_zscore2[epi_genes,iemt])
  mean_mes<-mean(gdsc_for_knn_zscore2[mes_genes,iemt])
  
  score_emt<-mean_mes-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

EMT_gdsc<-data.frame(gdsc_ID=colnames(gdsc_for_knn_zscore2)[-1],emt_score=all_scores_emt)
EMT_gdsc[,1]<-gsub(EMT_gdsc[,1],pattern="\\.",replacement="-")

#
#Second output created
#

gdsc_pseudotime_EMT<-merge(pseudotimegdscderivedMCF,EMT_gdsc,by.x="gdsc_ID")

#
# Plot EMT score gdsc samples along pseudotime
#


pdf("gdsc_EMT_score_along_MCFtrajectory.pdf")

p<-ggplot(gdsc_pseudotime_EMT,
          aes(x=all_pseudotime, y=emt_score))+ geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_point(shape=16,size=2,color="darkorange3",alpha=0.5)+scale_x_reverse()+
  geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+ 
  scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p)

dev.off()

write.table(file="GDSC_on_pseudotime_EMT_MCF_withEMTscores.txt",gdsc_pseudotime_EMT,row.names=F,quote=F,sep="\t")

