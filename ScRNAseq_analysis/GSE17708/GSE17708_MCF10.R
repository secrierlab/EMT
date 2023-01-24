########################
#POG MCF10
########################

library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)
library(biomaRt)
library(plyr)
library(dplyr)

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/MCF10/")
#setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/MCF10/")

source('functions_plot_pseudospace.R')
folder_analysis<-getwd()

current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

norm_expr_data<-function(FM,pseudo_expr){
  
  FM2 <- FM/estimateSizeFactorsForMatrix(FM)
  if (is.null(pseudo_expr)) 
    pseudo_expr <- 1
  FM2 <- FM2 + pseudo_expr
  FM2 <- log2(FM2)
  
  return(FM2)
}

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
#setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/Data/")


cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

mock <- Biobase::exprs(cds.list[["Mock"]])
tgfb <- Biobase::exprs(cds.list[["TGFB"]])

#This is the information about reported genes
mock_genes<-cds.list[["Mock"]]@featureData@data[,c(1,2)]
tgfb_genes<-cds.list[["TGFB"]]@featureData@data[,c(1,2)]

#This is the expression data
mock2<-data.frame(genes=mock_genes,as.matrix(mock))
tgfb2<-data.frame(genes=tgfb_genes,as.matrix(tgfb))

#This is the peseudotime data
mock_pData <- pData(cds.list[["Mock"]])[,c("cell","sample","Pseudotime")]
tgfb_pData <- pData(cds.list[["TGFB"]])[,c("cell","sample","Pseudotime")]
mock_pData$Pseudotime<-(mock_pData$Pseudotime/max(mock_pData$Pseudotime))*100
mock_pData$Pseudotime <- 100 - mock_pData$Pseudotime
tgfb_pData$Pseudotime<-(tgfb_pData$Pseudotime/max(tgfb_pData$Pseudotime))*100
tgfb_pData$Pseudotime <- 100 - tgfb_pData$Pseudotime

# 
# Normalization, see the documentaiton of monocle for the type of procedures, and monocle:::normalize_expr_data
# 

mock_norm <- norm_expr_data(mock2[,-c(1:2)],1)
mock_norm2 <- data.frame(mock2[,c(1:2)],mock_norm)

tgfb_norm <- norm_expr_data(tgfb2[,-c(1:2)],1)
tgfb_norm2 <- data.frame(tgfb2[,c(1:2)],tgfb_norm)


# 
# Check scRNA-seq data for VIM
#   


mock_norm2_gene<-t(mock_norm2[mock_norm2[,2]=='FN1',-c(1:2)])
mock_norm2_gene<-data.frame(sample=rownames(mock_norm2_gene),mock_norm2_gene)
mock_norm2_gene[,1]<-gsub(mock_norm2_gene[,1],pattern='\\.',replacement='-')
mock_norm2_gene<-merge(mock_norm2_gene,mock_pData,by.x='sample',by.y='cell')[,c(2,3,4)]
colnames(mock_norm2_gene)[1]<-'value'
colnames(mock_norm2_gene)[2]<-'status'

setwd(output_dir)

png('MOCK_FN1_original_scRNA_seq.png', height=6,width=8,units='in',res=300)
ggplot(mock_norm2_gene,aes(x=Pseudotime,y=value,color=status,shape=status))+geom_point(aes(alpha=0.9))+scale_color_manual(values=c('#E69F00','#56B4E9'))+geom_smooth(method='loess',fullrange=T,aes(linetype=status),color='black')
dev.off()

tgfb_norm2_gene<-t(tgfb_norm2[tgfb_norm2[,2]=='FN1',-c(1:2)])
tgfb_norm2_gene<-data.frame(sample=rownames(tgfb_norm2_gene),tgfb_norm2_gene)
tgfb_norm2_gene[,1]<-gsub(tgfb_norm2_gene[,1],pattern='\\.',replacement='-')
tgfb_norm2_gene<-merge(tgfb_norm2_gene,tgfb_pData,by.x='sample',by.y='cell')[,c(2,3,4)]
colnames(tgfb_norm2_gene)[1]<-'value'
colnames(tgfb_norm2_gene)[2]<-'status'

png('TGFB_FN1_original_scRNA_seq.png', height=6,width=8,units='in',res=300)
ggplot(tgfb_norm2_gene,aes(x=Pseudotime,y=value,color=status,shape=status))+geom_point(aes(alpha=0.9))+scale_color_manual(values=c('#E69F00','#56B4E9'))+geom_smooth(method='loess',fullrange=T,aes(linetype=status),color='black')
dev.off()


setwd(input_dir)
save(mock_pData, file = "mock_pData.RData")
save(tgfb_pData, file = "tgfb_pData.RData")
save(mock_norm2, file = "mock_norm2.RData")
save(tgfb_norm2, file = "tgfb_norm2.RData")



load("mock_norm2.RData")
load("mock_pData.RData")
load("tgfb_norm2.RData")
load("tgfb_pData.RData")


#Load POG dataset:

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/")
#setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/")
load('GSE17708_expr_data.RData')
tcga_samples <- colnames(tcga)
tcga_samples <- tcga_samples[-1]


input_mock<-merge(tcga,mock_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
input_tgfb<-merge(tcga,tgfb_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')

# 
# Combat Correction
# 
library(sva)

batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(mock_norm2[,-c(1:2)]))))
dat <- as.matrix(input_mock[,-1])
batch <- as.factor(batch)
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
input_mock <- input_mock[keep.rows,]
batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(mock_norm2[,-c(1:2)]))))
combat_mock = ComBat(dat=as.matrix(input_mock[,-1]), batch=batch)
combat_mock2<-data.frame(genes=input_mock[,1],combat_mock)
combat_mock2[,1]<-as.character(combat_mock2[,1])
#combat_mock2[combat_mock2<0]<-0


batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(tgfb_norm2[,-c(1:2)]))))
dat <- as.matrix(input_tgfb[,-1])
batch <- as.factor(batch)
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
input_tgfb <- input_tgfb[keep.rows,]
batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(tgfb_norm2[,-c(1:2)]))))
combat_tgfb = ComBat(dat=as.matrix(input_tgfb[,-1]), batch=batch)
combat_tgfb2<-data.frame(genes=input_tgfb[,1],combat_tgfb)
combat_tgfb2[,1]<-as.character(combat_tgfb2[,1])
#combat_tgfb2[combat_tgfb2<0]<-0


# 
# Filter samples with low variance  
# 

vardata<-apply(combat_mock2[,-1],1,var)
combat_mock2$vardata<-vardata
combat_mock2<-combat_mock2[order(combat_mock2$vardata,decreasing=T),]
combat_mock_filter<- combat_mock2[1:round((nrow(combat_mock2)*30)/100),-ncol(combat_mock2)]


vardata<-apply(combat_tgfb2[,-1],1,var)
combat_tgfb2$vardata<-vardata
combat_tgfb2<-combat_tgfb2[order(combat_tgfb2$vardata,decreasing=T),]
combat_tgfb2_filter<- combat_tgfb2[1:round((nrow(combat_tgfb2)*30)/100),-ncol(combat_tgfb2)]



library(factoextra)

knnProject<-function(FM,knn_k = 10, num_dim = 25,...){
  
  message("Performing dimensionality reduction with PCA")
  
  irlba_res <- prcomp(t(FM[,-1]))
  
  irlba_pca_res <- irlba_res$x[,c(1:num_dim)]
  rownames(irlba_pca_res)<-colnames(FM[,-1])
  
  tcga_pca<-irlba_pca_res[rownames(irlba_pca_res) %in% tcga_samples,]
  scrna_pca<-irlba_pca_res[!(rownames(irlba_pca_res) %in% tcga_samples),]
  
  knn.list <- list()
  
  for(patient in 1:length(rownames(tcga_pca))){
    print(patient)
    
    patient2<-rownames(tcga_pca)[patient]
    
    tcga.mtx<-tcga_pca[row.names(tcga_pca)==patient2,]
    
    res<-get.knnx(scrna_pca, t(data.frame(tcga.mtx)), k=knn_k, algorithm="kd_tree")
    
    knn.list[[patient2]]<-res
    row.names(knn.list[[patient2]]$nn.index)<-patient2
    colnames(knn.list[[patient2]]$nn.index)<-rownames(scrna_pca)[res$nn.index]
    row.names(knn.list[[patient2]]$nn.dist) <- patient2
    
    leninner<-length(grep(rownames(scrna_pca)[res$nn.index],pattern='inner'))
    lenouter<-length(grep(rownames(scrna_pca)[res$nn.index],pattern='outer'))
    
  }
  return(knn.list)
  
}

setwd(output_dir)

combat_mock_tcga_knn<-knnProject(combat_mock_filter,num_dim=25)

combat_tgfb2_tcga_knn<-knnProject(combat_tgfb2_filter,num_dim=25)


knn_mock_pseudospace <- NULL
ratio_mock_outer_inner<- NULL 
knn_mock_corr<-NULL

for(patient in names(combat_mock_tcga_knn)){
  
  scRNAseq<-gsub(colnames(combat_mock_tcga_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(mock_pData[mock_pData[,1]%in%scRNAseq,'Pseudotime'])
  knn_mock_pseudospace<-c(knn_mock_pseudospace,knn_pseudospace)
  
  tab<-table( mock_pData[mock_pData[,1]%in%scRNAseq,'sample'])
  
  if(length(tab)>1){
    
    ratio<-table(mock_pData[mock_pData[,1]%in%scRNAseq,'sample'])[2]/table( mock_pData[mock_pData[,1]%in%scRNAseq,'sample'])[1]
    
    if(ratio!=1){
      
      ratio_mock_outer_inner<-c(ratio_mock_outer_inner,ifelse(ratio>1,'outer','inner'))
      
    } else{
      
      ratio_mock_outer_inner<-c(ratio_mock_outer_inner,'neutral')
      
    }
    
  } else {
    
    ratio_mock_outer_inner<-c(ratio_mock_outer_inner,sapply(strsplit(names(tab),split='_'),'[[',2))
    
  }
  
  idx_tcga<-grep(colnames(combat_mock_filter),pattern=patient)
  exp_profile_tcga<-combat_mock_filter[,idx_tcga]
  
  scRNAseq<-colnames(combat_mock_tcga_knn[[patient]]$nn.index)
  idx_sc<-grep(colnames(combat_mock_filter),pattern=paste(scRNAseq,collapse='|'))
  exp_profile_sc<-combat_mock_filter[,idx_sc]
  knn_mock_corr<-c(knn_mock_corr,mean(cor(exp_profile_tcga,exp_profile_sc,method='spearman')))
  
}

knn_df_tcga_mock<-data.frame(patients=names(combat_mock_tcga_knn),pseudospace=knn_mock_pseudospace,ratio=ratio_mock_outer_inner,cor_mock=knn_mock_corr)

knn_tgfb_pseudospace <- NULL
ratio_tgfb_outer_inner<- NULL
knn_tgfb_corr<-NULL

for(patient in names(combat_tgfb2_tcga_knn)){
  
  scRNAseq<-gsub(colnames(combat_tgfb2_tcga_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(tgfb_pData[tgfb_pData[,1]%in%scRNAseq,'Pseudotime'])
  knn_tgfb_pseudospace<-c(knn_tgfb_pseudospace,knn_pseudospace)
  
  tab<-table(tgfb_pData[tgfb_pData[,1]%in%scRNAseq,'sample'])
  
  if(length(tab)>1){
    
    ratio<-table(tgfb_pData[tgfb_pData[,1]%in%scRNAseq,'sample'])[2]/table( tgfb_pData[tgfb_pData[,1]%in%scRNAseq,'sample'])[1]
    
    if(ratio!=1){	
      
      ratio_tgfb_outer_inner<-c(ratio_tgfb_outer_inner,ifelse(ratio>1,'outer','inner'))
      
    } else {
      
      ratio_tgfb_outer_inner<-c(ratio_tgfb_outer_inner,'neutral')
      
    }
    
  } else {
    
    ratio_tgfb_outer_inner<-c(ratio_tgfb_outer_inner,sapply(strsplit(names(tab),split='_'),'[[',2))
    
  }
  
  idx_tcga<-grep(colnames(combat_tgfb2_filter),pattern=patient)
  exp_profile_tcga<-combat_tgfb2_filter[,idx_tcga]
  
  scRNAseq<-colnames(combat_tgfb2_tcga_knn[[patient]]$nn.index)
  idx_sc<-grep(colnames(combat_tgfb2_filter),pattern=paste(scRNAseq,collapse='|'))
  exp_profile_sc<-combat_tgfb2_filter[,idx_sc]
  knn_tgfb_corr<-c(knn_tgfb_corr,mean(cor(exp_profile_tcga,exp_profile_sc,method='spearman')))
}

knn_df_tcga_tgfb<-data.frame(patients=names(combat_tgfb2_tcga_knn),pseudospace=knn_tgfb_pseudospace,ratio=ratio_tgfb_outer_inner,cor_tgfb=knn_tgfb_corr)


setwd(output_dir)

write.table(knn_df_tcga_tgfb,file='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells_no_correction_primarytumor_withzeros.txt',row.names=F,quote=F,sep='\t')
write.table(knn_df_tcga_mock,file='KNN_projection_TCGA_to_MCF10A_mock_treated_cells_no_correction_primarytumor_withzeros.txt',row.names=F,quote=F,sep='\t')
input_cp<-data.frame(knn_df_tcga_mock,knn_df_tcga_tgfb[,c(2:4)])
colnames(input_cp)[2:7]<-c('mock','mock_ratio','cor_mock','tgfb','tgfb_ratio','cor_tgfb')
write.table(input_cp,file="proj_pseudospace_withzeros.txt",quote=F,row.names=F,sep="\t")

