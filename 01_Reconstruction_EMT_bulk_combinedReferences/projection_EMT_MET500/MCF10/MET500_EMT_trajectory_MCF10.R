################################################
#MET500 MCF10 reference
################################################


library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

#setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/MCF10/")
setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/MCF10")

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


#setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/Data")


cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

mock <- Biobase::exprs(cds.list[["Mock"]])
tgfb <- Biobase::exprs(cds.list[["TGFB"]])

mock_genes<-cds.list[["Mock"]]@featureData@data[,c(1,2)]
tgfb_genes<-cds.list[["TGFB"]]@featureData@data[,c(1,2)]

mock2<-data.frame(genes=mock_genes,as.matrix(mock))
tgfb2<-data.frame(genes=tgfb_genes,as.matrix(tgfb))

mock_pData <- pData(cds.list[["Mock"]])[,c("cell","sample","Pseudotime")]
tgfb_pData <- pData(cds.list[["TGFB"]])[,c("cell","sample","Pseudotime")]

# 
# Normalization, see the documentaiton of monocle for the type of procedures, and monocle:::normalize_expr_data
# 

mock_norm <- norm_expr_data(mock2[,-c(1:2)],1)
mock_norm2 <- data.frame(mock2[,c(1:2)],mock_norm)

tgfb_norm <- norm_expr_data(tgfb2[,-c(1:2)],1)
tgfb_norm2 <- data.frame(tgfb2[,c(1:2)],tgfb_norm)

print(range(mock_norm2[,-c(1:2)]))
print(range(tgfb_norm2[,-c(1:2)]))

setwd(input_dir)
save(mock_norm2, file = "mock_norm2.RData")
save(tgfb_norm2, file = "tgfb_norm2.RData")
save(mock_pData, file = "mock_pData.RData")
save(tgfb_pData, file = "tgfb_pData.RData")


setwd(input_dir)
load("mock_norm2.RData")
load("mock_pData.RData")
load("tgfb_norm2.RData")
load("tgfb_pData.RData")
mock_pData$Pseudotime<-(mock_pData$Pseudotime/max(mock_pData$Pseudotime))*100
mock_pData$Pseudotime <- 100 - mock_pData$Pseudotime
tgfb_pData$Pseudotime<-(tgfb_pData$Pseudotime/max(tgfb_pData$Pseudotime))*100
tgfb_pData$Pseudotime <- 100 - tgfb_pData$Pseudotime

# 
# upload TCGA
# 
#setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/Data")

load('scaled_expr_data.RData')
gene_symbol <- rownames(scaled.data)
gene_symbol <- data.frame(gene_symbol)
scaled.data <- cbind(gene_symbol, scaled.data)
tcga<-scaled.data
rm(scaled.data)


# UPLOAD MET500 DATA-SET
load('MET500_rna_seq_log2_fpkm.ALLGENES.RData')
met500_unlog<-2^met500_rnaseq2[,-c(1:2)]

fpkmToTpm <- function(fpkm) {
  
  (fpkm / sum(fpkm)  *(1e6))
  
}

met500_convert<-log2(apply(met500_unlog[,-c(1:2)],1,fpkmToTpm)+1)
colnames(met500_convert)<-met500_rnaseq2[,2]
met500_convert<-t(met500_convert)
print(range(met500_convert))

# get the commons gene-names between TCGA and MET500
tcga_genes<-tcga[,1]
met500_genes<-rownames(met500_convert)
common_genes<-intersect(tcga_genes,met500_genes)

tcga<-tcga[which(tcga[,1]%in%common_genes),]

met500_convert2<-met500_convert[which(rownames(met500_convert)%in%common_genes),]
colnames(met500_convert2)<-as.character(sapply(strsplit(colnames(met500_convert2),split='\\.'),'[[',1))
met500_convert2<-data.frame(genes=rownames(met500_convert2),met500_convert2)

met500_convert2<-as.data.frame(setDT(met500_convert2)[, lapply(.SD, mean), by = .(genes)])
tcga<-as.data.frame(setDT(tcga)[, lapply(.SD, mean), by = .(gene_symbol)])
#tcga[,-1]<-log(tcga[,-1]+1,2)

print(range(tcga[,-1]))

input_ge<-merge(tcga,met500_convert2,by.x='gene_symbol',by.y='genes')
rownames(input_ge)<-input_ge[,1]

input_mock<-merge(input_ge,mock_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
#input_mock[input_mock<0]<-0 # trim negative values
input_tgfb<-merge(input_ge,tgfb_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
#input_tgfb[input_tgfb<0]<-0 # trim negative values

# 
# Combat Correction 1) TCGA 2) MET500 3)scRNA-seq
# 
library(sva)

n_tcga<-length(grep(colnames(tcga),pattern='TCGA'))
n_met500<-length(colnames(met500_convert2)[-1])
n_scrnaseq<-ncol(mock_norm2[,-c(1:2)])

batch<-factor(c(rep(1,n_tcga),rep(2,n_met500),rep(3,n_scrnaseq)))

combat_mock = ComBat(dat=as.matrix(input_mock[,-1]), batch=batch)
combat_mock2<-data.frame(genes=input_mock[,1],combat_mock)
combat_mock2[,1]<-as.character(combat_mock2[,1])

n_tcga<-length(grep(colnames(tcga),pattern='TCGA'))
n_met500<-length(colnames(met500_convert2)[-1])
n_scrnaseq<-ncol(tgfb_norm2[,-c(1:2)])

batch<-factor(c(rep(1,n_tcga),rep(2,n_met500),rep(3,n_scrnaseq)))

combat_tgfb = ComBat(dat=as.matrix(input_tgfb[,-1]), batch=batch)
combat_tgfb2<-data.frame(genes=input_tgfb[,1],combat_tgfb)
combat_tgfb2[,1]<-as.character(combat_tgfb2[,1])


# 
# Filter samples with low variance  
# 

vardata<-apply(combat_mock2[,-1],1,var)
idx<-which(vardata>quantile(vardata)[4])
combat_mock_filter<-combat_mock2[idx,]

vardata<-apply(combat_tgfb2[,-1],1,var)
idx<-which(vardata>quantile(vardata)[4])
combat_tgfb2_filter<-combat_tgfb2[idx,]

setwd(output_dir)

library(factoextra)



knnProject<-function(FM,knn_k = 10, num_dim = 25,...){
  
  message("Performing dimensionality reduction with PCA")
  
  irlba_res <- prcomp(t(FM[,-1]))
  
  irlba_pca_res <- irlba_res$x[,c(1:num_dim)]
  rownames(irlba_pca_res)<-colnames(FM[,-1])
  
  tcga_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='agregated',invert=T),]
  print(dim(tcga_pca))
  
  scrna_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='agregated'),]
  print(dim(scrna_pca))
  
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
    
  }
  
  return(knn.list)
  
}

combat_mock_tcga_knn<-knnProject(combat_mock_filter,knn=10,num_dim=30)

combat_tgfb2_tcga_knn<-knnProject(combat_tgfb2_filter,knn=10,num_dim=30)


knn_mock_pseudospace <- NULL

for(patient in names(combat_mock_tcga_knn)){
  
  scRNAseq<-gsub(colnames(combat_mock_tcga_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(mock_pData[mock_pData[,1]%in%scRNAseq,'Pseudotime'])
  knn_mock_pseudospace<-c(knn_mock_pseudospace,knn_pseudospace)
}

tumors<-sapply(strsplit(names(combat_mock_tcga_knn),split='\\.'),'[[',1)
knn_df_tcga_mock<-data.frame(tumors=tumors,patients=names(combat_mock_tcga_knn),pseudospace=knn_mock_pseudospace)

knn_tgfb_pseudospace <- NULL

for(patient in names(combat_tgfb2_tcga_knn)){
  
  scRNAseq<-gsub(colnames(combat_tgfb2_tcga_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(tgfb_pData[tgfb_pData[,1]%in%scRNAseq,'Pseudotime'])
  knn_tgfb_pseudospace<-c(knn_tgfb_pseudospace,knn_pseudospace)
}

tumors<-sapply(strsplit(names(combat_tgfb2_tcga_knn),split='\\.'),'[[',1)
knn_df_tcga_tgfb<-data.frame(tumors=tumors,patients=names(combat_tgfb2_tcga_knn),pseudospace=knn_tgfb_pseudospace)


knn_df_tcga_met500_mock<-knn_df_tcga_mock
colnames(knn_df_tcga_met500_mock)[ncol(knn_df_tcga_met500_mock)]<-"mock_pseudospace"

knn_df_tcga_met500_tgfb<-knn_df_tcga_tgfb
colnames(knn_df_tcga_met500_tgfb)[ncol(knn_df_tcga_met500_tgfb)]<-"tgfb_pseudospace"

input_cp<-cbind(knn_df_tcga_met500_mock,knn_df_tcga_met500_tgfb)

epithelial_cancer<-c("ACC","BLCA","BRCA","CESC","CHOL",
                     "COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC",
                     "LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM",
                     "STAD","THYM","THCA","UCS","UCEC","UVM")

input_cp_epi<-input_cp[which(input_cp$tumors%in%epithelial_cancer),]
input_cp_met500<-input_cp[grep(input_cp$patients,pattern="TCGA",invert=T),]
input_cp_final<-rbind(input_cp_epi,input_cp_met500)
write.table(input_cp_final,file='KNN_projection_TCGA_MET500_to_MCF10A_treated_cells.txt',row.names=F,quote=F,sep='\t')



