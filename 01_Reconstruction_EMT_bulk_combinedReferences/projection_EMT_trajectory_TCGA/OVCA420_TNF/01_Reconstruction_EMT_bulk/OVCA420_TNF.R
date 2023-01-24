library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

i <- "OVCA420_TNF"

setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/",i,"/01_Reconstruction_EMT_bulk/",sep = ""))
#setwd(paste("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/",i,"/01_Reconstruction_EMT_bulk/",sep = "))

source('functions_plot_pseudospace.R')
folder_analysis<-getwd()

current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
#setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/Data/")

#Load data:
input_string<-paste(i,".rds",sep="")
input <-readRDS(input_string)
no_sti_emt<-rownames(input@meta.data[grep(input@meta.data$Time,pattern="_rm",invert=T),])
matrix_temp<-input@assays$RNA@scale.data
print(dim(matrix_temp))
mock_norm2 <-data.frame(gene_symbol=rownames(input@assays$RNA@scale.data),matrix_temp)
colnames(mock_norm2)[1]<-"genes.gene_short_name"
mock_pData <-data.frame(cell=rownames(input@meta.data),sample=i,Pseudotime=as.numeric(input@meta.data$Pseudotime))
mock_norm2 <- mock_norm2[,colnames(mock_norm2) %in% c("gene_symbol",no_sti_emt)]
mock_pData <- mock_pData[mock_pData$cell %in% no_sti_emt,]
mock_pData$Pseudotime<-(mock_pData$Pseudotime/max(mock_pData$Pseudotime))*100



# 
# upload TCGA
# 

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
#setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/Data/")


load('scaled_expr_data.RData')
gene_symbol <- rownames(scaled.data)
gene_symbol <- data.frame(gene_symbol)
scaled.data <- cbind(gene_symbol, scaled.data)
tcga<-scaled.data
rm(scaled.data)
tcga_samples <- colnames(tcga)
tcga_samples <- tcga_samples[-1]



input_mock<-merge(tcga,mock_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')

# 
# Combat Correction
# 
library(sva)

batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(mock_norm2[,-c(1)]))))
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
batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(mock_norm2[,-c(1)]))))
combat_mock = ComBat(dat=as.matrix(input_mock[,-1]), batch=batch)
combat_mock2<-data.frame(genes=input_mock[,1],combat_mock)
combat_mock2[,1]<-as.character(combat_mock2[,1])
combat_mock2[combat_mock2<0]<-0


# 
# Filter samples with low variance  
# 

vardata<-apply(combat_mock2[,-1],1,var)
combat_mock2$vardata<-vardata
combat_mock2<-combat_mock2[order(combat_mock2$vardata,decreasing=T),]
combat_mock_filter<- combat_mock2[1:round((nrow(combat_mock2)*30)/100),-ncol(combat_mock2)]


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

knn_mock_pseudospace <- NULL


for(patient in names(combat_mock_tcga_knn)){
  
  scRNAseq<-gsub(colnames(combat_mock_tcga_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(mock_pData[mock_pData[,1]%in%scRNAseq,'Pseudotime'])
  knn_mock_pseudospace<-c(knn_mock_pseudospace,knn_pseudospace)
  
  
}

knn_df_tcga_mock<-data.frame(patients=names(combat_mock_tcga_knn),pseudospace=knn_mock_pseudospace)

setwd(output_dir)

write.table(knn_df_tcga_mock,file=paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ""),row.names=F,quote=F,sep='\t')
