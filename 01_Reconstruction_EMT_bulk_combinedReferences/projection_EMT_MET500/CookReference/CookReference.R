################################################
#Cook et at reference
################################################

library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

i <- "OVCA420_TNF"

setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/CookReference")

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

#Load reference data:
setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/Data")
input_string<-paste(i,".rds",sep="")
input <-readRDS(input_string)
no_sti_emt<-rownames(input@meta.data[grep(input@meta.data$Time,pattern="_rm",invert=T),])
matrix_temp<-input@assays$RNA@scale.data
print(dim(matrix_temp))
mock_norm2 <-data.frame(gene_symbol=rownames(input@assays$RNA@scale.data),matrix_temp)
colnames(mock_norm2)[1]<-"genes.gene_short_name"
mock_pData <-data.frame(cell=rownames(input@meta.data),sample=i,Pseudotime=as.numeric(input@meta.data$Pseudotime))
mock_norm2 <- mock_norm2[,colnames(mock_norm2) %in% c("genes.gene_short_name",no_sti_emt)]
mock_pData <- mock_pData[mock_pData$cell %in% no_sti_emt,]
mock_pData$Pseudotime<-(mock_pData$Pseudotime/max(mock_pData$Pseudotime))*100


# 
# upload TCGA
# 
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

input_mock<-merge(input_ge,mock_norm2,by.x='gene_symbol',by.y='genes.gene_short_name')
#input_mock[input_mock<0]<-0 # trim negative values
samples <- colnames(input_ge)[-1]

# 
# Combat Correction 1) TCGA 2) MET500 3)scRNA-seq
# 
library(sva)

n_tcga<-length(grep(colnames(tcga),pattern='TCGA'))
n_met500<-length(colnames(met500_convert2)[-1])
n_scrnaseq<-ncol(mock_norm2[,-c(1)])

batch<-factor(c(rep(1,n_tcga),rep(2,n_met500),rep(3,n_scrnaseq)))

combat_mock = ComBat(dat=as.matrix(input_mock[,-1]), batch=batch)
combat_mock2<-data.frame(genes=input_mock[,1],combat_mock)
combat_mock2[,1]<-as.character(combat_mock2[,1])


# 
# Filter samples with low variance  
# 

vardata<-apply(combat_mock2[,-1],1,var)
idx<-which(vardata>quantile(vardata)[4])
combat_mock_filter<-combat_mock2[idx,]


setwd(output_dir)

library(factoextra)

knnProject<-function(FM,knn_k = 10, num_dim = 25,...){
  
  message("Performing dimensionality reduction with PCA")
  
  irlba_res <- prcomp(t(FM[,-1]))
  
  irlba_pca_res <- irlba_res$x[,c(1:num_dim)]
  rownames(irlba_pca_res)<-colnames(FM[,-1])
  
  tcga_pca<-irlba_pca_res[rownames(irlba_pca_res) %in% samples,]
  scrna_pca<-irlba_pca_res[!(rownames(irlba_pca_res) %in% samples),]
  
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


knn_mock_pseudospace <- NULL

for(patient in names(combat_mock_tcga_knn)){
  
  scRNAseq<-gsub(colnames(combat_mock_tcga_knn[[patient]]$nn.index),pattern='\\.',replacement='-')
  knn_pseudospace<-mean(mock_pData[mock_pData[,1]%in%scRNAseq,'Pseudotime'])
  knn_mock_pseudospace<-c(knn_mock_pseudospace,knn_pseudospace)
}

tumors<-sapply(strsplit(names(combat_mock_tcga_knn),split='\\.'),'[[',1)
knn_df_tcga_mock<-data.frame(tumors=tumors,patients=names(combat_mock_tcga_knn),pseudospace=knn_mock_pseudospace)


knn_df_tcga_met500_mock<-knn_df_tcga_mock
colnames(knn_df_tcga_met500_mock)[ncol(knn_df_tcga_met500_mock)]<-"mock_pseudospace"

epithelial_cancer<-c("ACC","BLCA","BRCA","CESC","CHOL",
                     "COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC",
                     "LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM",
                     "STAD","THYM","THCA","UCS","UCEC","UVM")

input_cp <- knn_df_tcga_met500_mock
input_cp_epi<-input_cp[which(input_cp$tumors%in%epithelial_cancer),]
input_cp_met500<-input_cp[grep(input_cp$patients,pattern="TCGA",invert=T),]
input_cp_final<-rbind(input_cp_epi,input_cp_met500)
write.table(input_cp_final,file=paste("KNN_projection_TCGA_MET500_to_",i,"_treated_cells.txt",sep = ""),row.names=F,quote=F,sep='\t')

