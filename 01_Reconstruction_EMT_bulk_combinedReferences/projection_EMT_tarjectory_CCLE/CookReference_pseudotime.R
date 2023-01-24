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

i <- "OVCA420_TNF"

setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/")

source("knnForCCLE_MCF.R")
source("petalChartGMT.R")

setwd("/home/annawiecek/EMT/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/CookReference")
folder_analysis<-getwd()
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

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



setwd(input_dir)

# 
# List of metastatic genes
#
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
markers_genes<-markers_genes_read[,1]

# 
# Read metmap data
#


setwd(input_dir)
brain_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 1)
colnames(brain_metmap500)[1]<-"cell_line"

lung_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 2)
colnames(lung_metmap500)[1]<-"cell_line"

liver_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 3)
colnames(liver_metmap500)[1]<-"cell_line"

bone_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 4)
colnames(bone_metmap500)[1]<-"cell_line"

kidney_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 5)
colnames(kidney_metmap500)[1]<-"cell_line"

met_potential_list<-vector(mode="list",length=5)
met_potential_list[[1]]<-brain_metmap500
met_potential_list[[2]]<-lung_metmap500
met_potential_list[[3]]<-liver_metmap500
met_potential_list[[4]]<-bone_metmap500
met_potential_list[[5]]<-kidney_metmap500
names(met_potential_list)<-c("brain","lung","liver","bone","kidney")

# 
# Analysis ccle
#

tab_input<-fread(file="CCLE_RNAseq_rsem_genes_tpm_20180929.txt",data.table=F)

#
# Step1: convert ensembl to gene-symbol
#
tab_input[,1]<-sapply(strsplit(tab_input[,1],split="\\."),"[[",1)
load("gene_IDs.RData")
tab_input2<-merge(gene_IDs,tab_input,by.x="ensembl_gene_id",by.y="gene_id")

type1<-sapply(strsplit(colnames(tab_input2),split="_"),"[",2)
type2<-sapply(strsplit(colnames(tab_input2),split="_"),"[",3)

cells_to_use<-grep(grep(grep(grep(colnames(tab_input2),pattern="HAEMATOPOIETIC",invert=T,value=T),pattern="CENTRAL_NERVOUS",invert=T,value=T),pattern="FIBROBLAST",invert=T,value=T),pattern="GANGLIA",invert=T,value=T)

#
# Now map the bulk data onto the CCLE data and compute the metastatic potential for each tissue
#
setwd(output_dir)

ccle_for_knn<-tab_input2[,which(colnames(tab_input2)%in% c("hgnc_symbol",cells_to_use))]
ccle_for_knn$ensembl_gene_id <- NULL
ccle_for_knn$transcript_ids <- NULL
colnames(ccle_for_knn)[-1]<-paste("CCLE",colnames(ccle_for_knn)[-1],sep="_")

ccle_for_knn[,-1]<-log(ccle_for_knn[,-1]+1,2)

combined_MCF_ccle<-merge(mock_norm2,ccle_for_knn,by.x="genes.gene_short_name",by.y="hgnc_symbol")

irlba_res <- prcomp(t(combined_MCF_ccle[,-c(1)]))


#
# Combat
#
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-1),rep("CCLE",ncol(ccle_for_knn)-1)))

dat = combined_MCF_ccle[,-c(1)]

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

combined_MCF_ccle2<-combined_MCF_ccle[keep.rows,]

combat_MCF_ccle = ComBat(dat=as.matrix(combined_MCF_ccle2[,-c(1)]), batch=groups)
rownames(combat_MCF_ccle)<-combined_MCF_ccle2[,1]

irlba_res_combat <- prcomp(t(combat_MCF_ccle))

# Get only the most variable genes (top 30%)
#
vardata<-apply(combat_MCF_ccle,1,var)
combat_MCF_ccle2<-cbind(combat_MCF_ccle,vardata=vardata)
combat_MCF_ccle2<-combat_MCF_ccle2[order(combat_MCF_ccle2[,"vardata"],decreasing=T),]
combat_MCF_ccle_filter<- combat_MCF_ccle2[1:round((nrow(combat_MCF_ccle2)*30)/100),-ncol(combat_MCF_ccle2)]

#determine number of components to use

pca_for_components <- prcomp(t(combat_MCF_ccle_filter))

pdf(paste("MCF_CCLE_PrincipalComponents_",i,".pdf",sep = ""),width=15)
fviz_eig(pca_for_components,ncp=50,addlabels=T)
#with 30 PC i explain 33.% of variance 30 is also the point in which the variance become 0.1
variance_explained<-sum(get_eig(pca_for_components)[1:30,2])
dev.off()

FM = pca_for_components
num_dim = 30
knn_k = 10
ann_pseudotime = mock_pData
list_ann = met_potential_list

#
#First output created
#

pseudotimeCCLEderivedMCF<-knnForCCLE_MCF(FM, knn_k=knn_k,num_dim=num_dim,ann_pseudotime,list_ann)
colnames(pseudotimeCCLEderivedMCF)[1]<-"CCLE_ID"
save(pseudotimeCCLEderivedMCF, file = paste("Pseudotime_",i,".RData",sep = ""))





