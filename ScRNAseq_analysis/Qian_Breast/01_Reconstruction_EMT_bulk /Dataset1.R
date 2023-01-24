library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

setwd("~/Documents/EMT_paper_revision/01_Reconstruction_EMT_bulk")
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

setwd(input_dir)

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

# 
# upload signatures
# 

# 
# upload TCGA
# 
setwd("~/Documents/EMT_paper_revision/Data")
load('TPM_data.RData')

tcga<-TPM.data
tcga[,-1]<-log(tcga[,-1]+1,2) #Data has to be log normalised
colnames(tcga)[1] <- "gene_symbol"

#Merge the datasets
input_mock<-merge(tcga,mock_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
input_tgfb<-merge(tcga,tgfb_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')

# 
# Combat Correction
# 
library(sva)

batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(mock_norm2[,-c(1:2)]))))

rm(cds.list)
rm(TPM.data)

combat_mock = ComBat(dat=as.matrix(input_mock[,-1]), batch=batch)
combat_mock2<-data.frame(genes=input_mock[,1],combat_mock)
combat_mock2[,1]<-as.character(combat_mock2[,1])
combat_mock2[combat_mock2<0]<-0
rm(input_mock)

batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(tgfb_norm2[,-c(1:2)]))))

combat_tgfb = ComBat(dat=as.matrix(input_tgfb[,-1]), batch=batch)
combat_tgfb2<-data.frame(genes=input_tgfb[,1],combat_tgfb)
combat_tgfb2[,1]<-as.character(combat_tgfb2[,1])
combat_tgfb2[combat_tgfb2<0]<-0
rm(input_tgfb)

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


#Save combat treated data:
save(combat_tgfb2_filter, file = "combat_tgfb2_filter.RData")
save(combat_mock_filter, file = "combat_mock_filter.RData")
tcga_samples <- colnames(tcga)
tcga_samples <- tcga_samples[-1]
save(tcga_samples, file = "tcga_samples.RData")
save(mock_pData, file = "mock_pData.RData")
save(tgfb_pData, file = "tgfb_pData.RData")

load("combat_tgfb2_filter.RData")
load("combat_mock_filter.RData")
load("tcga_samples.RData")
load("mock_pData.RData")
load("tgfb_pData.RData")

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

write.table(knn_df_tcga_tgfb,file='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells_no_correction_primarytumor.txt',row.names=F,quote=F,sep='\t')
write.table(knn_df_tcga_mock,file='KNN_projection_TCGA_to_MCF10A_mock_treated_cells_no_correction_primarytumor.txt',row.names=F,quote=F,sep='\t')

knn_df_tcga_tgfb<-read.delim(file="KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells_no_correction_primarytumor.txt")
knn_df_tcga_mock<-read.delim(file="KNN_projection_TCGA_to_MCF10A_mock_treated_cells_no_correction_primarytumor.txt")

input_cp<-data.frame(knn_df_tcga_mock,knn_df_tcga_tgfb[,c(2:4)])
colnames(input_cp)[2:7]<-c('mock','mock_ratio','cor_mock','tgfb','tgfb_ratio','cor_tgfb')
pseudospace_input_sort<-pseudospace_input[order(pseudospace_input$mock,decreasing=F),]

input_cp$mock<-(input_cp$mock/max(input_cp$mock))*100
input_cp$tgfb<-(input_cp$tgfb/max(input_cp$tgfb))*100
write.table(input_cp,file="proj_pseudospace.txt",quote=F,row.names=F,sep="\t")






########################
#Extra analysis
setwd("~/Documents/EMT_paper_revision/GSE75688/01_Reconstruction_EMT_bulk/output_dir/")
pseudospace_input<-read.delim(file="proj_pseudospace.txt")
pseudospace_input$Patient <- sapply(pseudospace_input$patients, function(x)
  strsplit(x,"_")[[1]][1])
pseudospace_input$SampleType <- sapply(pseudospace_input$patients, function(x)
  strsplit(x,"_")[[1]][2])
table(pseudospace_input$Patient)
table(pseudospace_input$SampleType)
pseudospace_input$Stage <- sapply(pseudospace_input$Patient, function(x)
  ifelse(x %in% c("BC03LN","BC07LN"),"Metastases","Primary"))
pseudospace_input$patients2 <- pseudospace_input$patients
pseudospace_input_sort<-pseudospace_input[order(pseudospace_input$mock,decreasing=F),]

#Calculate EMT scores:
setwd("~/Documents/EMT_paper_revision/GSE75688/02_HMM_macrostases_EMT")
epi_genes <- c("CDH1","DSP","OCLN")
mes_genes <- c("VIM","CDH2","FOXC2","SNAI1","SNAI2","TWIST1","FN1","ITGB6","MMP2","MMP3","MMP9","SOX10","GCS")
setwd("~/Documents/EMT_paper_revision/GSE75688/02_HMM_macrostases_EMT/data")
markers <- read.table("EMT_and_pEMT_markers.txt", header = TRUE)
markers <- markers$genes
markers <- c(markers, epi_genes, mes_genes)
markers_genes <- unique(markers)
setwd("~/Documents/EMT_paper_revision/GSE75688/Data")
load('TPM_data.RData')
TCGA_GEXP_ALL <- TPM.data
TCGA_GEXP_ALL_rid<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,]
pseudospace_input_sort<-pseudospace_input_sort[which(pseudospace_input_sort$patients2%in%colnames(TCGA_GEXP_ALL_rid)),]
TCGA_GEXP_ALL_sort<-log(t(TCGA_GEXP_ALL_rid[,match(pseudospace_input_sort$patients2,colnames(TCGA_GEXP_ALL_rid))])+1,2)
colnames(TCGA_GEXP_ALL_sort)<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,1]

TCGA_GEXP_ALL_sort2<-TCGA_GEXP_ALL_sort
zscore<-function(x){(x-mean(x))/sd(x)}
TCGA_GEXP_ALL_sort2<-t(apply(TCGA_GEXP_ALL_sort2,2,zscore))
all_scores_emt<-NULL
mes_genes <- mes_genes[mes_genes %in% rownames(TCGA_GEXP_ALL_sort2)]
epi_genes <- epi_genes[epi_genes %in% rownames(TCGA_GEXP_ALL_sort2)]

for(iemt in 1:ncol(TCGA_GEXP_ALL_sort2)){
  
  mean_epi<-mean(TCGA_GEXP_ALL_sort2[epi_genes,iemt])
  mes_epi<-mean(TCGA_GEXP_ALL_sort2[mes_genes,iemt])
  
  score_emt<-mes_epi-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

pseudospace_input_sort$EMT_score <- all_scores_emt
setwd("~/Documents/EMT_paper_revision/GSE75688/02_HMM_macrostases_EMT")
save(pseudospace_input_sort, file = "pseudotime_input_sort_EMT_score.RData")



#Compare between metastatic and non-metastatic patients
setwd("~/Documents/EMT_paper_revision/GSE75688/01_Reconstruction_EMT_bulk/Figures")
table(pseudospace_input_sort$Stage)
pdf("Dataset1_metastatic_vs_primary_pseudotime.pdf",height = 5,width = 5)
my_comparisons <- list( c("Metastases", "Primary"))
p <- ggplot(pseudospace_input_sort, aes(x=Stage, y=mock, fill = Stage)) + 
  geom_boxplot()+
  labs(x="Cancer Stage", y = "Pseudotime")+
  theme_classic()
p +  stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
dev.off()

pdf("Dataset1_metastatic_vs_primary_EMT.pdf",height = 5,width = 5)
my_comparisons <- list( c("Metastases", "Primary"))
p <- ggplot(pseudospace_input_sort, aes(x=Stage, y=EMT_score, fill = Stage)) + 
  geom_boxplot()+
  labs(x="Cancer Stage", y = "EMT score")+
  theme_classic()
p +  stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
dev.off()



#sc vs bulk
table(pseudospace_input_sort$SampleType)
bulk_patients <- pseudospace_input_sort[pseudospace_input_sort$SampleType %in% "Pooled",]
bulk_patients <- unique(as.character(bulk_patients$Patient))
mean_sc <- NULL
bulk <- NULL 
for (i in bulk_patients) {
  
  print(i)
  test <- pseudospace_input_sort[pseudospace_input_sort$Patient %in% i,]
  test1 <- test[!(test$SampleType %in% "Pooled"),]
  test1 <- mean(test1$mock)
  test2 <- test[test$SampleType %in% "Pooled",]
  test2 <- mean(test2$mock)
  mean_sc <- c(mean_sc, test1)
  bulk <- c(bulk, test2)
  
}
summary <- data.frame(mean_sc, bulk)
pdf("Dataset1_bulk_vs_scRNAseq_pseudotime.pdf", width = 5, height = 5)
p <- ggscatter(summary, x = "bulk", y = "mean_sc",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               cor.coef = TRUE, alpha = 0.1, size = 2)
p + labs(x = "Bulk sample peseudotime", y = "mean scRNAseq sample pseudotime") 
dev.off()

mean_sc <- NULL
bulk <- NULL 
for (i in bulk_patients) {
  
  print(i)
  test <- pseudospace_input_sort[pseudospace_input_sort$Patient %in% i,]
  test1 <- test[!(test$SampleType %in% "Pooled"),]
  test1 <- mean(test1$EMT_score)
  test2 <- test[test$SampleType %in% "Pooled",]
  test2 <- mean(test2$EMT_score)
  mean_sc <- c(mean_sc, test1)
  bulk <- c(bulk, test2)
  
}
summary <- data.frame(mean_sc, bulk)
pdf("Dataset1_bulk_vs_scRNAseq_emt.pdf", width = 5, height = 5)
p <- ggscatter(summary, x = "bulk", y = "mean_sc",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               cor.coef = TRUE, alpha = 0.1, size = 2)
p + labs(x = "Bulk sample emt score", y = "mean scRNAseq sample emt score") 
dev.off()


#Sort into EMT categories:
pseudospace_input_sort$EMT_category <- sapply(pseudospace_input_sort$EMT_score, function(x)
  ifelse(x >1, "High",
         ifelse(x < 0, "Low","Mid")))
pseudospace_input_sort$EMT_category <- factor(pseudospace_input_sort$EMT_category, levels = c("Low","Mid","High"))
table(pseudospace_input_sort$EMT_category)
#Bocplot comparison of how pseudotime differs across groups:
pdf("EMT_category_vs_pseudotime.pdf",height = 5,width = 5)
my_comparisons <- list( c("High", "Low"),c("High","Mid"),c("Mid","Low"))
p <- ggplot(pseudospace_input_sort, aes(x=EMT_category, y=mock, fill = EMT_category)) + 
  geom_boxplot()+
  labs(x="EMT class", y = "Pseudotime")+
  theme_classic()
p +  stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
dev.off()

#Barplot composition:
table(pseudospace_input_sort$SampleType)
bulk_patients <- pseudospace_input_sort[pseudospace_input_sort$SampleType %in% "Pooled",]
bulk_patients <- unique(as.character(bulk_patients$Patient))
composition <- NULL
for (i in bulk_patients) {
  
  print(i)
  test <- pseudospace_input_sort[pseudospace_input_sort$Patient %in% i,]
  test <- test[!(test$SampleType %in% "Pooled"),]
  test.low <- test[test$EMT_category %in% "Low",]
  test.low <- dim(test.low)[1]
  test.mid <- test[test$EMT_category %in% "Mid",]
  test.mid <- dim(test.mid)[1]
  test.high <- test[test$EMT_category %in% "High",]
  test.high <- dim(test.high)[1]
  composition <- c(composition, test.low, test.mid, test.high)
  
}
Sample <- c(rep(bulk_patients[1],3),rep(bulk_patients[2],3),rep(bulk_patients[3],3),rep(bulk_patients[4],3),rep(bulk_patients[5],3),rep(bulk_patients[6],3),rep(bulk_patients[7],3),rep(bulk_patients[8],3),rep(bulk_patients[9],3))
N <- composition
composition <- rep(c("Low","Mid","High"),9)
summary <- data.frame(Sample, composition, N)
summary <- summary[order(summary$Sample),]
summary$composition <- factor(summary$composition, levels = c("High","Mid","Low"))
pdf("Barplot_composition.pdf", height = 5, width = 4)
p <- ggplot(summary, aes(fill=composition, y=N, x=Sample, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45)
dev.off()





#Pseudotime vs emt score
pdf("Dataset1_pseudotime_vs_EMT.pdf", width = 5, height = 5)
p <- ggscatter(pseudospace_input_sort, x = "mock", y = "EMT_score",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               cor.coef = TRUE, alpha = 0.1, size = 2)
p + labs(x = "Pseudotime", y = "EMT score") 
dev.off()

#Heatmap:
library(pheatmap)
library(RColorBrewer)
anno <- pseudospace_input_sort$mock
anno <- data.frame(anno)
rownames(anno) <- pseudospace_input_sort$patients
## Plot shap overall metrics
colfunc <- colorRampPalette(c("#1B7837","grey95","#762A83"))
colfunc(326)
myCols <- colfunc(326)


library(pheatmap)
p<-pheatmap(TCGA_GEXP_ALL_sort2)
pdf("Heatmap.pdf", height = 3,width = 6)
p
dev.off()



#PLOT
pseudotime <- pseudospace_input_sort$mock
anno <- data.frame(pseudotime)
rownames(anno) <- pseudospace_input_sort$patients

## Plot shap overall metrics
colfunc <- colorRampPalette(c("#1B7837","grey95","#762A83"))
colfunc(326)
myCols <- colfunc(326)


library(pheatmap)
p<-pheatmap(mat = t(TCGA_GEXP_ALL_sort), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, annotation_col = anno)
pdf("Heatmap.pdf", height = 5,width = 6)
p
dev.off()








