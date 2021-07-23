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

setwd("/home/data/pseudospace/CCLE")
source("knnForCCLE_MCF.R")
source("petalChartGMT.R")

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

setwd('/home/data/pseudospace/input_pseudospace')

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

setwd("/home/data/pseudospace/HMM")
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
markers_genes<-markers_genes_read[,1]

# 
# Read metmap data
#

brain_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 1)
colnames(brain_metmap500)[1]<-"cell_line"

lung_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 2)
colnames(lung_metmap500)[1]<-"cell_line"

liver_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 3)
colnames(liver_metmap500)[1]<-"cell_line"

bone_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 4)
colnames(bone_metmap500)[1]<-"cell_line"

kidney_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 5)
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

tab_input<-fread(file="/home/data/pseudospace/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929.txt",data.table=F)

#
# Step1: convert ensembl to gene-symbol
#
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

tab_input[,1]<-sapply(strsplit(tab_input[,1],split="\\."),"[[",1)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = tab_input[,1], mart= mart)
tab_input2<-merge(gene_IDs,tab_input,by.x="ensembl_gene_id",by.y="gene_id")

type1<-sapply(strsplit(colnames(tab_input2_emt),split="_"),"[",2)
type2<-sapply(strsplit(colnames(tab_input2_emt),split="_"),"[",3)

cells_to_use<-grep(grep(grep(grep(colnames(tab_input2_emt),pattern="HAEMATOPOIETIC",invert=T,value=T),pattern="CENTRAL_NERVOUS",invert=T,value=T),pattern="FIBROBLAST",invert=T,value=T),pattern="GANGLIA",invert=T,value=T)

#
# Now map the bulk data onto the CCLE data and compute the metastatic potential for each tissue
#

ccle_for_knn<-tab_input2[,which(colnames(tab_input2)%in% c("hgnc_symbol",cells_to_use))]
colnames(ccle_for_knn)[-1]<-paste("CCLE",colnames(ccle_for_knn)[-1],sep="_")

ccle_for_knn[,-1]<-log(ccle_for_knn[,-1]+1,2)

combined_MCF_ccle<-merge(mock_norm2,ccle_for_knn,by.x="genes.gene_short_name",by.y="hgnc_symbol")

irlba_res <- prcomp(t(combined_MCF_ccle[,-c(1:2)]))

#
# PCA pre-combat
#
setwd("/home/data/pseudospace/CCLE/")

pdf('pca_MCF_ccle_nocombat.pdf')
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-2),rep("CCLE",ncol(ccle_for_knn)-1)))
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
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-2),rep("CCLE",ncol(ccle_for_knn)-1)))

dat = combined_MCF_ccle[,-c(1:2)]

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

combat_MCF_ccle = ComBat(dat=as.matrix(combined_MCF_ccle2[,-c(1:2)]), batch=groups)
rownames(combat_MCF_ccle)<-combined_MCF_ccle2[,1]

irlba_res_combat <- prcomp(t(combat_MCF_ccle))

pdf('pca_MCF_ccle_combat.pdf')
groups<-as.factor(c(rep("MCF",ncol(mock_norm2)-2),rep("CCLE",ncol(ccle_for_knn)-1)))
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
vardata<-apply(combat_MCF_ccle,1,var)
combat_MCF_ccle2<-cbind(combat_MCF_ccle,vardata=vardata)
combat_MCF_ccle2<-combat_MCF_ccle2[order(combat_MCF_ccle2[,"vardata"],decreasing=T),]
combat_MCF_ccle_filter<- combat_MCF_ccle2[1:round((nrow(combat_MCF_ccle2)*30)/100),-ncol(combat_MCF_ccle2)]

#determine number of components to use

pca_for_components <- prcomp(t(combat_MCF_ccle_filter))

pdf('MCF_CCLE_PrincipalComponents.pdf',width=15)
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

#
# QC: markers and along pseudotime
#
combat_MCF_ccle2_markers<-ccle_for_knn[ccle_for_knn[,1]%in%c("CDH1","CRB3","DSP","CDH2","FN1","VIM"),]
combat_MCF_ccle2_markers_melt<-melt(combat_MCF_ccle2_markers)
colnames(combat_MCF_ccle2_markers_melt)[2]<-"CCLE_ID"
combat_MCF_ccle2_markers_melt2<-merge(combat_MCF_ccle2_markers_melt,pseudotimeCCLEderivedMCF,by.x="CCLE_ID",by.y="CCLE_ID")

combat_MCF_ccle2_markers_melt2$hgnc_symbol<-factor(combat_MCF_ccle2_markers_melt2$hgnc_symbol,levels=c("CDH1","CRB3","DSP","CDH2","FN1","VIM"))

setwd("/home/data/pseudospace/CCLE/")

pdf("CCLE_EMT_markers_along_MCFtrajectory.pdf")

p<-ggplot(combat_MCF_ccle2_markers_melt2,
          aes(x=all_pseudotime, y=value))+
  geom_point(shape=16,size=2,color="darkorchid3",alpha=0.5)+ylim(0,max(combat_MCF_ccle2_markers_melt2$value))+scale_x_reverse()+
  facet_wrap(.~hgnc_symbol,nrow=3,ncol=3)+ 
  geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+ 
  scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p)

dev.off()

#
# Compute the EMT of samples in each CCLE samples 
#

epi_genes<-intersect(ccle_for_knn[,1], markers_genes_read[markers_genes_read$status=="Epithelial_marker",1])
mes_genes<-intersect(ccle_for_knn[,1], markers_genes_read[markers_genes_read$status=="Mesenchymal_marker",1])

zscore<-function(x){(x-mean(x))/sd(x)}

ccle_for_knn_zscore<-data.frame(gene_symbol=ccle_for_knn[,1],t(apply(ccle_for_knn[,-1],1,zscore)))

ccle_for_knn_zscore2<-ccle_for_knn_zscore[ccle_for_knn_zscore[,1]%in%c(epi_genes,mes_genes),]
rownames(ccle_for_knn_zscore2)<-ccle_for_knn_zscore[ccle_for_knn_zscore[,1]%in%c(epi_genes,mes_genes),1]
  
all_scores_emt<-NULL

for(iemt in 2:ncol(ccle_for_knn_zscore2)){
  
  mean_epi<-mean(ccle_for_knn_zscore2[epi_genes,iemt])
  mean_mes<-mean(ccle_for_knn_zscore2[mes_genes,iemt])
  
  score_emt<-mean_mes-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

EMT_CCLE<-data.frame(CCLE_ID=colnames(ccle_for_knn_zscore2)[-1],emt_score=all_scores_emt)

#
#Second output created
#

CCLE_pseudotime_EMT<-merge(pseudotimeCCLEderivedMCF,EMT_CCLE,by.x="CCLE_ID")

#
# Plot EMT score CCLE samples along pseudotime
#

pdf("CCLE_EMT_markers_along_MCFtrajectory.pdf")

p<-ggplot(combat_MCF_ccle2_markers_melt2,
          aes(x=all_pseudotime, y=value))+
          geom_point(shape=16,size=2,color="darkorchid3",alpha=0.5)+ylim(0,max(combat_MCF_ccle2_markers_melt2$value))+scale_x_reverse()+
          facet_wrap(.~hgnc_symbol,nrow=3,ncol=3)+ 
          geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+ 
          scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p)

dev.off()

pdf("CCLE_EMT_score_along_MCFtrajectory.pdf")

p<-ggplot(CCLE_pseudotime_EMT,
          aes(x=all_pseudotime, y=emt_score))+ geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_point(shape=16,size=2,color="midnightblue",alpha=0.5)+scale_x_reverse()+
  geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+ 
  scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p)

dev.off()

df_all_metastatic_potential<-do.call(rbind,met_potential_list)
df_all_metastatic_potential2<-data.frame(target_tissue=rownames(df_all_metastatic_potential),df_all_metastatic_potential)
df_all_metastatic_potential2[,1]<-strsplit(df_all_metastatic_potential2[,1])
df_all_metastatic_potential2[,1]<-sapply(strsplit(as.character(df_all_metastatic_potential2[,1]),split="\\."),"[[",1)

CCLE_pseudotime_EMT_temp<-CCLE_pseudotime_EMT
CCLE_pseudotime_EMT_temp[,1]<-gsub(CCLE_pseudotime_EMT_temp[,1],pattern="CCLE_",replacement="")

CCLE_pseudotime_EMT_metpot<-merge(CCLE_pseudotime_EMT_temp,df_all_metastatic_potential2,by.x="CCLE_ID",by.y="cell_line")
CCLE_pseudotime_EMT_metpot$status_metpot<-rep("NO",nrow(CCLE_pseudotime_EMT_metpot))

CCLE_pseudotime_EMT_metpot[which(CCLE_pseudotime_EMT_metpot$mean<=-4),"status_metpot"]<-"non_metastatic"
CCLE_pseudotime_EMT_metpot[which(CCLE_pseudotime_EMT_metpot$mean > -4 & CCLE_pseudotime_EMT_metpot$mean < -2),"status_metpot"]<-"weakly_metastatic"
CCLE_pseudotime_EMT_metpot[which(CCLE_pseudotime_EMT_metpot$mean >= -2),"status_metpot"]<-"metastatic"

#
#Third output created
#

CCLE_pseudotime_EMT_metpot$status_metpot<-factor(CCLE_pseudotime_EMT_metpot$status_metpot,levels=c("non_metastatic","weakly_metastatic","metastatic"))

pdf("CCLE_EMT_stratified_for_metpotential.pdf")
my_comparisons <- list( c("non_metastatic", "weakly_metastatic"), c("non_metastatic", "metastatic"), c("weakly_metastatic", "metastatic") )

p1<-ggviolin(CCLE_pseudotime_EMT_metpot, "status_metpot", "emt_score", fill = "status_metpot",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"),xlab="Metastatic Potential")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
print(p1)
p2<-ggviolin(CCLE_pseudotime_EMT_metpot, "status_metpot", "emt_score", fill = "status_metpot",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="Metastatic Potential")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=0, linetype="dashed", color = "red")+facet_wrap(~target_tissue)
print(p2)
dev.off()

#
# Find States of EMT using depmixs4
#
genes_to_use<-markers_genes

ccle_HMM<-ccle_for_knn[ccle_for_knn[,1]%in%genes_to_use,]
rownames(ccle_HMM)<-ccle_HMM[ccle_HMM[,1]%in%genes_to_use,1]

ccle_HMM2<-data.frame(ID=rownames(t(ccle_HMM[,-1])),t(ccle_HMM[,-1]))

input_HMM<-merge(CCLE_pseudotime_EMT,ccle_HMM2,by.x="CCLE_ID",by.y="ID")
input_HMM2<-input_HMM[order(input_HMM$all_pseudotime,decreasing=T),]

list_HMM<- lapply(as.list(paste(colnames(input_HMM2[,4:ncol(input_HMM2)]),"~1")),as.formula)

list_families<-vector(mode="list",length(list_HMM))

for(lf in 1:length(list_families)){
  
  list_families[[lf]]<-gaussian()
}

nstates=3

set.seed(123)

HMM<-depmix(list_HMM,input_HMM2[,-c(1:3)],nstates=nstates,family=list_families)

HMMfit<-fit(HMM, verbose = FALSE) #fit our model to the data set

#use also summary summary(HMMfit, which = "transition") to get the transitions

idx_trans<-which(names(getpars(HMMfit))=="")

transition_matrix<-matrix(getpars(HMMfit)[idx_trans],byrow=T,nrow=nstates)
rownames(transition_matrix)<-paste("from",1:nstates,sep="")
colnames(transition_matrix)<-paste("to",1:nstates,sep="")

HMMpost<-posterior(HMMfit)
output_HMM<-data.frame(samples=input_HMM2[,1],HMM_states=HMMpost)

#
#forth output created
#

CCLE_pseudotime_EMT_HMM<-merge(CCLE_pseudotime_EMT,output_HMM[,c(1:2)],by.x="CCLE_ID",by.y="samples")
colnames(CCLE_pseudotime_EMT_HMM)[ncol(CCLE_pseudotime_EMT_HMM)]<-"hmm_states"
  
my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )

pdf("CCLE_EMT_stratified_for_HMM.pdf")
p1<-ggviolin(CCLE_pseudotime_EMT_HMM, "hmm_states", "emt_score", fill = "hmm_states",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
print(p1)
dev.off()

library(diagram)

pdf(paste("HMM_CCLE_graph_transition_mock_v2",nstates,"pdf",sep="."))
plotmat(round(transition_matrix,3),
        lwd = 1,
        box.lwd = 2,
        cex.txt = 0.7,
        box.size=0.1,
        box.prop=0.5,
        box.type = "circle",
        arr.length=.1,
        arr.width=.1,
        self.cex = .3,
        self.shifty = -.01,
        self.shiftx = .15,
        main = "",
        relsize=0.9,
        box.col=c('#e69f00','#56b4e9','#ee0c0c'))
dev.off()

CCLE_pseudotime_EMT_HMM$hmm_states<-as.factor(CCLE_pseudotime_EMT_HMM$hmm_states)

CCLE_pseudotime_EMT_HMM_temp<-CCLE_pseudotime_EMT_HMM
CCLE_pseudotime_EMT_HMM_temp[,1]<-gsub(CCLE_pseudotime_EMT_HMM_temp[,1],pattern="CCLE_",replacement="")

CCLE_pseudotime_EMT_metpot_HMM<-merge(CCLE_pseudotime_EMT_HMM_temp,df_all_metastatic_potential2,by.x="CCLE_ID",by.y="cell_line")
CCLE_pseudotime_EMT_metpot_HMM$status_metpot<-rep("NO",nrow(CCLE_pseudotime_EMT_metpot_HMM))

CCLE_pseudotime_EMT_metpot_HMM[which(CCLE_pseudotime_EMT_metpot_HMM$mean<=-4),"status_metpot"]<-"non_metastatic"
CCLE_pseudotime_EMT_metpot_HMM[which(CCLE_pseudotime_EMT_metpot_HMM$mean > -4 & CCLE_pseudotime_EMT_metpot_HMM$mean < -2),"status_metpot"]<-"weakly_metastatic"
CCLE_pseudotime_EMT_metpot_HMM[which(CCLE_pseudotime_EMT_metpot_HMM$mean >= -2),"status_metpot"]<-"metastatic"
#
#fifth output created
#

CCLE_pseudotime_EMT_metpot_HMM$status_metpot<-factor(CCLE_pseudotime_EMT_metpot_HMM$status_metpot,levels=c("non_metastatic","weakly_metastatic","metastatic"))

pdf("CCLE_METPOT_stratified_for_HMM.pdf")

my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )

p1<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "hmm_states", "mean", fill = "hmm_states",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=c(-2), linetype="dashed", color = "red")
print(p1+ylab("Metastatic potential"))
p2<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "hmm_states", "mean", fill = "hmm_states",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=c(-2), linetype="dashed", color = "red")+facet_wrap(~target_tissue)
print(p2+ylab("Metastatic potential"))
dev.off()

pdf("CCLE_PENETRANCE_stratified_for_HMM.pdf")

p1<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "hmm_states", "penetrance", fill = "hmm_states",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=c(0), linetype="dashed", color = "red")
print(p1)
p2<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "hmm_states", "penetrance", fill = "hmm_states",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=c(0), linetype="dashed", color = "red")+facet_wrap(~target_tissue)
print(p2)
dev.off()

my_comparisons <- list( c("non_metastatic", "weakly_metastatic"), c("non_metastatic", "metastatic"), c("weakly_metastatic", "metastatic") )

pdf("CCLE_PENETRANCE_stratified_for_METPOT.pdf")

p1<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "status_metpot", "penetrance", fill = "status_metpot",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=c(0), linetype="dashed", color = "red")
print(p1)
p2<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "status_metpot", "penetrance", fill = "status_metpot",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")+ geom_hline(yintercept=c(0), linetype="dashed", color = "red")+facet_wrap(~target_tissue)
print(p2)
dev.off()

save(list=c("pseudotimeCCLEderivedMCF","CCLE_pseudotime_EMT","CCLE_pseudotime_EMT_metpot","CCLE_pseudotime_EMT_HMM","CCLE_pseudotime_EMT_metpot_HMM"),file="CCLE_pseudotime_MCF.RData")

library("vcd")
# plot just a subset of the table

pdf("CCLE_on_MCF10_stats_groups.pdf")
mtx_stat<-table(CCLE_pseudotime_EMT_metpot_HMM$hmm_states,CCLE_pseudotime_EMT_metpot_HMM$status_metpot)
assoc(mtx_stat, shade = TRUE)
chisq <- chisq.test(mtx_stat)
library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE)
contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
corrplot(contrib, is.cor = FALSE)
dev.off()

pdf("CCLE_on_MCF10_Pseudotime_MetPot.pdf")

my_comparisons <- list( c("non_metastatic", "weakly_metastatic"), c("non_metastatic", "metastatic"), c("weakly_metastatic", "metastatic") )

CCLE_pseudotime_EMT_metpot_HMM$pseudospace2<-max(CCLE_pseudotime_EMT_metpot_HMM$all_pseudotime)-CCLE_pseudotime_EMT_metpot_HMM$all_pseudotime

p1<-ggviolin(CCLE_pseudotime_EMT_metpot_HMM, "status_metpot", "pseudospace2", fill = "status_metpot",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"),xlab="HMM states")+ stat_compare_means(comparisons = my_comparisons,method="wilcox")

print(p1)

dev.off()



pdf("CCLE_EMT_METPOT_PENETRANCE_alonpseudotime_v2.pdf")

p1<-ggplot(CCLE_pseudotime_EMT_metpot_HMM,
          aes(x=all_pseudotime, y=emt_score,colour = status_metpot))+ geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_point(shape=16,size=2,alpha=0.5)+scale_x_reverse()+
  geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+
  scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p1)

p2<-ggplot(CCLE_pseudotime_EMT_metpot_HMM,
          aes(x=all_pseudotime, y=mean,colour = status_metpot))+ geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_point(shape=16,size=2,alpha=0.5)+scale_x_reverse()+
  geom_smooth(method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+
  scale_linetype_manual(values=c("solid"))+theme_classic(base_size=15)

print(p2)

dev.off()

petalChartGMT(CCLE_pseudotime_EMT_metpot_HMM,var_group="hmm_states",title_string="Met.Potential HMM states",output="CCLE_petalchart_MetPot_for_HMMstates.pdf")
petalChartGMT_EMT_for_METPOT(CCLE_pseudotime_EMT_metpot_HMM,var_group="status_metpot",title_string="EMT scores Met.Potential",output="CCLE_petalchart_EMTscores_for_MetPot.pdf")

