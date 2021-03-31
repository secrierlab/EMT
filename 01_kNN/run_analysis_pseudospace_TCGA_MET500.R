library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)


# 
# setwd('/mnt/data/lab/datasets/GSE114687')

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html

# sce<-read10xCounts("/mnt/data/lab/datasets/GSE114687")

# single_cell_data<-as.matrix(counts(sce))
# colnames(single_cell_data)<-sce@colData@listData[[2]]

# ens<-(sce@rowRanges@elementMetadata@listData[[1]])
# gene_symbol<-sce@rowRanges@elementMetadata@listData[[2]]

# SC2<-data.frame(ens=ens,gene_symbol=gene_symbol,single_cell_data)
# colnames(SC2)<-gsub(colnames(SC2),pattern='\\.',replacement='-')

# setwd('/mnt/data/lab/datasets/GSE114687')
# metadata<-read.delim(file='metadata.tsv',header=F,stringsAsFactors=F)
# metadata$new_column<-sapply(strsplit(metadata[,1],split='_'),'[[',1)

# samples<-sapply(strsplit(metadata[,1],split='_'),'[[',1)

# treatment<-sapply(strsplit(metadata[,1],split='_'),'[[',3)
# description<-metadata[,3]

# new_column<-paste(treatment,description,sep='-')

# SC3<-data.frame(SC2[,c(1:2)],SC2[,colnames(SC2)%in%samples])
# colnames(SC3)<-gsub(colnames(SC3),pattern='\\.',replacement='-')


# new_colnames<-metadata[match(colnames(SC3)[-c(1:2)],metadata$new_column),3]

# colnames(SC3)[3:ncol(SC3)]<-new_colnames

# mock<-data.frame(SC3[,c(1:2)],SC3[,grep(colnames(SC3),pattern='Mock')])
# tgfb<-data.frame(SC3[,c(1:2)],SC3[,grep(colnames(SC3),pattern='TGFB')])

setwd('/home/guidantoniomt/pseudospace')

source('functions_plot_pseudospace.R')

norm_expr_data<-function(FM,pseudo_expr){
           
	    FM2 <- FM/estimateSizeFactorsForMatrix(FM)
            if (is.null(pseudo_expr)) 
                pseudo_expr <- 1
            FM2 <- FM2 + pseudo_expr
            FM2 <- log2(FM2)
	
	return(FM2)
}

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

# 
# upload signatures
# 

# mock_emt_gs<-read.delim(file='spontaneous_EMT_trajectory_test.txt',stringsAsFactors=F)[,1]
# mock_tgfb_gs<-read.delim(file='TGFB_EMT_trajectory_test.txt',stringsAsFactors=F)[,1]

# mock_norm2_gs<-mock_norm2[mock_norm2[,1]%in%mock_emt_gs,]
# tgfb_norm2_gs<-tgfb_norm2[tgfb_norm2[,1]%in%mock_tgfb_gs,]

# mock_symbol_gs<-as.character(mock_norm2_gs[,2])
# tgfb_symbol_gs<-as.character(tgfb_norm2_gs[,2])

# 
# upload TCGA
# 
setwd('/home/guidantoniomt/analysis_results')
tcga<-data.frame(fread('TCGA_combat_tumor_type_correction.txt'))

# UPLOAD MET500 DATA-SET
setwd('/home/guidantoniomt/analysis_results')
load('MET500_rna_seq_log2_fpkm.ALLGENES.RData')

# look here: https://www.biostars.org/p/273370/
# and here: https://www.biostars.org/p/322114/
# https://www.biostars.org/p/160989/

met500_unlog<-2^met500_rnaseq2[,-c(1:2)]

fpkmToTpm <- function(fpkm) {

(fpkm / sum(fpkm)  *(1e6))

}


met500_convert<-log2(apply(met500_unlog[,-c(1:2)],1,fpkmToTpm)+1)
colnames(met500_convert)<-met500_rnaseq2[,2]


met500_convert<-log2(apply(met500_unlog[,-c(1:2)],1,fpkmToTpm)+1)
colnames(met500_convert)<-met500_rnaseq2[,2]
met500_convert<-t(met500_convert)

# get the commons gene-names between TCGA and MET500
tcga_genes<-tcga[,1]
met500_genes<-colnames(met500_convert)
common_genes<-intersect(tcga_genes,met500_genes)

tcga<-tcga[which(tcga[,1]%in%common_genes),]

met500_convert2<-met500_convert[which(rownames(met500_convert)%in%common_genes),]
colnames(met500_convert2)<-as.character(sapply(strsplit(colnames(met500_convert2),split='\\.'),'[[',1))
met500_convert2<-data.frame(genes=rownames(met500_convert2),met500_convert2)

met500_convert2<-aggregate(.~genes,met500_convert2,mean)
tcga<-as.data.frame(setDT(tcga)[, lapply(.SD, mean), by = .(genes)])

input_ge<-merge(tcga,met500_convert2,by='genes')
rownames(input_ge)<-input_ge[,1]

# tcga_mock<-tcga[tcga[,1]%in%mock_symbol_gs,]
# tcga_tgfb<-tcga[tcga[,1]%in%tgfb_symbol_gs,]

input_mock<-merge(input_ge,mock_norm2[,-1],by.x='genes',by.y='genes.gene_short_name')
input_mock[input_mock<0]<-0 # trim negative values
input_tgfb<-merge(input_ge,tgfb_norm2[,-1],by.x='genes',by.y='genes.gene_short_name')
input_tgfb[input_tgfb<0]<-0 # trim negative values



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

		leninner<-length(grep(rownames(scrna_pca)[res$nn.index],pattern='inner'))
		lenouter<-length(grep(rownames(scrna_pca)[res$nn.index],pattern='outer'))

		}
return(knn.list)
			
}

combat_mock_tcga_knn<-knnProject(combat_mock_filter,knn=10,num_dim=25)

combat_tgfb2_tcga_knn<-knnProject(combat_tgfb2_filter,knn=10,num_dim=25)


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
knn_df_tcga_met500_tgfb<-knn_df_tcga_tgfb

setwd('/home/guidantoniomt/pseudospace/TCGA_MET500')

write.table(knn_df_tcga_met500_tgfb,file='KNN_projection_TCGA_MET500_to_MCF10A_tgfb_treated_cells.txt',row.names=F,quote=F,sep='\t')
write.table(knn_df_tcga_met500_mock,file='KNN_projection_TCGA_MET500_to_MCF10A_mock_treated_cells.txt',row.names=F,quote=F,sep='\t')


save.image(file="KNN_projection_TCGA_MET500_to_MCF10A.RData")

setwd('/home/guidantoniomt/analysis_results')

clinical_data<-read.delim(file='clinical_table_TCGA_MET500_with_RECLASSIFICATION_withmedian.txt',stringsAsFactors=F)
clinical_data_sub<-clinical_data[,c('samples','type')]
colnames(clinical_data_sub)[1:2]<-c('patients','tumors')

knn_df_tcga_met500_mock$patients<-as.character(knn_df_tcga_met500_mock$patients)
idx<- grep(knn_df_tcga_met500_mock$patients,pattern='TCGA')
patients_mock<-grep(knn_df_tcga_met500_mock$patients,pattern='TCGA',value=T)
knn_df_tcga_met500_mock[idx,'patients']<-gsub(make.names(sapply(strsplit(as.character(patients_mock),split='\\.'),'[',4),unique=T),pattern='\\.',replacement='_')

knn_df_tcga_met500_tgfb$patients<-as.character(knn_df_tcga_met500_tgfb$patients)
idx<- grep(knn_df_tcga_met500_tgfb$patients,pattern='TCGA')
patients_tgfb<-grep(knn_df_tcga_met500_tgfb$patients,pattern='TCGA',value=T)
knn_df_tcga_met500_tgfb[idx,'patients']<-gsub(make.names(sapply(strsplit(as.character(patients_tgfb),split='\\.'),'[',4),unique=T),pattern='\\.',replacement='_')

knn_df_tcga_met500_mock<-merge(knn_df_tcga_met500_mock,clinical_data_sub,by='patients')
knn_df_tcga_met500_tgfb<-merge(knn_df_tcga_met500_tgfb,clinical_data_sub,by='patients')
colnames(knn_df_tcga_met500_mock)[c(4)]<-c('tumors')
colnames(knn_df_tcga_met500_tgfb)[c(4)]<-c('tumors')

setwd('/home/guidantoniomt/pseudospace/TCGA_MET500')

 pdf("KNN_projection_TCGA_MET500_to_MCF10A_mock_treated_cells.pdf",height=20)
 p<-ggplot(knn_df_tcga_met500_mock, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")+facet_wrap(~tumors, ncol = 1, scales = "free_y") +xlab("Spontaneous EMT") +
 ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 11),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts() 
 print(p)
 dev.off()

 pdf("KNN_projection_TCGA_MET500_to_MCF10A_tgfb_treated_cells.pdf",height=20)
 p<-ggplot(knn_df_tcga_met500_tgfb, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C") +facet_wrap(~tumors, ncol = 1, scales = "free_y") +xlab("TGFB EMT") +
 ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 11),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts()
 print(p)
 dev.off()

pdf("KNN_projection_TCGA_MET500_to_MCF10A_mock_treated_cells.ALLCANCER.pdf",height=5)
p<-ggplot(knn_df_tcga_met500_mock, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")+xlab("Spontaneous EMT") +
ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 15),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts()
print(p)
dev.off()

pdf("KNN_projection_TCGA_MET500to_MCF10A_tgfb_treated_cells.ALLCANCER.pdf",height=5)
p<-ggplot(knn_df_tcga_met500_tgfb, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C")+xlab("TGFB EMT") +
ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 15),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts()
print(p)
dev.off()

colnames(knn_df_tcga_met500_mock)[c(3)]<-c('mock')
colnames(knn_df_tcga_met500_tgfb)[c(3)]<-c('tgfb')

# Create a countour plot
input_cp<-data.frame(knn_df_tcga_met500_mock,knn_df_tcga_met500_tgfb)

input_cp$groups<-(rep(0,nrow(input_cp)))
input_cp[grep(input_cp$tumors,pattern='TCGA'),'groups']<-'TCGA'
input_cp[grep(input_cp$tumors,pattern='TCGA',invert=T),'groups']<-'MET500'

input_cp$groups<-as.factor(input_cp$groups)

library(gridExtra)
library(ggforce)

pdf('KNN_projection_TCGA_MET500_to_MCF10A_mock_vs_tgfb_v2.pdf',width=15)
#mod5 <- densityMclust(input_cp[,c('mock','tgfb')])
#input_plot<-input_cp[,c(3,4)]
#plot(mod5, what = "density", type = "hdr",grid=200)
mock<-ggplot(knn_df_tcga_mock, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")
tgfb<-ggplot(knn_df_tcga_tgfb, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C")+ coord_flip()

sp <- ggplot(input_cp, aes(x=mock, y=tgfb,col=groups))+geom_point()+facet_wrap(~groups)

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.border = element_blank(),
   panel.background = element_blank(),
   axis.title.x = element_blank(),
   axis.title.y = element_blank(),
   axis.text.x = element_blank(),
   axis.text.y = element_blank(),
   axis.ticks = element_blank()
     )

grid.arrange(mock,blankPlot,sp, tgfb,
        ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

dev.off()

pdf('densities_mock_vs_tgfb_all_cancers.pdf')
for(i in 1:6){
sp <- ggplot(input_cp, aes(x=mock, y=tgfb)) +geom_point()+stat_density_2d(binwidth=5e-04)+facet_wrap_paginate(~tumors, ncol = 3,nrow=3,page=i)
print(sp)
}
dev.off()






# NOTE: our definition of "reference" and "query" are reversed w.r.t what DTW uses.
align_cells <- function(query_cell_matrix, ref_cell_matrix,
                        step_pattern=rabinerJuangStepPattern(3, "d", smoothed=TRUE),
                        open.end=FALSE,
                        open.begin=FALSE)
{
  # we want the query and reference to have the same number of cells, so downsample the larger
  # one to match the smaller one

  # GMT: get the common genes between ref_cell_matrix and query_cell_matrix
  common_genes<-intersect(ref_cell_matrix[,1],query_cell_matrix[,1])
  
  ref_cell_matrix<-ref_cell_matrix[ref_cell_matrix[,1]%in%common_genes,]
  query_cell_matrix<-query_cell_matrix[query_cell_matrix[,1]%in%common_genes,]  
  
  # GMT: the number of genes che be different because they are repeated
  ref_cell_matrix<-as.data.frame(setDT(ref_cell_matrix)[, lapply(.SD, mean), by = .(genes)])
  query_cell_matrix<-as.data.frame(setDT(query_cell_matrix)[, lapply(.SD, mean), by = .(genes)])
  
  # GMT: rownames are gene names
  rownames(ref_cell_matrix)<-ref_cell_matrix[,1]
  rownames(query_cell_matrix)<-query_cell_matrix[,1]

  # GMT: genes on the rows, samples on the columns
  ref_cell_matrix<-ref_cell_matrix[,-1]
  query_cell_matrix<-query_cell_matrix[,-1]

  if (ncol(ref_cell_matrix) < ncol(query_cell_matrix)){
    ref_cell_matrix_sub <- ref_cell_matrix
    sampled_idxs <- c(1, sort(sample(seq(2,ncol(query_cell_matrix)-1), ncol(ref_cell_matrix)-2)), ncol(query_cell_matrix))
    query_cell_matrix_sub <- query_cell_matrix[,sampled_idxs]
  }else if (ncol(ref_cell_matrix) > ncol(query_cell_matrix)){
    query_cell_matrix_sub <- query_cell_matrix
    sampled_idxs <- c(1, sort(sample(seq(2,ncol(ref_cell_matrix)-1), ncol(query_cell_matrix)-2)), ncol(ref_cell_matrix))
    ref_cell_matrix_sub <- ref_cell_matrix[,sampled_idxs]
  }else{
    query_cell_matrix_sub <- query_cell_matrix
    ref_cell_matrix_sub <- ref_cell_matrix
  }

  ref_cell_matrix_sub <- t(scale(t(ref_cell_matrix_sub)))
  query_cell_matrix_sub <- t(scale(t(query_cell_matrix_sub)))
  

  valid_genes <- rowSums(is.na(ref_cell_matrix_sub)) == 0 & rowSums(is.na(query_cell_matrix_sub)) == 0
  ref_cell_matrix_sub <- ref_cell_matrix_sub[valid_genes,]
  query_cell_matrix_sub <- query_cell_matrix_sub[valid_genes,]

  print (dim(ref_cell_matrix_sub))
  #dist_matrix <- (1 - cor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs"))^2
  print('correlation')
  #dist_matrix <- (1 - cor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs"))
 
  dist_matrix <- (1 - bigcor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs",size=50))

  cell_alignment_dtw <- dtw(dist_matrix,
                            step.pattern=step_pattern,
                            open.end=open.end,
                            open.begin=open.begin,
                            keep=T,
                            keep.internals=T)
  cell_alignment_dtw
}


warp_pseudotime <- function(ref_cds, query_cds, alignment_dtw)
{
 ref_pseudotimes <- pData(ref_cds)$Pseudotime
 names(ref_pseudotimes) <- row.names(pData(ref_cds))
 ref_pseudotimes <- sort(ref_pseudotimes)

 query_pseudotimes <- pData(query_cds)$Pseudotime
 names(query_pseudotimes) <- row.names(pData(query_cds))
 query_pseudotimes <- sort(query_pseudotimes)

 query_alignment_time <- warp(alignment_dtw, index.reference=T)

 #print (query_alignment_time)
 pData(ref_cds)$Alignment_Pseudotime <- pData(ref_cds)$Pseudotime
 pData(query_cds)$Alignment_Pseudotime <- NA

 #names(query_alignment_time) <- rownames(alignment_dtw$localCostMatrix)

 #pData(query_cds)[names(query_alignment_time),]$Alignment_Pseudotime <- query_alignment_time

 fc <- approxfun(seq(1, 101, by=1), query_alignment_time, rule=2)
 pData(query_cds)$Alignment_Pseudotime  <- fc(query_cds$Pseudotime)
 return(list(ref_cds=ref_cds, query_cds=query_cds))

 # BJ_MYO_aligned <- BJ_MYO_selected
 # pData(BJ_MYO_aligned)[colnames(smoothed_BJ_MYO_exprs),]$Pseudotime <- BJ_MYO_alignment_time
 # pData(BJ_MYO_selected)[colnames(smoothed_BJ_MYO_exprs),]$Pseudotime <- BJ_MYO_alignment_time
}


# 
# Alignment pseudospace
#
library(WGCNA)
library(data.table)
library(propagate)

load('KNN_projection_TCGA_to_MCF10A.RData')

ref_cell_matrix<-combat_mock_filter
ps<-knn_df_tcga_mock

query_cell_matrix<-combat_tgfb2_filter
ps<-knn_df_tcga_tgfb


res_alignment<-align_cells(query_cell_matrix, ref_cell_matrix,
                        step_pattern=rabinerJuangStepPattern(3, "d", smoothed=TRUE),
                        open.end=FALSE,
                        open.begin=FALSE)






























# Plot Mock 95 genes
load('KNN_projection_TCGA_to_MCF10A.RData')

setwd('/home/guidantoniomt/module_analysis_prediction_metastatic')
clusters<-read.delim(file='TCGA_SEGS_95genes.communities_clustering_methods.WITH_no_metastatic_GENES.txt')
gs<- gsub(unlist(strsplit(as.character(clusters[-c(13,14),2]),split='\\,')),pattern=' ',replacement='')


setwd('/home/guidantoniomt/analysis_results')
tab_clinical<-read.delim(file='clinical_table_TCGA_MET500_with_RECLASSIFICATION_withmedian.txt',stringsAsFactors=F)

mtx_exp<-combat_mock_filter
ps<-knn_df_tcga_mock

plotGenesWithPseudospace(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_mock_treated_cells.95genes.cancer')

mtx_exp<-combat_tgfb2_filter
ps<-knn_df_tcga_tgfb

plotGenesWithPseudospace(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.95genes.cancer')

# Plot DDR genes
setwd('/home/guidantoniomt/pathway_activities')

table_pathway<-read.delim(file='DDRpathways_full_table.txt',stringsAsFactors=F)
gs<-table_pathway[,6]

setwd("/home/guidantoniomt/pseudospace")

mtx_exp<-combat_mock_filter
ps<-knn_df_tcga_mock

plotGenesWithPseudospace(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_mock_treated_cells.DDRgenes.cancer')

mtx_exp<-combat_tgfb2_filter
ps<-knn_df_tcga_tgfb

plotGenesWithPseudospace(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.DDRgenes.cancer')


# Plot marker metastatic
gs<-c('CDH1','CRB3','DSP','CDH2','FN1','VIM')

mtx_exp<-combat_mock_filter
ps<-knn_df_tcga_mock

plotGenesWithPseudospace(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_mock_treated_cells.met_markersgenes.cancer', clusters=4)

mtx_exp<-combat_tgfb2_filter
ps<-knn_df_tcga_tgfb

plotGenesWithPseudospace(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.met_markersgenes.cancer',clusters=4)

setwd('/home/guidantoniomt/analysis_results')
load('TCGA_SEGDATA_with_status_Alison_and_CXI_17_11_2019.RData')

setwd('/home/guidantoniomt/module_analysis_prediction_metastatic')

GENCODE<-as.data.frame(fread(file="gencode.v31.annotation.gtf",header=F,skip='##', sep="\t",stringsAsFactors=F))
gencode_genes<-GENCODE[GENCODE[,3]=='gene',]
ENSG<-sapply(strsplit(sapply((strsplit(gencode_genes[,9],split=';')),'[[',1),split=' '),'[[',2)
GENE_NAME<-sapply(strsplit(sapply(strsplit(gsub(gencode_genes[,9],pattern="\"",replacement="",fixed=T),split=';'),'[[',3),' '),'[[',3)
gencode_genes_parse<-data.frame(chr=gencode_genes[,1],start=gencode_genes[,4],end=gencode_genes[,5],ENSG=ENSG,GENE_NAME=GENE_NAME,stringsAsFactors=F)
gencode_genes_parse[,1]<-gsub(gencode_genes_parse[,1],pattern='chr',replacement='')

annotation_file=gencode_genes_parse

library(GenomicRanges)
library(reshape2)

dictionary_copy_number<-data.frame(RES_ALL_CANCERS_with_effect3[,-c(1:2)],samples=RES_ALL_CANCERS_with_effect3[,1],cancer=RES_ALL_CANCERS_with_effect3[,2],stringsAsFactors=F)

# Plot Mock 95 genes

setwd('/home/guidantoniomt/module_analysis_prediction_metastatic')
clusters<-read.delim(file='TCGA_SEGS_95genes.communities_clustering_methods.WITH_no_metastatic_GENES.txt')
gs<- gsub(unlist(strsplit(as.character(clusters[-c(13,14),2]),split='\\,')),pattern=' ',replacement='')

setwd('/home/guidantoniomt/analysis_results')
tab_clinical<-read.delim(file='clinical_table_TCGA_MET500_with_RECLASSIFICATION_withmedian.txt',stringsAsFactors=F)

input_cnv<-createMatrixFromSegFile(dictionary_copy_number,annotation_file,gs)

mtx_exp<-input_cnv
ps<-knn_df_tcga_mock
ps<-ps[grep(ps[,2],pattern='TCGA'),]
ps[,2]<-sapply(strsplit(as.character(ps[,2]),'\\.'),'[[',4)

setwd("/home/guidantoniomt/pseudospace")

plotGenesWithPseudospaceCNV(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_mock_treated_cells.95genes.cnv.cancer',clusters=9)

ps<-knn_df_tcga_tgfb
ps<-ps[grep(ps[,2],pattern='TCGA'),]
ps[,2]<-sapply(strsplit(as.character(ps[,2]),'\\.'),'[[',4)

plotGenesWithPseudospaceCNV(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.95genes.cnv.cancer',clusters=9)

# Plot metastatic markers
gs<-c('CDH1','CRB3','DSP','CDH2','FN1','VIM')
setwd('/home/guidantoniomt/analysis_results')
tab_clinical<-read.delim(file='clinical_table_TCGA_MET500_with_RECLASSIFICATION_withmedian.txt',stringsAsFactors=F)

input_cnv<-createMatrixFromSegFile(dictionary_copy_number,annotation_file,gs)

mtx_exp<-input_cnv
ps<-knn_df_tcga_mock
ps<-ps[grep(ps[,2],pattern='TCGA'),]
ps[,2]<-sapply(strsplit(as.character(ps[,2]),'\\.'),'[[',4)

setwd("/home/guidantoniomt/pseudospace")

plotGenesWithPseudospaceCNV(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_mock_treated_cells.met_markersgenes.cnv.cancer',clusters=5)

ps<-knn_df_tcga_tgfb
ps<-ps[grep(ps[,2],pattern='TCGA'),]
ps[,2]<-sapply(strsplit(as.character(ps[,2]),'\\.'),'[[',4)

plotGenesWithPseudospaceCNV(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.met_markersgenes.cnv.cancer',clusters=5)

# Plot DDR genes
setwd('/home/guidantoniomt/pathway_activities')

table_pathway<-read.delim(file='DDRpathways_full_table.txt',stringsAsFactors=F)
gs<-table_pathway[,6]

input_cnv<-createMatrixFromSegFile(dictionary_copy_number,annotation_file,gs)

mtx_exp<-input_cnv
ps<-knn_df_tcga_mock
ps<-ps[grep(ps[,2],pattern='TCGA'),]
ps[,2]<-sapply(strsplit(as.character(ps[,2]),'\\.'),'[[',4)

setwd("/home/guidantoniomt/pseudospace")

plotGenesWithPseudospaceCNV(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_mock_treated_cells.DDRgenes.cnv.cancer',clusters=9)

ps<-knn_df_tcga_tgfb
ps<-ps[grep(ps[,2],pattern='TCGA'),]
ps[,2]<-sapply(strsplit(as.character(ps[,2]),'\\.'),'[[',4)

plotGenesWithPseudospaceCNV(mtx_exp=mtx_exp,gs=gs,ps=ps,tab_clinical=tab_clinical,clinical_feature='tumor_stage',trim=0.1,output='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.DDRgenes.cnv.cancer',clusters=9)

