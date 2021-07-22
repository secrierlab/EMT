library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)


setwd('/data/pseudospace/')

norm_expr_data<-function(FM,pseudo_expr){
           
	    FM2 <- FM/estimateSizeFactorsForMatrix(FM)
            if (is.null(pseudo_expr)) 
                pseudo_expr <- 1
            FM2 <- FM2 + pseudo_expr
            FM2 <- log2(FM2)
	
	return(FM2)
}

cds.list <- readRDS("/data/pseudospace/pseudospace_processed_trajectories_cds.list.rds")

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

# 
# upload TCGA
# 
load('/data/pseudospace/TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData')

# UPLOAD MET500 DATA-SET
load('/data/pseudospace/MET500_rna_seq_log2_fpkm.ALLGENES.RData')

# look here: https://www.biostars.org/p/273370/
# and here: https://www.biostars.org/p/322114/
# https://www.biostars.org/p/160989/

met500_unlog<-2^met500_rnaseq2[,-c(1:2)]

fpkmToTpm <- function(fpkm) {

(fpkm / sum(fpkm)  *(1e6))

}

met500_convert<-log2(apply(met500_unlog[,-c(1:2)],1,fpkmToTpm)+1)
colnames(met500_convert)<-met500_rnaseq2[,2]
met500_convert<-t(met500_convert)
print(range(met500_convert))

# get the commons gene-names between TCGA and MET500
tcga_genes<-TCGA_GEXP_ALL[,1]
met500_genes<-rownames(met500_convert)
common_genes<-intersect(tcga_genes,met500_genes)

tcga<-TCGA_GEXP_ALL[which(TCGA_GEXP_ALL[,1]%in%common_genes),]

met500_convert2<-met500_convert[which(rownames(met500_convert)%in%common_genes),]
colnames(met500_convert2)<-as.character(sapply(strsplit(colnames(met500_convert2),split='\\.'),'[[',1))
met500_convert2<-data.frame(genes=rownames(met500_convert2),met500_convert2)

met500_convert2<-as.data.frame(setDT(met500_convert2)[, lapply(.SD, mean), by = .(genes)])
tcga<-as.data.frame(setDT(tcga)[, lapply(.SD, mean), by = .(gene_symbol)])
tcga[,-1]<-log(tcga[,-1]+1,2)
  
print(range(tcga[,-1]))

input_ge<-merge(tcga,met500_convert2,by.x='gene_symbol',by.y='genes')
rownames(input_ge)<-input_ge[,1]

input_mock<-merge(input_ge,mock_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
input_mock[input_mock<0]<-0 # trim negative values
input_tgfb<-merge(input_ge,tgfb_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
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


library(factoextra)

pdf('TCGA_MET500_MOCK_PrincipalComponents.pdf',width=15)
pca_mock <- prcomp(t(combat_mock_filter[,-1]))
p1<-fviz_eig(pca_mock,ncp=50,addlabels=T)
print(p1)
#with 30 PC i explain 33.% of variance 30 is also the point in which the variance become 0.1
variance_explained_mock<-sum(get_eig(pca_mock)[1:30,2])
dev.off()

pdf('TCGA_MET500_TGFB_PrincipalComponents.pdf',width=15)
pca_tgfb <- prcomp(t(combat_tgfb2[,-1]))
p2<-fviz_eig(pca_tgfb,ncp=50,addlabels=T)
print(p2)
#with 30 PC i explain 33.% of variance 30 is also the point in which the variance become 0.1
variance_explained_tgfb<-sum(get_eig(pca_tgfb)[1:30,2])
dev.off()


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

setwd("/home/guidantoniomt/pseudospace")
write.table(input_cp_final,file='KNN_projection_TCGA_MET500_to_MCF10A_treated_cells.txt',row.names=F,quote=F,sep='\t')


# 
# Compute the EMT scores of the samples along the MOCK psuedotime
# 


setwd("/home/guidantoniomt/pseudospace/HMM")
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)

epi_genes<-intersect(input_ge[,1], markers_genes_read[markers_genes_read$status=="Epithelial_marker",1])
mes_genes<-intersect(input_ge[,1], markers_genes_read[markers_genes_read$status=="Mesenchymal_marker",1])

# Compute an EMT score of the patients
# VIM + CDH2 + FOXC2 + SNAI1 + SNAI2 + TWIST1 + FN1 + ITGB6 + MMP2 + MMP3 + MMP9 + SOX10 + GCS − CDH1 − DSP − OCLN
# https://www.nature.com/articles/s41598-018-21061-1

input_ge2<-t(input_ge[input_ge[,1]%in%c(epi_genes,mes_genes),which(colnames(input_ge)%in%input_cp_final$patients)])

zscore<-function(x){(x-mean(x))/sd(x)}
input_ge3<-t(apply(input_ge2,2,zscore))

all_scores_emt<-NULL

for(iemt in 1:ncol(input_ge3)){
  
  mean_epi<-mean(input_ge3[epi_genes,iemt])
  mean_mes<-mean(input_ge3[mes_genes,iemt])
  
  score_emt<-mean_mes-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

df_emt_samples<-data.frame(ID=colnames(input_ge3),EMT_scores=all_scores_emt)
df_emt_samples$ID<-as.character(df_emt_samples$ID)
input_cp_final$patients<-as.character(input_cp_final$patients)

input_cp_final2<-merge(x=input_cp_final[,-c(4:5)],y=df_emt_samples,by.x="patients",by.y="ID")

write.table(input_cp_final2,file='KNN_projection_TCGA_MET500_to_MCF10A_treated_cells_withEMT.txt',row.names=F,quote=F,sep='\t')

input_cp_final2$new_ann<-rep("no",nrow(input_cp_final2))
input_cp_final2$new_ann[grep(input_cp_final2$patients,pattern="TCGA")]<-"TCGA"
input_cp_final2$new_ann[grep(input_cp_final2$patients,pattern="TCGA",invert=T)]<-"MET500"
input_cp_final2$new_ann<-as.factor(input_cp_final2$new_ann)


input_cp_final2$emt_status<-ifelse(input_cp_final2$EMT_scores>0,"High_EMT","Low_EMT")
input_cp_final2$pseudospace_status<-ifelse(input_cp_final2$mock_pseudospace<=15,"Late_Pseudotime","Early_Pseudotime")

library(ggplot2)
library(ggridges)
library(corrplot)

setwd("/data/pseudospace")

pdf(paste("TCGA_MET500_pseudotime_vs_emt_scores_mock.pdf",sep="."),width=12,pointsize=12)
p<-ggplot(input_cp_final2, 
          aes(x=mock_pseudospace, y=EMT_scores,color=new_ann)) +
          geom_point(shape=18,size=2,alpha=0.5)+scale_x_reverse()+ 
          scale_color_manual(values=c('#ee0c0c','#B5B7FF'))+
          geom_smooth(aes(group=new_ann),method = "lm", formula = y ~ poly(x, 10),se=TRUE, linetype="dashed",color="blue")+
          geom_hline(yintercept=0, linetype="dashed", color = "black")+theme_classic()
print(p)

input_cp_final3<-input_cp_final2[input_cp_final2$new_ann%in%"MET500",]
chisq <- chisq.test(table(input_cp_final3$pseudospace_status,input_cp_final3$emt_status)[2:1,])
contrib <- 100*chisq$residuals^2/chisq$statistic
corrplot(chisq$residuals, is.cor = FALSE)
corrplot(contrib, is.cor = FALSE)
dev.off()


# 
# Compare the EMT scores of the MET500 samples and the samples defined using HMM
# 

setwd("/data/pseudospace/HMM")
tab_hmm<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)
samples_with_HMM_and_MET500<-merge(tab_hmm,input_cp_final2,by.x="samples",by.y="patients",all.x=T,all.y=T)
samples_with_HMM_and_MET500$HMM_states[which(is.na(samples_with_HMM_and_MET500$HMM_states))]<-"MET500"
samples_with_HMM_and_MET500$HMM_states<-gsub(samples_with_HMM_and_MET500$HMM_states,pattern="2",replacement="epi")
samples_with_HMM_and_MET500$HMM_states<-gsub(samples_with_HMM_and_MET500$HMM_states,pattern="1",replacement="pEMT")
samples_with_HMM_and_MET500$HMM_states<-gsub(samples_with_HMM_and_MET500$HMM_states,pattern="3",replacement="mes")
samples_with_HMM_and_MET500$HMM_states<-factor(samples_with_HMM_and_MET500$HMM_states,levels=c("epi","pEMT","mes","MET500"))

setwd("/data/pseudospace")
library(ggpubr)

pdf("Compare_EMT_scores_primary_HMMstates_with_MET500.pdf")
p<-ggviolin(samples_with_HMM_and_MET500,
            x = "HMM_states",
            y = "EMT_scores",
            fill = "HMM_states",
            palette  = c("#56b4e9","#E69f00","#EE0C0C","#A10070"),
            add = "boxplot",
            add.params = list(fill = "white"))+
  stat_compare_means()+geom_hline(yintercept=0)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

dev.off()

# 
# Compare the pseudotime values for each category HMM
# 

setwd("/data/pseudospace/HMM")
tab_hmm_test<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)

#tab_hmm_test$pseudospace<-max(tab_hmm_test$pseudospace)-tab_hmm_test$pseudospace
tab_hmm_test$pseudospace<-1-(tab_hmm_test$pseudospace/max(tab_hmm_test$pseudospace))

# 
# Compare the hypoxia scores of the MET500 samples with the hypoxia scores of the TCGA samples
# 

setwd("/data/pseudospace/pseudospace")

pdf("Pseudotime_values_for_HMMstate.pdf")
p<-ggviolin(tab_hmm_test,
            x = "biological_states",
            y = "pseudospace",
            fill = "biological_states",
            palette  = c("#56b4e9","#E69f00","#EE0C0C","#A10070"),
            add = c("boxplot"),
            add.params = list(fill = "white"))+
  stat_compare_means()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()
