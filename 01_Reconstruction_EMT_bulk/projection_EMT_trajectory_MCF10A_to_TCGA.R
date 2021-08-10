library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

source('functions_plot_pseudospace.R')
current_dir<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

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
setwd(input_dir)
load('TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData')

tcga<-TCGA_GEXP_ALL
tcga[,-1]<-log(tcga[,-1]+1,2)

input_mock<-merge(tcga,mock_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')
input_tgfb<-merge(tcga,tgfb_norm2[,-1],by.x='gene_symbol',by.y='genes.gene_short_name')

# 
# Combat Correction
# 
library(sva)

batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(mock_norm2[,-c(1:2)]))))

combat_mock = ComBat(dat=as.matrix(input_mock[,-1]), batch=batch)
combat_mock2<-data.frame(genes=input_mock[,1],combat_mock)
combat_mock2[,1]<-as.character(combat_mock2[,1])
combat_mock2[combat_mock2<0]<-0

batch<-factor(c(rep(1,ncol(tcga[,-1])),rep(2,ncol(tgfb_norm2[,-c(1:2)]))))

combat_tgfb = ComBat(dat=as.matrix(input_tgfb[,-1]), batch=batch)
combat_tgfb2<-data.frame(genes=input_tgfb[,1],combat_tgfb)
combat_tgfb2[,1]<-as.character(combat_tgfb2[,1])
combat_tgfb2[combat_tgfb2<0]<-0

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
	
		tcga_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='\\.TCGA'),]
		scrna_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='\\.TCGA',invert=T),]

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
irlba_res <- prcomp(t(combat_mock_filter[,-1]))

pdf('pca_mock.pdf',width=15)
fviz_eig(irlba_res,ncp=30,addlabels=T)
dev.off()

sum(get_eig(irlba_res)[1:25,2])

irlba_res <- prcomp(t(combat_tgfb2_filter[,-1]))

pdf('pca_tgfb.pdf')
fviz_eig(irlba_res,ncp=30,addlabels=T)
dev.off()

sum(get_eig(irlba_res)[1:25,2])


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

tumors<-sapply(strsplit(names(combat_mock_tcga_knn),split='\\.'),'[[',1)
knn_df_tcga_mock<-data.frame(tumors=tumors,patients=names(combat_mock_tcga_knn),pseudospace=knn_mock_pseudospace,ratio=ratio_mock_outer_inner,cor_mock=knn_mock_corr)

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

tumors<-sapply(strsplit(names(combat_tgfb2_tcga_knn),split='\\.'),'[[',1)
knn_df_tcga_tgfb<-data.frame(tumors=tumors,patients=names(combat_tgfb2_tcga_knn),pseudospace=knn_tgfb_pseudospace,ratio=ratio_tgfb_outer_inner,cor_tgfb=knn_tgfb_corr)

setwd(output_dir)

pdf("KNN_projection_TCGA_to_MCF10A_mock_treated_cells.pdf",height=15)
# Generate plots of the density of projected HSNCC cells across pseudospace
p<-ggplot(knn_df_tcga_mock, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2") +facet_wrap(~tumors, ncol = 1, scales = "free_y") +xlab("Spontaneous EMT") +
ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 11),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts() 
print(p)
dev.off()

pdf("KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.pdf",height=15)
# Generate plots of the density of projected HSNCC cells across pseudospace
p<-ggplot(knn_df_tcga_tgfb, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C") +facet_wrap(~tumors, ncol = 1, scales = "free_y") +xlab("TGFB EMT") +
ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 11),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts()
print(p)
dev.off()

pdf("KNN_projection_TCGA_to_MCF10A_mock_treated_cells.ALLCANCER.pdf",height=5)
# Generate plots of the density of projected HSNCC cells across pseudospace
p<-ggplot(knn_df_tcga_mock, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")+xlab("Spontaneous EMT") +
ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 15),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts()
print(p)
dev.off()

pdf("KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells.ALLCANCER.pdf",height=5)
# Generate plots of the density of projected HSNCC cells across pseudospace
p<-ggplot(knn_df_tcga_tgfb, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C")+xlab("TGFB EMT") +
ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 15),axis.text.y = element_text(size = 3)) + monocle:::monocle_theme_opts()
print(p)
dev.off()

setwd(output_dir)

write.table(knn_df_tcga_tgfb,file='KNN_projection_TCGA_to_MCF10A_tgfb_treated_cells_no_correction_primarytumor.txt',row.names=F,quote=F,sep='\t')
write.table(knn_df_tcga_mock,file='KNN_projection_TCGA_to_MCF10A_mock_treated_cells_no_correction_primarytumor.txt',row.names=F,quote=F,sep='\t')


# Create a countour plot

input_cp<-data.frame(knn_df_tcga_mock,knn_df_tcga_tgfb[,c(3:5)])
colnames(input_cp)[3:8]<-c('mock','mock_ratio','cor_mock','tgfb','tgfb_ratio','cor_tgfb')


setwd(input_dir)
tab_clinical<-read.delim(file='clinical_table_TCGA_MET500_with_RECLASSIFICATION_withmedian.txt',stringsAsFactors=F)
tab_clinical<-tab_clinical[tab_clinical$assay=='TCGA',]

sub_clinical2<-tab_clinical[,c('samples','tumor_stage','metastatic_annotation','NEW_METASTATIC')]
sub_clinical2[which(sub_clinical2$metastatic_annotation==1),'tumor_stage']<-'metastatic'

high_stage<-which(sub_clinical2$tumor_stage=='stage iii' | sub_clinical2$tumor_stage=='stage iiia'| sub_clinical2$tumor_stage=='stage iiib'| sub_clinical2$tumor_stage=='stage iiic' | sub_clinical2$tumor_stage=='stage iv'| sub_clinical2$tumor_stage=='stage iva'| sub_clinical2$tumor_stage=='stage ivb'| sub_clinical2$tumor_stage=='stage ivc')

sub_clinical2[high_stage,'tumor_stage']<-'high_stage'
sub_clinical2[-which(sub_clinical2$tumor_stage%in%c('metastatic','high_stage','not reported','--','is')),'tumor_stage']<-'low_stage'

sub_clinical2$tumor_stage<-gsub(sub_clinical2$tumor_stage,pattern='is',replacement='not reported')
sub_clinical2$tumor_stage<-gsub(sub_clinical2$tumor_stage,pattern='\\--',replacement='not reported')

colnames(sub_clinical2)[1]<-'patients'

input_cp$patients2<-input_cp$patients
input_cp$patients<-gsub(make.names(sapply(strsplit(as.character(input_cp$patients),split='\\.'),'[',4),unique=T),pattern='\\.',replacement='_')
input_cp<-merge(input_cp,sub_clinical2[,c(1,2)],by='patients')

input_cp$tumor_stage<-as.factor(input_cp$tumor_stage)

setwd('/home/data/pseudospace/tcga_nocorrect_primary_tumors/23_04_2020_threegroups')

save.image(file="KNN_projection_TCGA_to_MCF10A_no_correction_primarytumor.RData")

setwd(output_dir)
library(ggpubr)

png('MOCK_correlation_TCGA_samples_with_scRNA_seq.png', height=6,width=15,units='in',res=300)
ggboxplot(input_cp,x='tumor_stage',y='cor_mock',col='tumor_stage',palette=c('#E69F00','#56B4E9','#999999','#0DB275'))+facet_wrap(~mock_ratio)
dev.off()


png('TGFB_correlation_TCGA_samples_with_scRNA_seq.png', height=6,width=15,units='in',res=300)
ggboxplot(input_cp,x='tumor_stage',y='cor_tgfb',col='tumor_stage',palette=c('#E69F00','#56B4E9','#999999','#0DB275'))+facet_wrap(~tgfb_ratio)
dev.off()

library(rootSolve)

set.seed(123)
mixmdl = normalmixEM((input_cp$tgfb/max(input_cp$tgfb))*100,k=3,fast=TRUE,maxit=500000)
mixmdl$mu

f <- function(x, m1, sd1, m2, sd2, p1, p2){
	                  dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

pdf("TCGA_TGFB_normalmix.pdf")
plot(mixmdl,which=2)
lines(density((input_cp$tgfb/max(input_cp$tgfb))*100), lty=2, lwd=2)

x1_tgfb<-uniroot.all(f, lower=0, upper=100, m1=mixmdl$mu[1], sd1=mixmdl$sigma[1], m2=mixmdl$mu[2], sd2=mixmdl$sigma[2], p1=mixmdl$lambda[1], p2=mixmdl$lambda[2])
x2_tgfb<-uniroot.all(f, lower=0, upper=100, m1=mixmdl$mu[2], sd1=mixmdl$sigma[2], m2=mixmdl$mu[3], sd2=mixmdl$sigma[3], p1=mixmdl$lambda[2], p2=mixmdl$lambda[3])

abline(v=c(x1_tgfb,x2_tgfb),col="orange",lty=2,lwd=2)

dev.off()


set.seed(123)
mixmdl = normalmixEM((input_cp$mock/max(input_cp$mock))*100,k=3,fast=TRUE,maxit=500000)
mixmdl$mu

pdf("TCGA_MOCK_normalmix.pdf")
plot(mixmdl,which=2)
lines(density((input_cp$mock/max(input_cp$mock))*100), lty=2, lwd=2)

x1_mock<-uniroot.all(f, lower=0, upper=100, m1=mixmdl$mu[1], sd1=mixmdl$sigma[1], m2=mixmdl$mu[2], sd2=mixmdl$sigma[2], p1=mixmdl$lambda[1], p2=mixmdl$lambda[2])
x2_mock<-uniroot.all(f, lower=0, upper=100, m1=mixmdl$mu[2], sd1=mixmdl$sigma[2], m2=mixmdl$mu[3], sd2=mixmdl$sigma[3], p1=mixmdl$lambda[2], p2=mixmdl$lambda[3])

abline(v=c(x1_mock,x2_mock),col="orange",lty=2,lwd=2)

dev.off()

input_cp$mock<-(input_cp$mock/max(input_cp$mock))*100
input_cp$tgfb<-(input_cp$tgfb/max(input_cp$tgfb))*100


knn_df_tcga_mock$pseudospace<-(knn_df_tcga_mock$pseudospace/max(knn_df_tcga_mock$pseudospace))*100
knn_df_tcga_tgfb$pseudospace<-(knn_df_tcga_tgfb$pseudospace/max(knn_df_tcga_tgfb$pseudospace))*100


library(gridExtra)
library(ggforce)


pdf('KNN_projection_TCGA_to_MCF10A_mock_vs_tgfb_v2.pdf')

mock<-ggplot(knn_df_tcga_mock, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")
tgfb<-ggplot(knn_df_tcga_tgfb, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C")+ coord_flip()

sp <- ggplot(input_cp, aes(x=mock, y=tgfb,col=groups)) +geom_point(alpha = 5/10)+geom_hline(yintercept=c(x1_tgfb,x2_tgfb), linetype="dashed", color = "red")+geom_vline(xintercept=c(x1_mock,x2_mock),linetype="dashed",color="red")

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

#
# Define groups in the pseudospace using mesenchymal and epithelial markers
#
library(ggplot2)
library(gridExtra)
library(cluster)

#source("/home/data/pseudospace/tcga_nocorrect_primary_tumors/explore_markers_in_pseudotime.R")
source(paste(current_dir,"explore_markers_in_pseudotime_epithelial.R",sep="/"))

#
# Save results PAM
#

input_cp_pemt<-list_results_clusters_markers[["no_pemt"]]
input_cp_pemt$clusters2<-rep("null",nrow(input_cp_pemt))

input_cp_pemt$clusters2[which(input_cp_pemt$clusters==1)]<-"epithelial"
input_cp_pemt$clusters2[which(input_cp_pemt$clusters==2)]<-"mixed"
input_cp_pemt$clusters2[which(input_cp_pemt$clusters==3)]<-"Mesenchymal"

setwd(output_dir)
write.table(input_cp_pemt,file="proj_pseudospace_mcf_tgfb_mock_pam3_without_pemt.txt",quote=F,row.names=F,sep="\t")

#
# Save results mclust
#

source(paste(current_dir,"explore_markers_in_pseudotime_epithelial.R",sep="/"))

input_cp_pemt<-list_results_clusters_markers[["with_pemt"]]
input_cp_pemt$clusters2<-rep("null",nrow(input_cp_pemt))

input_cp_pemt$clusters2[which(input_cp_pemt$clusters==1)]<-"mixed"
input_cp_pemt$clusters2[which(input_cp_pemt$clusters==2)]<-"epithelial"
input_cp_pemt$clusters2[which(input_cp_pemt$clusters==3)]<-"Mesenchymal"

write.table(input_cp_pemt,file="proj_pseudospace_mcf_tgfb_mock_Mclust3_with_pemt.txt",quote=F,row.names=F,sep="\t")