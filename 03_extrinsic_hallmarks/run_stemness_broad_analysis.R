library(readxl)
library(data.table)
library(GSVA)
library(corrplot)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

#
# Upload the TCGA data
#
setwd(input_dir)
load("TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")
tcga<-TCGA_GEXP_ALL
tcga[,-1]<-log(tcga[,-1]+1,2)
tcga<-setDT(tcga)

tcga2 <- data.frame(tcga[, lapply(.SD, mean), by = gene_symbol])
rownames(tcga2)<-tcga2[,1]

#
# upload the catalogue with the stemness markers
#
stemness_catalog<-as.data.frame(read_excel("SuppStemnessCatalog.xlsx",col_names=F))
stemness_catalog_list <- split(stemness_catalog[,-c(1,2)], seq(nrow(stemness_catalog)))
stemness_catalog_list<-lapply(stemness_catalog_list,FUN=function(X){as.character(na.omit(as.character(X)))})
names(stemness_catalog_list)<-stemness_catalog[,2]

#
# upload segmentation data
#
knn_df_tcga<-read.delim(file="HMM_results_nstates_tumors_for_states3.withEMT.txt",stringsAsFactors=F)
knn_df_tcga$samples2<-knn_df_tcga$samples
knn_df_tcga$samples<- unlist(lapply(strsplit(as.character(knn_df_tcga$samples),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")}))

knn_df_tcga$states<-gsub(knn_df_tcga$states,pattern=1,replacement="pEMT")
knn_df_tcga$states<-gsub(knn_df_tcga$states,pattern=2,replacement="epi")
knn_df_tcga$states<-gsub(knn_df_tcga$states,pattern=3,replacement="mes")

knn_df_tcga$states<-as.factor(knn_df_tcga$states)
knn_df_tcga$states<- factor(knn_df_tcga$states, levels = c("epi", "pEMT", "mes"))

setwd(output_dir)

list_cancers<-unique(sapply(strsplit(colnames(tcga2)[-1],split="\\."),"[[",1))
stemness_list_results<-vector(mode="list",length=length(list_cancers))
names(stemness_list_results)<-list_cancers

for(i in 1:length(list_cancers)){

  tum=list_cancers[i]

  idx_tum<-grep(unlist(sapply(strsplit(colnames(tcga2),split="\\."),"[[",1)),pattern=tum)

  matrix_tum<-tcga2[,idx_tum]
  matrix_tum$var<-apply(matrix_tum,1,var)
  matrix_tum2<-matrix_tum[order(matrix_tum$var,decreasing=T),-c(ncol(matrix_tum))][1:round(nrow(matrix_tum)*30/100),]

  res<-gsva(as.matrix(matrix_tum2),stemness_catalog_list,method="ssgsea",parallel.sz=60)

  stemness_list_results[[i]]<-res
}

setwd(output_dir)
save(stemness_list_results,file="PancancerStemnessMarkers.ssgsea.RData")

input_stemness<-data.frame(matrix(0,nrow=1,ncol=length(names(stemness_catalog_list))))
colnames(input_stemness)<-names(stemness_catalog_list)

for(i in 1:length(stemness_list_results)){

  sub_stemness<-t(stemness_list_results[[i]])
  input_stemness<-rbind(input_stemness,sub_stemness)

}

input_stemness<-data.frame(samplesID=as.character(unlist(lapply(stemness_list_results,colnames))),input_stemness[-1,])
input_stemness[,1]<-unlist(lapply(strsplit(as.character(input_stemness$samplesID),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")}))

matrix_pseudotime_stemness<-merge(knn_df_tcga,input_stemness,by.x="samples",by.y="samplesID")

library(ggpubr)

pdf("stemness_scores_HMMstates.pdf")

for(i in 6:ncol(matrix_pseudotime_stemness)){

print(i)

sub_matrix_stemness<-matrix_pseudotime_stemness[,c(2,4,i)]

bp2<-ggviolin(sub_matrix_stemness,
              x="states",
              y=colnames(sub_matrix_stemness)[3],
              add="boxplot",
              add.params = list(fill = "white"),
              font.label="bold",
              fill="states",palette=c(pEMT="orange2",epi="steelblue2",mes="red2"))+stat_compare_means(comparisons=list(c('epi','mes'),c('epi','pEMT'),c('pEMT','mes')),method='wilcox.test')

print(bp2)

}

dev.off()

#
# Explore more the relations between the stemness scores and HMM states
#

# add the data from Malta et al., with a pre-computed stemness score
setwd(input_dir)

stemness_tab<-read.delim(file="/home/guidantoniomt/datasets/TCGA/TCGA_supp/stemness_scores/StemnessScores_RNAexp.txt")
stemness_tab$TCGAlong.id<-unlist(lapply(strsplit(as.character(stemness_tab$TCGAlong.id),split="\\-"),FUN=function(X){paste(X[1:4],collapse="-")}))
colnames(stemness_tab)[1]<-"samples"
stemness_tab[,1]<-as.factor(stemness_tab[,1])
stemness_tab<-data.frame(stemness_tab %>% dplyr::group_by(samples) %>% dplyr::summarize(mRNAsi=mean(mRNAsi)))

stemness_full<-merge(matrix_pseudotime_stemness,stemness_tab,by.x="samples",by.y="samples")

colnames(stemness_full)[6]<-"Miranda_et_al_PNAS"

setwd(output_dir)

pdf("comparisons_stemness_scores_for_states.pdf",height=10,width=10)

for(i in c("epi","pEMT","mes")){

print(i)

sub_matrix_stemness<-stemness_full[stemness_full$states%in%i,c(6:ncol(stemness_full))]

corrplot(cor(as.matrix(sub_matrix_stemness)), method="color",col=colorRampPalette(c("blue","white","red"))(200))

stemness_melt_state<-melt(sub_matrix_stemness)
colnames(stemness_melt_state)<-c("stemness_scores","values")

bp2<-ggviolin(stemness_melt_state,
              x="stemness_scores",
              y="values",
              add="boxplot",
              add.params = list(fill = "white"),
              font.label="bold",
              fill="stemness_scores",ylim=c(-0.5,1))+rotate_x_text(angle = 90)+stat_compare_means()
print(bp2+ggtitle(i))

#p<-stemness_melt_state %>%
#  ggplot(aes(x=reorder(stemness_scores,values), y=values)) +
#  geom_boxplot() + coord_cartesian(ylim = c(-0.5, 1))+
#  labs(y="Stemness scores", x="Signatures")+ggtitle(i)

#print(p)

}

dev.off()

#
# Create an heatmap with the median of the stemness scores for each HMM states
#

list_stemness_hmm<-split(stemness_full,stemness_full$states)

list_stemness_hmm_tissues<-lapply(list_stemness_hmm,FUN=function(X){split(X,X$tumors)})

aggregate_stemness_scores_hmm<-lapply(list_stemness_hmm_tissues,FUN=function(X){lapply(X,FUN=function(X){apply(X[,-c(1:5)],2,quantile,0.75)})})

epi_raw<-do.call(rbind,aggregate_stemness_scores_hmm[["epi"]])
epi_clean<-melt(data.frame(status="epi",tissues=rownames(epi_raw),epi_raw))

pemt_raw<-do.call(rbind,aggregate_stemness_scores_hmm[["pEMT"]])
pemt_clean<-melt(data.frame(status="pEMT",tissues=rownames(pemt_raw),pemt_raw))

mes_raw<-do.call(rbind,aggregate_stemness_scores_hmm[["mes"]])
mes_clean<-melt(data.frame(status="mes",tissues=rownames(mes_raw),mes_raw))

full_matrix<-rbind(epi_clean,pemt_clean,mes_clean)

pdf("heatmap_stemness_for_states.pdf",width=16,height=4)
ggplot(full_matrix, aes(tissues, variable, fill= value)) +
        geom_tile(color = "black")+scale_fill_gradient2(low = "white", mid = "orange2", high = "blue2", midpoint = .25)+facet_grid(tissues~status)+theme_bw()
dev.off()

pdf("heatmap_stemness_for_states.pdf",width=16,height=4)
ggplot(full_matrix, aes(tissues, variable, fill= value)) +
        geom_tile(color = "black")+facet_wrap(.~status)+scale_fill_gradient2(low = "white", mid = "orange2", high = "blue2", midpoint = .25)+theme_bw()
dev.off()

