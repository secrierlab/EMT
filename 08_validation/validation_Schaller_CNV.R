library(openxlsx)
library(biomaRt)
library(patchwork)
library(ggplot2)

#
# Upload the data from Swiatkowska
#

setwd("/home/guidantoniomt/pseudospace/Schaller_mouse_validation")
tab_hits<-read.xlsx("mmc4.xlsx",sheet=3)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV1 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = as.character(tab_hits[,1]), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

tab_hits2<-merge(genesV1,tab_hits,by.x="MGI.symbol",by.y="SYMBOL")

tab_hits_for_split<-tab_hits2[,-c(1:4)]
tab_hits_genes_ann<-tab_hits2[,c(1:4)]

tab_hit_list<-split.default(tab_hits_for_split, rep(1:48, each = 3))

#
# Upload the biomarkers
#

setwd("/home/guidantoniomt/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

getFeaturesSS<-function(res_lasso,ss=800){
  
  df_coef2<-do.call(rbind,res_lasso)
  
  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]
  
  features<-data.frame(table(df_coef2[,1]))
  
  features_freq_df<-features[order(features[,2],decreasing=T),]
  
  features_to_select<-features_freq_df[features_freq_df[,2]>=ss,1]
  
  return(features_to_select)
}


getCoefSS<-function(res_lasso,features_to_select){
  
  df_coef2<-do.call(rbind,res_lasso)
  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]
  
  features<-data.frame(table(df_coef2[,1]))
  
  input_for_boxplot<-df_coef2[which(df_coef2[,1]%in%features_to_select),]
  colnames(input_for_boxplot)[2]<-"values"
  
  return(input_for_boxplot)
  
}

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
cnv_markers_mes_vs_epi<-gsub(grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
cnv_markers_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T)
coef_mes_vs_epi<-getCoefSS(res_lasso,c(cnv_markers_mes_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_mes_vs_epi2<-merge(coef_mes_vs_epi,mu,by="coef")
coef_mes_vs_epi2[,1]<-gsub(coef_mes_vs_epi2[,1],pattern="^X",replacement="")

cnv_markers_mes_vs_epi<-intersect(cnv_markers_mes_vs_epi,unique(coef_mes_vs_epi2[coef_mes_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
cnv_markers_pemt_vs_epi<-gsub(grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
cnv_markers_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T)
coef_pemt_vs_epi<-getCoefSS(res_lasso,c(cnv_markers_pemt_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_pemt_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_pemt_vs_epi2<-merge(coef_pemt_vs_epi,mu,by="coef")
coef_coef_pemt_vs_epi2[,1]<-gsub(coef_coef_pemt_vs_epi2[,1],pattern="^X",replacement="")

cnv_markers_pemt_vs_epi<-intersect(cnv_markers_pemt_vs_epi,unique(coef_coef_pemt_vs_epi2[coef_coef_pemt_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_mix<-getFeaturesSS(res_lasso,ss=800)
cnv_markers_mes_vs_mix<-gsub(grep(markers_mes_vs_mix,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
cnv_markers_mes_vs_mix_raw<-grep(markers_mes_vs_mix,pattern="dndscv",value=T,invert=T)
coef_mes_vs_mix<-getCoefSS(res_lasso,c(cnv_markers_mes_vs_mix_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_mix, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_mes_vs_mix2<-merge(coef_mes_vs_mix,mu,by="coef")
coef_coef_mes_vs_mix2[,1]<-gsub(coef_coef_mes_vs_mix2[,1],pattern="^X",replacement="")

cnv_markers_mes_vs_mix<-intersect(cnv_markers_mes_vs_mix,unique(coef_coef_mes_vs_mix2[coef_coef_mes_vs_mix2$color%in%"positive",1]))

raw_markers_specifics_mes_vs_epi<-setdiff(c(cnv_markers_mes_vs_epi),c(cnv_markers_pemt_vs_epi,cnv_markers_mes_vs_mix))
raw_markers_specifics_pemt_vs_epi<-setdiff(y=c(cnv_markers_mes_vs_epi,cnv_markers_mes_vs_mix),x=c(cnv_markers_pemt_vs_epi))
raw_markers_common<-intersect(c(cnv_markers_mes_vs_epi),c(cnv_markers_pemt_vs_epi))

raw_markers_specifics_mes_vs_mix<-setdiff(x=c(cnv_markers_mes_vs_mix),y=c(cnv_markers_mes_vs_epi,cnv_markers_pemt_vs_epi))

list_common_all<-list(c(cnv_markers_mes_vs_epi),c(cnv_markers_mes_vs_mix),c(cnv_markers_pemt_vs_epi))

common_between_three_groups<-Reduce(intersect,list_common_all)

#
# See the effect of the CNA in CCLE and on the metastatic potential
#

list_features<-vector(mode="list",6)
list_features[[1]]<-raw_markers_specifics_mes_vs_epi
list_features[[2]]<-raw_markers_specifics_pemt_vs_epi
list_features[[3]]<-raw_markers_specifics_mes_vs_mix
list_features[[4]]<-raw_markers_common
list_features[[5]]<-c("NONO","MPL","KIAA1459","ETNK1","FIP1L1","CNTNAP2","CREB3L2")
list_features[[6]]<-c("PRDM2")

names(list_features)<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common","luad","breast")

df_features<-data.frame(do.call(c,list_features))
df_features2<-data.frame(comp=rownames(df_features),df_features)

df_features2[,1]<-gsub(df_features2[,1],pattern="[0-9]",replacement="")
colnames(df_features2)[2]<-"symbols"
df_features2$symbols<-gsub(df_features2$symbols,pattern="_focal",replacement="")
df_features2<-df_features2[grep(df_features2$symbols,pattern="arm",invert=T),]

ALL_TF_RES<-data.frame()

for(thl in 1:length(tab_hit_list)){
  
  current_experiment<-cbind(tab_hits_genes_ann,tab_hit_list[[thl]])
  int_biomarkers_hits<-merge(df_features2,current_experiment,by.x="symbols",by.y="HGNC.symbol")
  int_biomarkers_hits_sig<-int_biomarkers_hits[which(abs(int_biomarkers_hits$log2FoldChange)>=1 & int_biomarkers_hits$padj<=0.05),]
  int_biomarkers_hits_sig2<-data.frame(TF=unique(int_biomarkers_hits_sig[,6]),int_biomarkers_hits_sig)
  
  int_biomarkers_hits_sig2<-int_biomarkers_hits_sig2[,c(1,2,3,8)]
  colnames(int_biomarkers_hits_sig2)[4]<-"log2FoldChange"
  
  
  ALL_TF_RES<-rbind(ALL_TF_RES,int_biomarkers_hits_sig2)
}

ALL_TF_RES2<-ALL_TF_RES[-which(ALL_TF_RES$TF%in%c("0dCtr","4dCtr")),]
ALL_TF_RES2$TF<-toupper(ALL_TF_RES2$TF)

setwd("/home/guidantoniomt/pseudospace/Schaller_mouse_validation")

pdf("Schaller_mouse_validation_CNV.pdf")

ALL_TF_RES2_ml<-ALL_TF_RES2[-which(ALL_TF_RES2$comp%in%c("luad","breast")),]

ALL_TF_RES2_ml$comp<-factor(ALL_TF_RES2_ml$comp,levels=c("hemt_vs_epi","common"))

ph<-ggplot(ALL_TF_RES2_ml, aes(TF, symbols,fill=log2FoldChange)) +
  geom_tile()+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown2", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="bottom")

ph2<-ggplot(ALL_TF_RES2_ml, aes(comp, symbols,fill=comp)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="bottom")

print(ph+ph2+plot_layout(widths = c(0.8, 0.1)))

ALL_TF_RES2_ml<-ALL_TF_RES2[-which(ALL_TF_RES2$comp%in%c("luad","breast")),]

#
# Tissue specifics markers
#

ALL_TF_RES2_ts<-ALL_TF_RES2[which(ALL_TF_RES2$comp%in%c("luad","breast")),]

ph<-ggplot(ALL_TF_RES2_ts, aes(TF, symbols,fill=log2FoldChange)) +
  geom_tile()+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown2", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="bottom")

ph2<-ggplot(ALL_TF_RES2_ts, aes(comp, symbols,fill=comp)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="bottom")

print(ph+ph2+plot_layout(widths = c(0.8, 0.1)))

dev.off()

