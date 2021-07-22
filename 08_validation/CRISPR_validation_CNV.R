library(data.table)
library(openxlsx)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(vcd)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

RadarTheme<-theme(
  panel.background = element_rect(fill = "white", colour = "white", linetype = "black"),
  panel.grid.major = element_line(size = 0.25, linetype = 'black', colour = "white"), 
  panel.grid.minor = element_line(size = 0, linetype = 'black',colour = "black"))


setwd("/home/data/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

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

raw_markers_specifics_mes_vs_epi<-setdiff(c(cnv_markers_mes_vs_epi),c(cnv_markers_pemt_vs_epi))
raw_markers_specifics_pemt_vs_epi<-setdiff(y=c(cnv_markers_mes_vs_epi),x=c(cnv_markers_pemt_vs_epi))
raw_markers_common<-intersect(c(cnv_markers_mes_vs_epi),c(cnv_markers_pemt_vs_epi))

raw_markers_specifics_mes_vs_mix<-setdiff(x=c(cnv_markers_mes_vs_mix),y=c(cnv_markers_mes_vs_epi,cnv_markers_pemt_vs_epi))

list_common_all<-list(c(cnv_markers_mes_vs_epi),c(cnv_markers_mes_vs_mix),c(cnv_markers_pemt_vs_epi))

common_between_three_groups<-Reduce(intersect,list_common_all)

#
# See the effect of the CNA in CCLE and on the metastatic potential
#

list_features<-vector(mode="list",4)
list_features[[1]]<-raw_markers_specifics_mes_vs_epi
list_features[[2]]<-raw_markers_specifics_pemt_vs_epi
list_features[[3]]<-raw_markers_specifics_mes_vs_mix
list_features[[4]]<-raw_markers_common
#list_features[[5]]<-c("NONO_focal","MPL_focal","KIAA1459_focal","ETNK1_focal","FIP1L1_focal","CNTNAP2_focal","CREB3L2_focal")
#list_features[[6]]<-c("PRDM2_focal")

list_features[[5]]<-c("RUNX1_focal","NCOR1_focal","EPHA7_focal","TSHR_focal","CUX1_focal")
list_features[[6]]<-c("MDS2_focal","PAX7_focal","PRDM2_focal","ID3_focal","ARID1A_focal","BAP1_focal","PBRM1_focal")

names(list_features)<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common","luad","brca")

#
# Upload the CRISPR data
#

setwd("/home/data/pseudospace/CCLE")

tab_sampleinfo<-fread("/home/data/pseudospace/depmap/primary-screen-cell-line-info.csv",data.table=F)

# (DepMap 21Q2 Public+Score, CERES)

tab_crispr<-fread("CRISPR_gene_effect.csv",fill=T,data.table=F)
tab_crispr2<-merge(tab_sampleinfo,tab_crispr,by.x="depmap_id",by.y="DepMap_ID")
colnames(tab_crispr2)<-sapply(strsplit(colnames(tab_crispr2),split=" "),"[[",1)

#
# start the analysis with CRISPR data
#
tab_cripr_all_comparison<-data.frame()

for(i in 1:length(list_features)){
  
  names_comparison<-names(list_features)[i]
  
  genes_to_use<-gsub(grep(list_features[[i]],pattern="_focal",value=T),pattern="_focal",replacement="")
  
  if(length(genes_to_use)!=0){
    
  tab_crispr_current_gene<-tab_crispr2[,which(colnames(tab_crispr2)%in%c("depmap_id","primary_tissue","ccle_name",genes_to_use))]
  
  tab_crispr_current_gene_melt<-data.frame(ID=names_comparison,melt(tab_crispr_current_gene))
  
  tab_cripr_all_comparison<-rbind(tab_cripr_all_comparison,tab_crispr_current_gene_melt)
  }
  
}


groups_to_consider<-vector(mode="list",2)
groups_to_consider[[1]]<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common_all")
groups_to_consider[[2]]<-c("luad","brca")
names(groups_to_consider)<-c("ML","LUAD_BREAST")

for(gtc in 1:length(groups_to_consider)){
  
  groups_to_use<-groups_to_consider[[gtc]]
  output_string<-names(groups_to_consider)[[gtc]]
  
  custom_quantile<-function(X){quantile(X)[2]}
  
  tab_crips_subgroup<-tab_cripr_all_comparison[which(tab_cripr_all_comparison$ID %in% groups_to_use),]
  
  tab_quant<-aggregate(value ~ primary_tissue + variable, data=tab_crips_subgroup, custom_quantile)
  
  median_tissue<-aggregate(value ~ primary_tissue, data=tab_crips_subgroup, median)
  median_tissue_ordered<-median_tissue[order(median_tissue$value),]
  
  median_genes<-aggregate(value ~ variable, data=tab_crips_subgroup, median)
  median_genes_ordered<-median_genes[order(median_genes$value),]
  
  tab_crips_subgroup$primary_tissue<-factor(tab_crips_subgroup$primary_tissue,levels=median_tissue_ordered[,1])
  tab_crips_subgroup$variable<-factor(tab_crips_subgroup$variable,levels=median_genes_ordered[,1])
  
  ptissue<-ggplot(tab_crips_subgroup, aes(x=primary_tissue,variable, y=value)) + 
    geom_boxplot()+ylab("CERES")+ theme_bw()
  
  tab_quant$primary_tissue<-factor(tab_quant$primary_tissue,levels=median_tissue_ordered[,1])
  tab_quant$variable<-factor(tab_quant$variable,levels=median_genes_ordered[,1])
  
  ph<-ggplot(tab_quant, aes(primary_tissue, variable,fill=value)) +
    geom_tile()+
    scale_fill_gradient2(low = "blue2", mid = "white", high = "brown2", midpoint = 0)+theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="bottom")
  
  pgenes<-ggplot(tab_crips_subgroup, aes(x=variable, y=value)) + 
    geom_boxplot()+coord_flip()+ylab("CERES")+ theme_bw()
  
  p_annotation<-ggplot(tab_crips_subgroup, aes(1, variable))+geom_tile(aes(fill=ID), colour="white")+
    coord_equal()+
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.6))
  
  
  pdf(paste(output_string,"_comparisons_CRISPR_gene_effect.CNV.July.pdf",sep=""),width=10)
  
  # print((ptissue/ph)+plot_layout(heights=c(0.20,0.80)))
  print(ph|p_annotation|pgenes)+plot_layout(widths=c(0.80,0.10,0.10))
  print(pgenes)
  print(ptissue)
  
  dev.off()
  
}



