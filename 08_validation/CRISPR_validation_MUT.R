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
mut_genes_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
mut_genes_mes_vs_epi<-sapply(strsplit(grep(markers_mes_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)
coef_mes_vs_epi<-getCoefSS(res_lasso,c(mut_genes_mes_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_mes_vs_epi2<-merge(coef_mes_vs_epi,mu,by="coef")
coef_mes_vs_epi2[,1]<-gsub(coef_mes_vs_epi2[,1],pattern="^X",replacement="")

mut_genes_mes_vs_epi_raw<-intersect(mut_genes_mes_vs_epi_raw,unique(coef_mes_vs_epi2[coef_mes_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
mut_genes_pemt_vs_epi<-sapply(strsplit(grep(markers_pemt_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)
coef_pemt_vs_epi<-getCoefSS(res_lasso,c(mut_genes_pemt_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_pemt_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_pemt_vs_epi2<-merge(coef_pemt_vs_epi,mu,by="coef")
coef_coef_pemt_vs_epi2[,1]<-gsub(coef_coef_pemt_vs_epi2[,1],pattern="^X",replacement="")

mut_genes_pemt_vs_epi_raw<-intersect(mut_genes_pemt_vs_epi_raw,unique(coef_coef_pemt_vs_epi2[coef_coef_pemt_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_mix<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_mix_raw<-grep(markers_mes_vs_mix,pattern="dndscv",value=T)
mut_genes_mes_vs_mix<-sapply(strsplit(grep(markers_mes_vs_mix,pattern="dndscv",value=T),split="_"),"[[",1)
coef_mes_vs_mix<-getCoefSS(res_lasso,c(mut_genes_mes_vs_mix_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_mix, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_mes_vs_mix2<-merge(coef_mes_vs_mix,mu,by="coef")
coef_coef_mes_vs_mix2[,1]<-gsub(coef_coef_mes_vs_mix2[,1],pattern="^X",replacement="")

mut_genes_mes_vs_mix_raw<-intersect(mut_genes_mes_vs_mix_raw,unique(coef_coef_mes_vs_mix2[coef_coef_mes_vs_mix2$color%in%"positive",1]))


raw_markers_specifics_mes_vs_epi<-setdiff(c(mut_genes_mes_vs_epi_raw),c(mut_genes_pemt_vs_epi_raw,mut_genes_mes_vs_mix_raw))
raw_markers_specifics_pemt_vs_epi<-setdiff(y=c(mut_genes_mes_vs_epi_raw),x=c(mut_genes_pemt_vs_epi_raw,mut_genes_mes_vs_mix_raw))
raw_markers_common<-intersect(c(mut_genes_mes_vs_epi_raw),c(mut_genes_pemt_vs_epi_raw))

raw_markers_specifics_mes_vs_mix<-setdiff(x=c(mut_genes_mes_vs_mix_raw),y=c(mut_genes_mes_vs_epi_raw,mut_genes_pemt_vs_epi_raw))

list_common_all<-list(c(mut_genes_mes_vs_epi_raw),c(mut_genes_mes_vs_mix_raw),c(mut_genes_pemt_vs_epi_raw))

common_between_three_groups<-Reduce(intersect,list_common_all)

list_features<-vector(mode="list",7)
list_features[[1]]<-raw_markers_specifics_mes_vs_epi
list_features[[2]]<-raw_markers_specifics_pemt_vs_epi
list_features[[3]]<-raw_markers_specifics_mes_vs_mix
list_features[[4]]<-raw_markers_common
list_features[[5]]<-common_between_three_groups
list_features[[6]]<-c("EGFR_dndscv","SETD2_dndscv","RBM10_dndscv","CDKN2A_dndscv","RB1_dndscv","MGA_dndscv","ZIC1_dndscv","NF1_dndscv","REG3A_dndscv","ARID2_dndscv","ZFP36L1_dndscv","ARID1A_dndscv")
list_features[[7]]<-c("CBFB_dndscv","NCOR1_dndscv","AKT1_dndscv","SHISA4_dndscv","NF1_dndscv","FOXA1_dndscv")

names(list_features)<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common","common_all","luad","breast")

#
# Upload the CRISPR data
#

setwd("/home/guidantoniomt/pseudospace/CCLE")

tab_sampleinfo<-fread("/home/guidantoniomt/pseudospace/depmap/primary-screen-cell-line-info.csv",data.table=F)

# (DepMap 21Q2 Public+Score, CERES)

tab_crispr<-fread("CRISPR_gene_effect.csv",fill=T,data.table=F)
tab_crispr2<-merge(tab_sampleinfo,tab_crispr,by.x="depmap_id",by.y="DepMap_ID")
colnames(tab_crispr2)<-sapply(strsplit(colnames(tab_crispr2),split=" "),"[[",1)

tab_mut_ccle<-fread(file="/home/guidantoniomt/pseudospace/CCLE/CCLE_mutations.csv",header=T,fill = TRUE,data.table=F)
tab_mut_ccle<-tab_mut_ccle[which(tab_mut_ccle$Variant_annotation%in%"damaging"),]
colnames(tab_mut_ccle)[16]<-"Samples"

#
# start the analysis with CRISPR data
#
tab_cripr_all_comparison<-data.frame()
tab_mut_all_comparison<-data.frame()

for(i in 1:length(list_features)){
  
  names_comparison<-names(list_features)[i]
    
  genes_to_use<-gsub(list_features[[i]],pattern="_dndscv",replacement="")

  tab_crispr_current_gene<-tab_crispr2[,which(colnames(tab_crispr2)%in%c("depmap_id","primary_tissue","ccle_name",genes_to_use))]
    
  tab_crispr_current_gene_melt<-data.frame(ID=names_comparison,melt(tab_crispr_current_gene))
  
  tab_cripr_all_comparison<-rbind(tab_cripr_all_comparison,tab_crispr_current_gene_melt)

  tab_mut_genes<-tab_mut_ccle[which(tab_mut_ccle[,1]%in%genes_to_use),]
    
  tab_mut_genes_melt<-data.frame(ID=names_comparison,melt(t(table(tab_mut_genes$Hugo_Symbol,tab_mut_genes$Samples))))
  colnames(tab_mut_genes_melt)<-c("ID","depmap_ID","genes","mut_status")
  
  tab_mut_all_comparison<-rbind(tab_mut_all_comparison,tab_mut_genes_melt)
  
}


groups_to_consider<-vector(mode="list",2)
groups_to_consider[[1]]<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common_all")
groups_to_consider[[2]]<-c("luad","breast")
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
  
  
  pdf(paste(output_string,"_comparisons_CRISPR_gene_effect.MUT.pdf",sep=""),width=10)
  
  # print((ptissue/ph)+plot_layout(heights=c(0.20,0.80)))
  print(ph|p_annotation|pgenes)+plot_layout(widths=c(0.80,0.10,0.10))
  print(ptissue)
  print(pgenes)
  
  dev.off()
  
  #
  # Check if there are differences in the viability of the cell lines 
  #
  
  input_tot<-merge(tab_crips_subgroup,tab_mut_all_comparison,by.x="depmap_id",by.y="depmap_ID")
  input_tot$mut_status<-ifelse(input_tot$mut_status>=1,1,0)
  
  freq_alterations<-table(input_tot$genes,input_tot$mut_status)
  
  genes_to_use<-rownames(freq_alterations[freq_alterations[,1]>=10 & freq_alterations[,2]>=10,])
  
  input_tot2<-input_tot[which(input_tot$genes%in%genes_to_use),]
            
  stat.test_pancancer <- input_tot2 %>%
    group_by(genes) %>%
    wilcox_test(value ~ mut_status)
  
  stat.test_pancancer<-data.frame(stat.test_pancancer[stat.test_pancancer$p<=0.05,])
  
  all_tissues<-unique(input_tot2$primary_tissue)
  
  stats_tissues<-data.frame()
  
  for(at in all_tissues){
  
  input_tot_sub_tissue<-input_tot2[which(input_tot$primary_tissue%in%at),]
  freq_alterations_tissue<-as.data.frame.matrix(table(input_tot_sub_tissue$genes,input_tot_sub_tissue$mut_status))
  freq_alterations_tissue_df<-data.frame(genes=rownames(freq_alterations_tissue),freq_alterations_tissue)
    
  genes_to_use_tissue<-freq_alterations_tissue_df[which(freq_alterations_tissue_df[,2]>=10 & freq_alterations_tissue_df[,3]>=10),1]
  
  input_tot_sub_tissue<-input_tot_sub_tissue[which(input_tot_sub_tissue$genes%in%genes_to_use_tissue),]
  
  stat.single.tissue <- input_tot_sub_tissue %>%
    group_by(genes) %>%
    wilcox_test(value ~ mut_status)
  
  stat.single.tissue<-data.frame(stat.single.tissue)
  
  if(nrow(stat.single.tissue[stat.single.tissue$p<=0.05,])>=1){
    
  stat_res_tissue<-data.frame(tissue=at,stat.single.tissue[stat.single.tissue$p<=0.05,])
  
  stats_tissues<-rbind(stats_tissues,stat_res_tissue)
  }
  
  }
  
  if(nrow(stats_tissues)!=0){
  
  input_boxplot<-data.frame()
  
  for(stat in 1:nrow(stats_tissues)){
  
  tissue_stat<-stats_tissues[stat,1]
  genes_stat<-stats_tissues[stat,2]
    
  part_input_tot<-input_tot2[which(input_tot2$primary_tissue%in%tissue_stat & input_tot2$genes %in%genes_stat),]
  
  input_boxplot<-rbind(input_boxplot,part_input_tot)
  
  }
    
  input_boxplot$mut_status<-as.factor(input_boxplot$mut_status)
  
  p<-ggplot(input_boxplot, aes(x=genes, y=value, fill=mut_status)) +
    geom_boxplot()+ facet_wrap(~ primary_tissue)+ylab("CERES")
  
  }
  
}



