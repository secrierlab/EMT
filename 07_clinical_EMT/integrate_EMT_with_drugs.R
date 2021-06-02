library(openxlsx)
library(tidyr)
library(glmnet)
library(reshape)
library(ggrepel)
library(RColorBrewer)
library(wesanderson)
#
# Define some functions to use 
#


scatterplot_IC50_gf<-function(input_for_plot,xmin_legend=2,ymin_legend=3,min_y=2,max_y=10,list_tissues){
  #http://vrl.cs.brown.edu/color
  all_colors<-rev(c("#4f8c9d", "#ff69b4", "#800000", "#aedbf0", "#069668", "#ffff54", "#ff00ff", "#bfcd8e", "#7c440e", "#2524f9", "#613F75", "#ff8c00", "#ee0d0e"))
  input_for_plot$tissue_type<-as.factor(input_for_plot$tissue_type)
  
  names(all_colors)<-levels(unique(input_for_plot$tissue_type))
  
  input_for_plot$label2<-gsub(paste(input_for_plot$feature_name,input_for_plot$drug_name,sep="\n"),pattern="_mut",replacement="")
  
  input_for_plot[-which(input_for_plot$feature_delta_mean_ic50 > xmin_legend | input_for_plot$feature_delta_mean_ic50 <= -xmin_legend),"label2"]<-""
  input_for_plot[which(-log10(input_for_plot$feature_pval)<ymin_legend),"label2"]<-""
  
  p<-ggplot(input_for_plot,aes(x=feature_delta_mean_ic50, y=-log10(feature_pval),color=tissue_type,label=label2))  +
     geom_point(alpha=0.4)+ylim(c(min_y,max_y))+scale_colour_manual(name = "tissue_type",values = all_colors)+geom_text_repel(
       min.segment.length = 0, seed = 42, box.padding = 0.5) + geom_vline(xintercept = c(-xmin_legend,0,xmin_legend),linetype="dotted",size=0.5)+geom_hline(yintercept = ymin_legend,linetype="dotted",size=0.5)+
     theme_bw()+labs(y="-log10(Pvalue)",x="IC50 Effect") 
  
  print(p)
  
}

getIC50fromgenomicfeatures<-function(input_anova,markers_to_consider,pval_thr=0.05){
  
  markers_for_anova_mut<-grep(markers_to_consider,pattern="mut",value=T)
  markers_for_anova_cnv<-grep(markers_to_consider,pattern="mut",value=T,invert=T)
  
  anova_mut<-input_anova[which(input_anova$feature_name%in%markers_for_anova_mut),]
  anova_mut2<-anova_mut[anova_mut$feature_pval<=pval_thr,]
  
  return(anova_mut2)
  
}

#
# Integrate EMT with drugs response
#

gdsc_table_with_emt<-read.delim(file="/home/guidantoniomt/pseudospace/gdsc/GDSC_on_pseudotime_EMT_MCF_withEMTscores.txt")
table_cl_drugs_ic50<-read.xlsx("/home/guidantoniomt/pseudospace/gdsc/GDSC2_fitted_dose_response_25Feb20.xlsx",sheet=1)
table(table_cl_drugs_ic50$TCGA_DESC)

#get the cell line in common between those one used to compute the EMT score and IC50 data
table_cl_drugs_ic50_select<-table_cl_drugs_ic50[which(table_cl_drugs_ic50$CELL_LINE%in%gdsc_table_with_emt[,1]),which(colnames(table_cl_drugs_ic50)%in%c("CELL_LINE_NAME","DRUG_NAME","LN_IC50"))]

table_cl_drugs_ic50_select2<-aggregate(.~CELL_LINE_NAME+DRUG_NAME, table_cl_drugs_ic50_select, mean)

ic50matrix<-spread(table_cl_drugs_ic50_select2, key = DRUG_NAME, value = LN_IC50)

#
#boxplot emt by tissue
#
table_temp_drug<-table_cl_drugs_ic50[which(table_cl_drugs_ic50$CELL_LINE%in%gdsc_table_with_emt[,1]),which(colnames(table_cl_drugs_ic50)%in%c("CELL_LINE_NAME","TCGA_DESC","DRUG_NAME","LN_IC50"))]
emt_table<-merge(gdsc_table_with_emt[,-2],table_temp_drug,by.x="gdsc_ID",by.y="CELL_LINE_NAME",all.x=T,all.y=F)
emt_table2<-emt_table[,c(3,2)]

require(ggpubr)

pdf("boxplot_emt_score_bytissue.pdf",width=15)

all_colors<-rev(c("#4f8c9d", "#ff69b4", "#800000", "#aedbf0", "#069668", "#ffff54", "#ff00ff", "#bfcd8e", "#7c440e", "#2524f9", "#613F75", "#ff8c00", "#ee0d0e"))
colors2<-brewer.pal(n = 8, name = "Dark2")
colors3<-brewer.pal(n = 8, name = "Paired")

ggviolin(emt_table2, "TCGA_DESC", "emt_score", fill = "TCGA_DESC",
         palette = c(all_colors,colors2,colors3),
         add = "boxplot", add.params = list(fill = "white"))+theme(axis.text.x=element_text(angle = -90, hjust = 0))

dev.off()

#
#Correlation EMT score with IC50 of the drugs
#

emt_ic50<-merge(gdsc_table_with_emt[,-2],ic50matrix,by.x="gdsc_ID",by.y="CELL_LINE_NAME")
emt_ic50[is.na(emt_ic50)]<-0

emt_ic50_lm<-emt_ic50[,-1]

list_drugs_for_cor<-colnames(emt_ic50_lm[,-1])

rho_all<-NULL
pvalue_all<-NULL

for(current_drug in list_drugs_for_cor){
  
  df_estimate<-cor.test(emt_ic50_lm[,1],emt_ic50_lm[,colnames(emt_ic50_lm)%in%current_drug],method="spearman",exact=FALSE)$estimate
  df_pvalue<-cor.test(emt_ic50_lm[,1],emt_ic50_lm[,colnames(emt_ic50_lm)%in%current_drug],method="spearman",exact=FALSE)$p.value
  
  rho_all<-c(rho_all,df_estimate)
  pvalue_all<-c(pvalue_all,df_pvalue)
  
}

df_cor_drugs<-data.frame(drugs=list_drugs_for_cor,
                         cor=rho_all,
                         pvalue=pvalue_all,
                         padjust=p.adjust(pvalue_all,"BH"),
                         padjust_log=-log(p.adjust(pvalue_all,"BH"),10))

df_cor_drugs2<-df_cor_drugs[df_cor_drugs$padjust<=0.001,]
drugs_significantly_correlated<-df_cor_drugs2[,1]

melt_df_IC50<-melt(emt_ic50_lm[,-1])
melt_df_IC50_sig_drugs<-unique(melt_df_IC50[melt_df_IC50[,1]%in%drugs_significantly_correlated,])
colnames(melt_df_IC50_sig_drugs)[2]<-"LNIC50"

library(viridis)
library(ggplot2)
library(ggridges)
library(cowplot)
library("readxl")

#
# Correlation all drugs with EMT
#

pdf("correlation_emt_with_IC50.pdf",width=15)

df_cor_drugs2<-df_cor_drugs2[order(df_cor_drugs2$cor,decreasing=F),]
df_cor_drugs2$drugs<-factor(df_cor_drugs2$drugs,levels=df_cor_drugs2$drugs)

order_drugs<-df_cor_drugs2$drugs

cdat_sp <- ggplot(df_cor_drugs2, aes(x = 1, y = drugs,size = padjust_log,color=cor)) +
  geom_point() + scale_colour_viridis(option = "plasma",direction = -1)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("Correlation")

melt_df_IC50_sig_drugs$variable<-factor(melt_df_IC50_sig_drugs$variable,level=unique(df_cor_drugs2$drugs))

ggrid<-ggplot(melt_df_IC50_sig_drugs, aes(x = LNIC50, y = variable, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_continuous(type = "viridis",name="IC50") +theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), 
                                                             axis.line = element_line(colour = "black"),
                                                             axis.title.y=element_blank(),
                                                             axis.ticks.y=element_blank())

ggrid2<-ggplot(melt_df_IC50_sig_drugs, aes(x = LNIC50, y = variable, fill = factor(stat(quantile))))+stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 2, quantile_lines = TRUE) + scale_fill_viridis_d(name = "Quartiles",option="D",direction=-1)+theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), 
                                                             axis.line = element_line(colour = "black"),
                                                             axis.title.y=element_blank(),
                                                             axis.text.y=element_blank(),
                                                             axis.ticks.y=element_blank())

print(plot_grid(cdat_sp,ggrid,width=c(0.2,0.8),nrow=1))

print(plot_grid(cdat_sp,ggrid2,width=c(0.2,0.8),nrow=1))
# 
# melt_df_IC50_sig_drugs_binary<-data.frame()
# 
# for(drug_for_plot in unique(melt_df_IC50_sig_drugs[,1])){
#   
#   temp_df<-melt_df_IC50_sig_drugs[melt_df_IC50_sig_drugs[,1]%in%drug_for_plot,]
#   temp_df$binary<-ifelse(temp_df[,2]>=mean(temp_df[,2]),"r","s")
# 
#   # temp_df2<-cbind(drug_for_plot,data.frame(table(temp_df$binary)/sum(table(temp_df$binary))))
#   temp_df2<-cbind(drug_for_plot,data.frame(table(temp_df$binary)))
#   
#   colnames(temp_df2)<-c("drugs","status","Freq")
#                            
#   melt_df_IC50_sig_drugs_binary<-rbind(melt_df_IC50_sig_drugs_binary,temp_df2)  
# }
# 
# melt_df_IC50_sig_drugs_binary$Freq<-as.numeric(melt_df_IC50_sig_drugs_binary$Freq)
# 
# p<-ggplot(melt_df_IC50_sig_drugs_binary, aes(x = drugs, y = Freq, fill = status)) +
#   geom_col(position = "fill")
# 

dev.off()

#
# Generalized linear model with Y=emt and predictors are the drugs
#

# How to binarize the IC50 values:
# https://www.nature.com/articles/s41598-019-50720-0#Sec10
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5310525/


list_drugs_for_cor<-colnames(emt_ic50_lm[,-1])

drugs_all_glm<-NULL
pvalue_all_glm<-NULL
or_all_glm<-NULL
rsquared_all_glm<-NULL
cor_all_glm<-NULL

CI_all_glm<-data.frame()

for(current_drug in list_drugs_for_cor){
  
  current_drug_df<-emt_ic50_lm[,colnames(emt_ic50_lm)%in%c("emt_score",current_drug)]
  colnames(current_drug_df)[2]<-"drug"
  current_drug_df$drug<-ifelse(current_drug_df$drug>=median(current_drug_df$drug),1,0)
  # 
  # ggplot(current_drug_df, aes(x=emt_score, y=drug)) + geom_point() + 
  #   stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)
  
  plot(current_drug_df$drug ~ current_drug_df$emt_score )
  
  if(length(unique(current_drug_df$drug))==2){
  
  drugs_all_glm<-c(drugs_all_glm,current_drug)
    
  glm_formula<-as.formula(drug~emt_score)
  
  glm_model<-glm(glm_formula,current_drug_df, family = "binomial",control = list(maxit = 100))
  
  nullmod <- glm(drug~1,current_drug_df, family="binomial",control = list(maxit = 100))
  # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
  # https://stats.stackexchange.com/questions/8511/how-to-calculate-pseudo-r2-from-rs-logistic-regression
  rsquared<-as.numeric(1-logLik(glm_model)/logLik(nullmod))
  rsquared_all_glm<-c(rsquared_all_glm,rsquared)
  
  cor_all_glm<-c(cor_all_glm,cor(current_drug_df[,1],current_drug_df[,2],method="spearman"))
  
  pvalue_current_drugs<-summary(glm_model)$coefficients[,4][2]
  pvalue_all_glm<-c(pvalue_all_glm,pvalue_current_drugs)
  
  estimation_current_drugs<-summary(glm_model)$coefficients[,1][2]
  
  or<-exp(coef(glm_model)[2])
  confint<-exp(confint(glm_model)[2,])
  
  or_all_glm<-c(or_all_glm,or)
  CI_all_glm<-rbind(CI_all_glm,confint)
  
  }
  
}
colnames(CI_all_glm)<-c("CI_lower","CI_higher")

CI_all_glm_log<-log(CI_all_glm,2)
colnames(CI_all_glm_log)<-c("CI_lower_log","CI_higher_log")

results_glm<-data.frame(drugs=drugs_all_glm,
                        pvalue=pvalue_all_glm,
                        padjust=p.adjust(pvalue_all_glm,"BH"),
                        padjust_log=-log(p.adjust(pvalue_all_glm,"BH"),10),
                        rsquare=rsquared_all_glm,
                        cor=cor_all_glm,
                        OR=or_all_glm,
                        LOG_OR=log(or_all_glm,2),
                        CI_all_glm,
                        CI_all_glm_log)

results_glm2<-results_glm[results_glm$padjust<=0.001,]


#
# Now integrate the genomic features
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

setwd("/home/guidantoniomt/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
mut_genes_mes_vs_epi<-paste(sapply(strsplit(grep(markers_mes_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1),"_mut",sep="")
cnv_markers_mes_vs_epi<- sapply(strsplit(grep(sub(grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement=""),pattern="focal",value=T),split="_"),"[[",1)
all_markers_mes_vs_epi<-c(mut_genes_mes_vs_epi,cnv_markers_mes_vs_epi)

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
mut_genes_pemt_vs_epi<-paste(sapply(strsplit(grep(markers_pemt_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1),"_mut",sep="")
cnv_markers_pemt_vs_epi<-sapply(strsplit(grep(sub(grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement=""),pattern="focal",value=T),split="_"),"[[",1)
all_markers_pemt_vs_epi<-c(mut_genes_pemt_vs_epi,cnv_markers_pemt_vs_epi)

# Do not consider for the analysis these groups:
# "PANCANCER_ANOTA","GBM_ANOVA","LAML_ANOVA","DLBC_ANOVA","LGG_ANOVA"

list_tissues<-c("COREAD_ANOVA","BLCA_ANOVA","LIHC_ANOVA","KIRC_ANOVA","PAAD_ANOVA","OV_ANOVA","LUSC_ANOVA","THCA_ANOVA","ESCA_ANOVA","PRAD_ANOVA","BRCA_ANOVA","STAD_ANOVA","LUAD_ANOVA","SKCM_ANOVA","HNSC_ANOVA","UCEC_ANOVA")

setwd("/home/guidantoniomt/pseudospace/gdsc/")

all_tissues_ic50<-data.frame()

for(lt in 1:length(list_tissues)){
  
  print(lt)
  
  tab_temp_tissue_ic50<-read_excel("ANOVA_results_GDSC2_20Feb20.xlsx",sheet=list_tissues[lt])
    
  all_tissues_ic50<-rbind(all_tissues_ic50,tab_temp_tissue_ic50)

}

input_anova <- all_tissues_ic50
markers_to_consider_mes_vs_epi <- setdiff(all_markers_mes_vs_epi,all_markers_pemt_vs_epi)
emt_table <- gdsc_table_with_emt

tab_anova_ic50_mes_vs_epi<-getIC50fromgenomicfeatures(input_anova,markers_to_consider_mes_vs_epi,pval_thr=0.01)

pdf("IC50_mutations_mes_vs_epi.pval0.01.pdf")
scatterplot_IC50_gf(tab_anova_ic50_mes_vs_epi,xmin_legend=1,ymin_legend=2,min_y=2,max_y=4,list_tissues)
dev.off()

markers_to_consider_pEMT_vs_epi <- setdiff(all_markers_pemt_vs_epi,all_markers_mes_vs_epi)
emt_table <- gdsc_table_with_emt

tab_anova_ic50_pEMT_vs_epi<-getIC50fromgenomicfeatures(input_anova,markers_to_consider_pEMT_vs_epi,pval_thr=0.01)

pdf("IC50_mutations_pEMT_vs_epi.pval0.01.pdf")
scatterplot_IC50_gf(tab_anova_ic50_pEMT_vs_epi,xmin_legend=1,ymin_legend=2.5,min_y=2,max_y=4,list_tissues)
dev.off()

markers_to_consider_common <- intersect(all_markers_pemt_vs_epi,all_markers_mes_vs_epi)

tab_anova_ic50_common<-getIC50fromgenomicfeatures(input_anova,markers_to_consider_common,pval_thr=0.01)

pdf("IC50_mutations_common_mes_vs_epi_and_pemt_vs_epid.pval0.01.pdf")
scatterplot_IC50_gf(tab_anova_ic50_common,xmin_legend=1.5,ymin_legend=2.5,min_y=2,max_y=4,list_tissues)
dev.off()

