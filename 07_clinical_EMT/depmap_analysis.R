#
# Depmep analysis
#

# https://depmap.org/repurposing/
  
library(depmap)
library(data.table)

setwd("/home/guidantoniomt/pseudospace/CCLE")
load("CCLE_pseudotime_MCF.RData")

CCLE_pseudotime_EMT_HMM

setwd("/home/guidantoniomt/pseudospace/depmap")
annotation_drug<-read.csv(file="primary-screen-replicate-collapsed-treatment-info.csv",stringsAsFactors = F)
annotation_cellline<-read.csv(file="primary-screen-cell-line-info.csv",stringsAsFactors = F)
  
tab_input<-fread(file="primary-screen-replicate-collapsed-logfold-change.csv")
colnames(tab_input)[1]<-"code_cellline"

annotation_cellline2<-merge(annotation_cellline[,c(2:4)],tab_input,by.x="depmap_id",by.y="code_cellline")

CCLE_with_drugs<-merge(CCLE_pseudotime_EMT_metpot,annotation_cellline2,by.x="CCLE_ID",by.y="ccle_name",all.x=F,all.y=F)
CCLE_with_drugs2<-CCLE_with_drugs[,c(9,12:ncol(CCLE_with_drugs))]

list_comparisons<-vector(mode="list",3)

list_comparisons[[1]]<-c("non_metastatic","weakly_metastatic")
list_comparisons[[2]]<-c("non_metastatic","metastatic")
list_comparisons[[3]]<-c("weakly_metastatic","metastatic")

results_all_comparison<-data.frame()

for(comp_states in 1:length(list_comparisons)){
  
  conditions_available<-list_comparisons[[comp_states]]
  
  tab_drugs_current_comparison<-CCLE_with_drugs2[which(CCLE_with_drugs2[,1]%in%conditions_available),]
  
  drugs_all_lm<-NULL
  pvalue_all_lm<-NULL
  or_all_lm<-NULL
  rsquared_all_lm<-NULL
  cor_all_lm<-NULL
  CI_all_lm<-data.frame()
  mu1_all<-NULL
  mu2_all<-NULL
  fc_all<-NULL
  
  for(drug in 2:ncol(tab_drugs_current_comparison)){
    
    print(drug)
    
    current_drugs<-colnames(tab_drugs_current_comparison)[drug]
    drugs_all_lm<-c(drugs_all_lm,current_drugs)
    
    df_drugs<-tab_drugs_current_comparison[,colnames(tab_drugs_current_comparison)%in%c("status_metpot",current_drugs)]
    colnames(df_drugs)[2]<-"drug"
  
    df_drugs[,1]<-as.factor(df_drugs[,1])
    
    mu1<-mean(df_drugs[df_drugs[,1]%in%conditions_available[1],2],na.rm=T)
    mu2<-mean(df_drugs[df_drugs[,1]%in%conditions_available[2],2],na.rm=T)
    
    mu1_all<-c(mu1_all,mu1)
    mu2_all<-c(mu2_all,mu2)
    
    fc_all<-c(fc_all,mu2-mu1)
    
    formula_to_use<-as.formula(drug~status_metpot)
    
    # #the reference level is always the first one - see the lists (epi, epi, mix)
    df_drugs$status_metpot<-relevel(df_drugs$status_metpot, ref=conditions_available[1])
    
    lmmodel<-lm(formula_to_use,df_drugs)
    
    pvalue_current_drugss<-summary(lmmodel)$coefficients[2,4]
    
    pvalue_all_lm<-c(pvalue_all_lm,pvalue_current_drugss)
    
    or<-exp(coef(lmmodel)[2])
    # confint<-exp(confint(lmmodel)["protein",])
    
    or_all_lm<-c(or_all_lm,or)
    # CI_all_lm<-rbind(CI_all_lm,confint)
    
  }
  
  results_current_comparison<-data.frame(comparison=paste(conditions_available[1],conditions_available[2],sep="_"),
                                         drugs=drugs_all_lm,
                                         pvalue=pvalue_all_lm,
                                         padjust=p.adjust(pvalue_all_lm,"BH"),
                                         padjust_log=-log(p.adjust(pvalue_all_lm,"BH"),10),
                                         or=or_all_lm,
                                         logOR=log(or_all_lm),
                                         mu1=mu1_all,
                                         mu2=mu2_all,
                                         fc=fc_all)
  
  results_all_comparison<-rbind(results_all_comparison,results_current_comparison)
}

setwd("/home/guidantoniomt/pseudospace/depmap")

results_all_comparison$status<-rep("no_sig",nrow(results_all_comparison))
results_all_comparison$status[results_all_comparison$padjust<=0.05 & results_all_comparison$logOR>=0.25]<-"up"
results_all_comparison$status[results_all_comparison$padjust<=0.05 & results_all_comparison$logOR<=-0.25]<-"down"

write.table(results_all_comparison,file="depmap_lm_drugs_metpotential_padjust_0.05_logOR_0.5.txt",sep="\t",row.names=F,quote=F)


list_resDPA<-split(results_all_comparison,results_all_comparison$comparison)

library(ggpubr)
library(clusterProfiler)
library(ggrepel)

pdf("depmap_lm_drugs_metpotential_padjust_0.05_logOR_0.5.volcano.pdf")

for(lr in 1:length(list_resDPA)){
  
  current_comparison<-unique(list_resDPA[[lr]][,1])
  
  current_DPA<-list_resDPA[[lr]]
  
  current_DPA<-merge(current_DPA,annotation_drug,by.x="drugs",by.y="column_name",all.x=F,all.y=F)
    
  #do also a simple Volcano Plot
  p2<-ggplot(current_DPA, aes(x = log(or), y = padjust_log)) +
    geom_point(aes(color = status)) +
    scale_color_manual(values = c("blue2", "grey","red2")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel( data = subset(current_DPA, abs(log(or))>0.25),
                     aes(label = name),
                     size = 5,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines")
    )+labs(title=unique(list_resDPA[[lr]][,1]))
  
  print(p2)

}

dev.off()
