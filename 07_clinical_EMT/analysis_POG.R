library(data.table)
library(openxlsx)
library(biomaRt)
library(ggpmisc)
library(ggridges)
library(viridis)
library(cowplot)
library(ggpubr)
library(dplyr)
library(plyr)
library(survival)
library(survminer)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

plotSignificantDrugsGenes<-function(matrix_drugs_genes_pvalue_aov,input_nested_rid,class1,class2,therapy_table){
  
  for(i in 1:nrow(matrix_drugs_genes_pvalue_aov)){
    
    drug<-rownames(matrix_drugs_genes_pvalue_aov[i,])
    
    genes_current_drug<-matrix_drugs_genes_pvalue_aov[i,]
    genes_current_drug<-as.numeric(genes_current_drug)
    names(genes_current_drug)<-colnames(matrix_drugs_genes_pvalue_aov[i,])
    
    raw_genes_current_drug_sig<-genes_current_drug[which(genes_current_drug<=0.05)]
    genes_current_drug_sig<-names(genes_current_drug[which(genes_current_drug<=0.05)])
    
    if(length(genes_current_drug_sig)!=0){
      
      name_drug<-rownames(matrix_drugs_genes_pvalue_aov)[i]
      
      current_drug<-therapy_table[therapy_table$Drug_name%in%name_drug,c("Patient_ID","long_treatment")]
      
      input_nested3<-merge(current_drug,input_nested_rid,by.x="Patient_ID",by.y="POG_ID")
      
      input_nested4<-input_nested3[,colnames(input_nested3)%in%c("patient_id","long_treatment","hmm_states",genes_current_drug_sig)]
      
      ncol<-round(length(genes_current_drug_sig)/2)
      nrow<-round(length(genes_current_drug_sig)/2)
      
      list_plot<-vector(mode="list",length(genes_current_drug_sig))
      
      for(pg in 1:length(genes_current_drug_sig)){
        
        gene_plot<-genes_current_drug_sig[pg]
        pval_gene<-raw_genes_current_drug_sig[pg]
        
        input_rid<-input_nested4[,colnames(input_nested4)%in%c("hmm_states","long_treatment",gene_plot)]
        colnames(input_rid)[3]<-"gene"
        
        plot_list<-ggplot(input_rid, aes(x=factor(hmm_states), y=long_treatment, fill=factor(gene))) +
          geom_boxplot()+geom_point(aes(fill = factor(gene)), size = 2, shape = 21, position = position_jitterdodge())+ggtitle(paste(gene_plot,round(pval_gene,4),sep="-"))
        
        print(plot_list)
        
        list_plot[[pg]]<-plot_list
        
      }
      
      pdf(paste(drug,"aov",class1,class2,"pdf",sep="."))
      ml <- marrangeGrob(list_plot,nrow=2,ncol=2)
      print(ml)
      dev.off()
      
    }
    
  }
  
}

#
#  POG570 analysis
#
setwd(input_dir)
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
IDs <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "hgnc_symbol", values = markers_genes_read[,1],mart = mart)

markers_genes_read2<-merge(markers_genes_read,IDs,by.x="genes",by.y="hgnc_symbol")

setwd(input_dir)

tab_exp<-fread("POG570_TPM_expression.txt.gz",header=T)
tab_exp_markers<-merge(markers_genes_read2,tab_exp,by.x="ensembl_gene_id",by.y="genes")
tab_exp_markers[,-c(1:5)]<-log(tab_exp_markers[,-c(1:5)]+1,2)

# tab_exp_markers[,-c(1:5)]<-apply(tab_exp_markers[,-c(1:5)],2,FUN=function(X){(X-mean(X))/sd(X)})

all_scores_emt<-NULL

for(iemt in 6:ncol(tab_exp_markers)){
  
  mean_epi<-mean(tab_exp_markers[tab_exp_markers[,3]%in%"Epithelial_marker",iemt])
  mean_mes<-mean(tab_exp_markers[tab_exp_markers[,3]%in%"Mesenchymal_marker",iemt])
  
  score_emt<-mean_mes-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

emt_scores_df<-data.frame(Patient_ID=colnames(tab_exp_markers)[6:ncol(tab_exp_markers)],EMT_score=all_scores_emt)


dem_table <- read.xlsx(xlsxFile = "Table_S1_Demographics.xlsx", sheet = 1)
therapy_table<-read.xlsx(xlsxFile = "Table_S2_Treatment.xlsx", sheet = 1)
# convert the days of therapy in months
therapy_table$long_treatment<-(therapy_table$Therapy_end_or_biopsy_date-therapy_table$Therapy_start_date)/30
  
therapy_table2<-unique(therapy_table[,c(1,2,7,8,9,10)])

freq_therapies<-data.frame(table(therapy_table2[,2]))
drugs_to_consider<-freq_therapies[freq_therapies[,2]>=10,1]

therapy_table_to_consider<-therapy_table2[therapy_table2[,2]%in%drugs_to_consider,]

emt_with_dem<-merge(emt_scores_df,dem_table,by.x="Patient_ID",by.y="PATIENT_ID")


unique(emt_with_dem$PRIMARY_SITE)
# "Blood","Brain"

emt_with_dem2<-emt_with_dem[-which(emt_with_dem$PRIMARY_SITE%in%c("Blood","Brain")),]

input_for_chart<-merge(emt_with_dem2,therapy_table_to_consider,by.x="Patient_ID",by.y="Patient_ID")

input_for_chart$logmonth<-log(input_for_chart$long_treatment+1,2)

# input_for_chart$On_treatment_at_biopsy<-ifelse(input_for_chart$On_treatment_at_biopsy==1,"on_treatment","not_treatment")
# input_for_chart$On_treatment_at_biopsy<-as.factor(input_for_chart$On_treatment_at_biopsy)

# get the patient untreated at the enrollement
# df_pts_before_after_treatment[df_pts_before_after_treatment[,1]%in%30159,]

pts_no_treatment<-unique(input_for_chart[which(input_for_chart$On_treatment_at_biopsy%in%0),1])
df_pts_before_after_treatment<-input_for_chart[which(input_for_chart[,1]%in%pts_no_treatment),]

df_pts_before_after_treatment$On_treatment_at_biopsy<-ifelse(df_pts_before_after_treatment$On_treatment_at_biopsy==1,"after_treatment","before_treatment")
df_pts_before_after_treatment$On_treatment_at_biopsy<-factor(df_pts_before_after_treatment$On_treatment_at_biopsy,level=c("before_treatment","after_treatment"))

table(df_pts_before_after_treatment[,1],df_pts_before_after_treatment$On_treatment_at_biopsy)

setwd(output_dir)

pdf("boxplot_emt_before_after_treatment.pdf")

p2_sig<-ggplot(df_pts_before_after_treatment,aes(x=On_treatment_at_biopsy, y=EMT_score))+
geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")
p3<-p2_sig+stat_compare_means(comparisons = list(c("after_treatment","before_treatment")),method="wilcox.test")
print(p3)
  
list_drugs<-unique(df_pts_before_after_treatment$Drug_name)
  
for(lsds in list_drugs){
    
  sub_drug<-df_pts_before_after_treatment[df_pts_before_after_treatment$Drug_name%in%lsds,]
  
  p2_sig<-ggplot(sub_drug,aes(x=On_treatment_at_biopsy, y=EMT_score))+
  geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
  theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")
  p3<-p2_sig+stat_compare_means(comparisons = list(c("after_treatment","before_treatment")),method="wilcox.test")
  print(p3+ggtitle(lsds))
  
  }

dev.off()




list_drugs<-unique(input_for_chart$Drug_name)

cor_all<-NULL
pval_cor_all<-NULL

for(i in list_drugs){
  
  emt_drug<-input_for_chart[input_for_chart$Drug_name%in%i,c("EMT_score")]
  months_drug<-input_for_chart[input_for_chart$Drug_name%in%i,c("logmonth")]
  
  cor_save<-cor.test(emt_drug,months_drug,method="spearman")$estimate
  pval_save<-cor.test(emt_drug,months_drug,method="spearman")$p.value
  
  cor_all<-c(cor_all,cor_save)
  pval_cor_all<-c(pval_cor_all,pval_save)

}

dfcor<-data.frame(drugs=list_drugs,cor=cor_all,pval=pval_cor_all,pval_log=-log(pval_cor_all+1,10))
dfcor<-dfcor[order(dfcor$cor,decreasing=F),]
dfcor[,1]<-factor(dfcor[,1],levels=dfcor[,1])

dfcorsig<-dfcor[which(dfcor$pval<=0.05),]
dfcorsig<-dfcorsig[order(dfcorsig$cor,decreasing=F),]
dfcorsig[,1]<-factor(dfcorsig[,1],levels=dfcorsig[,1])

pdf("cor_POG_sig.pdf",width = 15,height=4)

# 
#  Plot correlation drugs significant
# 

cdat_sp_sig <- ggplot(dfcorsig, aes(x = 1, y = drugs,size = cor,color=pval_log)) +
geom_point() + scale_colour_viridis(option = "plasma",direction = -1)+
theme(legend.position="right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("Correlation")+ ggtitle("Correlation")

df_long_ther_cor<-input_for_chart[input_for_chart$Drug_name%in%dfcorsig[,1],]

df_long_ther_cor$Drug_name<-factor(df_long_ther_cor$Drug_name,levels=dfcorsig[,1])

p2_sig<-ggplot(df_long_ther_cor,aes(x=Drug_name, y=logmonth))+
geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
theme(text = element_text(size=8),legend.position="top")+coord_flip()+theme_bw()+ ylab("Months of Treatment (months), log(months+1,2)")+ ggtitle("Boxplot")

print(plot_grid(cdat_sp_sig,p2_sig,width=c(0.2,0.8),nrow=1))

dev.off()

# 
#  Plot correlation all drugs selected
# 
pdf("cor_POG.pdf",width = 15,height=10)

cdat_sp <- ggplot(dfcor, aes(x = 1, y = drugs,size = cor,color=pval_log)) +
geom_point() + scale_colour_viridis(option = "plasma",direction = -1)+
theme(legend.position="right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("Correlation")+ ggtitle("Correlation")

df_long_ther_cor<-input_for_chart[input_for_chart$Drug_name%in%dfcor[,1],]

df_long_ther_cor$Drug_name<-factor(df_long_ther_cor$Drug_name,levels=dfcor[,1])

p2<-ggplot(df_long_ther_cor,aes(x=Drug_name, y=logmonth))+
  geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
  theme(text = element_text(size=8),legend.position="top")+coord_flip()+theme_bw()+ ylab("Months of Treatment (months), log(months+1,2)")+ ggtitle("Boxplot")

print(plot_grid(cdat_sp,p2,width=c(0.2,0.8),nrow=1))

dev.off()

pdf("analysis_POG.pdf")

p1<-ggplot(input_for_chart, aes(x = long_treatment, y = Drug_name)) + geom_density_ridges()

print(p1)

p2<-ggplot(input_for_chart,aes(x=Drug_name, y=logmonth))+
geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
theme(text = element_text(size=8))+coord_flip()+theme_bw()+ theme(legend.position="top")+labs(title="Boxplots")

print(p2)


list_drugs<-unique(input_for_chart$Drug_name)

for(i in list_drugs){

model<-lm(EMT_score~long_treatment,input_for_chart[input_for_chart$Drug_name%in%i,])

pval<-round(coef(summary(model))[2,4],3)
rsquare<-round(summary(model)$r.squared,3)
  
psm<-ggplot(input_for_chart[input_for_chart$Drug_name%in%i,], aes(x=EMT_score, y=long_treatment))+
geom_point()+geom_smooth(method=lm, se = TRUE)+theme_bw()+ylab("months (log2)")+xlab(paste("EMT score"))+ggtitle(paste(i,",","R2:",rsquare,"p-value:",pval))

print(psm)

}
dev.off()


# 
#  Study the relations between drugs and HMM states of EMT
# 

load("Pseudotime_POG_on_MCF.RData")

list_drugs<-unique(input_for_chart$Drug_name)

pval_aov_all<-NULL
out_drugs<-NULL

for(i in list_drugs){
  
  df_id_with_timetreat<-input_for_chart[input_for_chart$Drug_name%in%i,c("Patient_ID","long_treatment")]
  
  emt_states<-POG_pseudotime_EMT_metpot_HMM[POG_pseudotime_EMT_metpot_HMM[,1]%in%df_id_with_timetreat[,1],c(1,4)]  
  
  input_aov<-merge(emt_states,df_id_with_timetreat,by.x="POG_ID",by.y="Patient_ID")
  
  if(length(unique(input_aov[,2]))!=1){
    
  res.aov <- aov(long_treatment ~ hmm_states, data = input_aov)
  pval<-summary(res.aov)[[1]][[1,"Pr(>F)"]]
    
  pval_aov_all<-c(pval_aov_all,pval)
  out_drugs<-c(out_drugs,i)
    
  }
  
}

dfaov<-data.frame(drugs=out_drugs,pval=pval_aov_all,qvalue=p.adjust(pval_aov_all))

significant_drugs<-dfaov[dfaov[,2]<=0.05,1]

df_selected_drug_treatment<-input_for_chart[input_for_chart$Drug_name%in%significant_drugs,c("Patient_ID","Drug_name","long_treatment","logmonth")]  
input_boxplot<-merge(POG_pseudotime_EMT_metpot_HMM[,c(1,4)],df_selected_drug_treatment,by.x="POG_ID",by.y="Patient_ID")

input_boxplot$hmm_states<-gsub(input_boxplot$hmm_states,pattern = "1",replacement="epi")
input_boxplot$hmm_states<-gsub(input_boxplot$hmm_states,pattern = "2",replacement="mes")
input_boxplot$hmm_states<-gsub(input_boxplot$hmm_states,pattern = "3",replacement="hEMT")

pdf("POG_HMM_significant_drugs.pdf")

for(si in significant_drugs){

my_comparisons<-list(c("mes","epi"),c("epi","hEMT"),c("mes","hEMT"))
  
p1<-ggviolin(input_boxplot[input_boxplot$Drug_name%in%si,], x = "hmm_states", y = "long_treatment", fill = "hmm_states",
  palette = c("#56b4e9", "#E69f00", "#EE0C0C"),
  add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test" ,label = "p.signif")+ggtitle(si)

print(p1)

}
dev.off()


# 
# Significant groups
# 

setwd(input_dir)

load("Pseudotime_POG_on_MCF.RData")

POG_pseudotime_EMT_metpot_HMM$hmm_states<-gsub(POG_pseudotime_EMT_metpot_HMM$hmm_states,pattern = "1",replacement="epi")
POG_pseudotime_EMT_metpot_HMM$hmm_states<-gsub(POG_pseudotime_EMT_metpot_HMM$hmm_states,pattern = "2",replacement="mes")
POG_pseudotime_EMT_metpot_HMM$hmm_states<-gsub(POG_pseudotime_EMT_metpot_HMM$hmm_states,pattern = "3",replacement="hEMT")

my_comparisons<-list(c("mes","epi"),c("epi","hEMT"),c("mes","hEMT"))

print(length(which(POG_pseudotime_EMT_metpot_HMM$PRIMARY_SITE%in%c("Blood","Brain"))))

setwd(output_dir)

for(myc in 1:length(my_comparisons)){

current_comparison<-my_comparisons[[myc]]

pdf(paste("POG_HMM_significant_drugs",current_comparison[1],current_comparison[2],"pdf",sep="."))

for(si in unique(input_for_chart$Drug_name)){

  df_id_with_timetreat<-input_for_chart[input_for_chart$Drug_name%in%si,c("Patient_ID","long_treatment")]
  
  POG_pseudotime_EMT_metpot_HMM2<-POG_pseudotime_EMT_metpot_HMM[which(POG_pseudotime_EMT_metpot_HMM$hmm_states %in%current_comparison),]
  
  emt_states<-POG_pseudotime_EMT_metpot_HMM2[POG_pseudotime_EMT_metpot_HMM2[,1]%in%df_id_with_timetreat[,1],c(1,4)]  
  
  input_comp<-merge(emt_states,df_id_with_timetreat,by.x="POG_ID",by.y="Patient_ID")
  
  print(current_comparison)
  
  groups_for_stat<-current_comparison
  
  input_comp$hmm_states<-as.factor(input_comp$hmm_states)
    
  p1<-ggplot(input_comp, aes(x=hmm_states, y=long_treatment, fill=hmm_states)) +
  geom_boxplot()+geom_point(aes(fill = factor(hmm_states)), size = 2, shape = 21, position = position_jitterdodge())+theme_bw()+ggtitle(si)
  
  print(p1+stat_compare_means(comparisons = list(groups_for_stat),method="wilcox.test"))
}
dev.off()

}

# 
#  Do the overall survival analysis only for the patients in the three groups
# 

setwd(input_dir)

surv_table<- read.xlsx(xlsxFile = "Table_S2_Treatment.xlsx", sheet = 3)
surv_table2<-merge(surv_table,POG_pseudotime_EMT_metpot_HMM,by.x="Patient_ID",by.y="POG_ID")

surv_table2$hmm_states<-gsub(surv_table2$hmm_states,pattern = "1",replacement="epi")
surv_table2$hmm_states<-gsub(surv_table2$hmm_states,pattern = "2",replacement="mes")
surv_table2$hmm_states<-gsub(surv_table2$hmm_states,pattern = "3",replacement="hEMT")

sfit <- survfit(Surv(as.numeric(surv_table2$Overall_survival_days)/365, event = surv_table2$Alive_0_Death_1) ~ surv_table2$hmm_states,data=surv_table2)

setwd(output_dir)

pdf("HMM_POG570.overallsurvival.pdf")
p4<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               legend.labs=c("epi", "mes","hEMT"), legend.title="biological_states",
               palette=c("dodgerblue2", "red2","orange2"),
               title="Kaplan-Meier Curve for epi, mes, hEMT",
               risk.table.height=0.3)
print(p4)
dev.off()

# 
#  Mutational data: look at biomarkers discovered
# 

setwd(input_dir)
cosmic_table<-read.csv(file="Census_allMon Aug 17 13_43_26 2020.csv",stringsAsFactors=F)
genes_to_consider<-cosmic_table[,1]

# Import the mutations in PO570
therapy_table<-read.xlsx(xlsxFile = "Table_S2_Treatment.xlsx", sheet = 1)
therapy_table$long_treatment<-(therapy_table$Therapy_end_or_biopsy_date-therapy_table$Therapy_start_date)/30

tab_mut<-fread("POG570_small_mutations.txt.gz",data.table=F)
tab_mut_cosmic<-tab_mut[which(tab_mut$gene_id%in%genes_to_consider),]
#remove introns variants
tab_mut_cosmic2<-tab_mut_cosmic[-which(tab_mut_cosmic$effect%in%"intron_variant"),]

# mut_matrix_drugs<-merge(tab_mut_cosmic2,therapy_table,by.x="patient_id",by.y="Patient_ID")

freq_therapies<-data.frame(table(therapy_table2[,2]))
drugs_to_consider<-freq_therapies[freq_therapies[,2]>=10,1]

# mut_matrix_drugs2<-mut_matrix_drugs[mut_matrix_drugs$Drug_name%in%drugs_to_consider,]

POG_pseudotime_EMT_metpot_HMM$hmm_states<-gsub(POG_pseudotime_EMT_metpot_HMM$hmm_states,pattern = "1",replacement="epi")
POG_pseudotime_EMT_metpot_HMM$hmm_states<-gsub(POG_pseudotime_EMT_metpot_HMM$hmm_states,pattern = "2",replacement="mes")
POG_pseudotime_EMT_metpot_HMM$hmm_states<-gsub(POG_pseudotime_EMT_metpot_HMM$hmm_states,pattern = "3",replacement="hEMT")

# mut_matrix_drugs2$patient_id<-as.character(mut_matrix_drugs2$patient_id)

# input_nested<-merge(mut_matrix_drugs2,POG_pseudotime_EMT_metpot_HMM,by.x="patient_id",by.y="POG_ID")

df_patients_with_genes<-unique(tab_mut_cosmic2[,c(11,12)])
matrix_patients_for_genes<-as.data.frame.matrix(table(df_patients_with_genes))
matrix_patients_for_genes2<-data.frame(patientID=rownames(matrix_patients_for_genes),matrix_patients_for_genes)
matrix_patients_for_genes2_with_hmm<-merge(POG_pseudotime_EMT_metpot_HMM[,c(1,4)],matrix_patients_for_genes2,by.x="POG_ID",by.y="patientID")

list_comparisons<-list(c("mes","epi"),c("hEMT","epi"))

setwd(output_dir)

for(i in 1:length(list_comparisons)){
  
  groups_to_compare<-list_comparisons[[i]]
  
  input_nested_rid<-matrix_patients_for_genes2_with_hmm[matrix_patients_for_genes2_with_hmm$hmm_states%in%groups_to_compare,]
  
  list_drugs_results<-vector(mode="list",length(drugs_to_consider))
  
  for(d in 1:length(drugs_to_consider)){
  
  d2<-drugs_to_consider[d]
    
  current_drug<-therapy_table[therapy_table$Drug_name%in%d2,c("Patient_ID","long_treatment")]
  
  input_nested3<-merge(current_drug,input_nested_rid,by.x="Patient_ID",by.y="POG_ID")
  
  genes_to_consider<-apply(input_nested3[,-c(1:3)],2,max)
  genes_to_consider<-names(genes_to_consider[which(genes_to_consider!=0)])
  
  input_nested4<-input_nested3[,colnames(input_nested3)%in%c("patient_id","long_treatment","hmm_states",genes_to_consider)]
  
    pvalue_current_drug<-NULL
    
    for(g in colnames(input_nested4)[-c(1:3)]){
      
    gene_var<-as.character(g)
    
    print(gene_var)
    
    if(length(unique(input_nested4$hmm_states))>1 & length(unique(input_nested3[,colnames(input_nested3)%in%c(gene_var)])) == 2){
      
    fq<-sum( table(input_nested3[,colnames(input_nested3)%in%c("hmm_states",gene_var)]))
    
    check_median<-input_nested3[,colnames(input_nested3)%in%c("hmm_states","long_treatment",as.character(gene_var))]
    check_median$new_id<-paste(check_median$hmm_states,check_median[,colnames(check_median)%in%gene_var],sep="_")
      
    custom_fun<-function(X){quantile(X)[4]}
    
    median_longtherapy<-aggregate(check_median$long_treatment, by=list(Category=check_median$new_id), FUN=custom_fun)
    
    idx_max<-which.max(median_longtherapy[,2])
    
    print(median_longtherapy)
    
    group_max<-median_longtherapy[idx_max,1]
    
    group_reference<-paste(groups_to_compare[1],"_1",sep="")
    
        if(group_max==group_reference & fq > 50){
        
            formula_to_use<-as.formula(paste("long_treatment ~ hmm_states /", gene_var))
              
            stat_aov<-aov(formula_to_use,input_nested3)
          
            pvalue<-as.numeric(summary(stat_aov)[[1]][2,]["Pr(>F)"])
            
            pvalue_current_drug<-c(pvalue_current_drug,pvalue)
        
        }else{
          
            pvalue_current_drug<-c(pvalue_current_drug,1)
        
        }
        
    
    }else{
      
      pvalue_current_drug<-c(pvalue_current_drug,1)
      
    }
    
    
    }
    
    names(pvalue_current_drug)<- colnames(input_nested4)[-c(1:3)]
    
    list_drugs_results[[d]]<-as.data.frame(t(data.frame(pvalue_current_drug)))
    
  
  }
  
   matrix_drugs_genes_pvalue_aov<-do.call(rbind.fill,list_drugs_results)
   rownames(matrix_drugs_genes_pvalue_aov)<- drugs_to_consider
  
   plotSignificantDrugsGenes(matrix_drugs_genes_pvalue_aov,input_nested_rid,groups_to_compare[1],groups_to_compare[2],therapy_table)
   
}


