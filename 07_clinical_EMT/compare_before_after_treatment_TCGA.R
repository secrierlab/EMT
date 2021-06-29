library(sva)
library(biomaRt)
library(ggpubr)
library(openxlsx)
library(data.table)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)

prepareCombat<-function(dat,batch){        
  
  dat <- as.matrix(dat)
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(dat[, batch == batch_level], 1,
                         function(x) {
                           var(x) == 0
                         })))
    }else {
      
      return(which(rep(1, 3) == 2))
      
    }       
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n",
                length(zero.rows)))
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
}
#
# Load canonical marker genes for EMT
#

setwd("/home/guidantoniomt/pseudospace/HMM")
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)

# 
# Compare treatment before and after with the TCGA data
#


knn_df_tcga<-read.delim(file="/home/guidantoniomt/pseudospace/HMM/HMM_results_nstates_tumors_for_states3.withEMT.txt",stringsAsFactors=F)
knn_df_tcga$samples2<-knn_df_tcga$samples
knn_df_tcga$samples<- unlist(lapply(strsplit(as.character(knn_df_tcga$samples),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")}))


# 
# Load the TCGA dataset
# 

setwd("/home/guidantoniomt/pseudospace/input_pseudospace")
load("TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")

tcga<-TCGA_GEXP_ALL
tcga[,-1]<-log(tcga[,-1]+1,2)
tcga<-setDT(tcga)

tcga2 <- data.frame(tcga[, lapply(.SD, mean), by = gene_symbol])
rownames(tcga2)<-tcga2[,1]

columns_to_use<-c("gene_symbol",knn_df_tcga$samples2)
tcga3<-tcga2[,which(colnames(tcga2)%in%columns_to_use)]

therapy_table_tcga<-fread("/home/guidantoniomt/pseudospace/survival_analysis/media-2.tsv",data.table=F)
therapy_table_tcga$patient_barcode<-toupper(therapy_table_tcga$patient_barcode)

therapy_table_tcga_uni<-unique(therapy_table_tcga[,c("patient_barcode","drug_name")])

# 
# Load the POG dataset
# 
dem_table <- read.xlsx(xlsxFile = "/home/guidantoniomt/pseudospace/POG570/Table_S1_Demographics.xlsx", sheet = 1)
samples_to_remove<-dem_table[which(dem_table$PRIMARY_SITE%in%c("Blood","Brain")),1]

therapy_table<-read.xlsx(xlsxFile = "/home/guidantoniomt/pseudospace/POG570/Table_S2_Treatment.xlsx", sheet = 1)

setwd("/home/guidantoniomt/pseudospace/POG570")
tab_exp<-fread("POG570_TPM_expression.txt.gz",header=T)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
IDs <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "ensembl_gene_id", values = tab_exp[,1],mart = mart)

matrix_pog<-merge(IDs,tab_exp,by.x="ensembl_gene_id",by.y="genes")
matrix_pog[,-c(1:2)]<-log(matrix_pog[,-c(1:2)]+1,2)

# 
# Combat: Merge TCGA with POG 
# 

matrix_tcga_pog<-setDT(merge(tcga3,matrix_pog[,-1],by.x="gene_symbol",by.y="hgnc_symbol"))

matrix_tcga_pog2 <- as.data.frame(matrix_tcga_pog[, lapply(.SD, mean), by = gene_symbol])

rownames(matrix_tcga_pog2)<-matrix_tcga_pog2[,1]

matrix_tcga_pog3<-as.matrix(matrix_tcga_pog2[,-1])

classes_for_combat<-c(rep("tcga",ncol(tcga3[,-1])),rep("POG",ncol(matrix_pog[,-c(1:2)])))

mat_input_combat<-prepareCombat(matrix_tcga_pog3, batch=classes_for_combat)

matcombat <- ComBat(mat_input_combat, batch=classes_for_combat)

zscore<-function(x){(x-mean(x))/sd(x)}

rownames_save<-rownames(matcombat)
colnames_save<-colnames(matcombat)

matcombat<-t(apply(matcombat,1,zscore))

rownames(matcombat)<-rownames_save

matcombat2<-data.frame(genes=rownames(matcombat),matcombat)

matcombat_emt<-merge(markers_genes_read,matcombat2,by.x="genes",by.y="genes")


#
# Now I can re-compute the EMT scores
#

all_scores_emt<-NULL

for(iemt in 5:ncol(matcombat_emt)){
  
  mean_epi<-mean(matcombat_emt[matcombat_emt[,2]%in%"Epithelial_marker",iemt])
  mean_mes<-mean(matcombat_emt[matcombat_emt[,2]%in%"Mesenchymal_marker",iemt])
  
  score_emt<-mean_mes-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

dfEMTscores<-data.frame(ID=colnames(matcombat_emt)[-c(1:4)],EMT_score=all_scores_emt)
dfEMTscores$status<-rep("before_treatment",nrow(dfEMTscores))
dfEMTscores[grep(dfEMTscores[,1],pattern="TCGA",invert=T),3]<-"after_treatment"

#merge with the HMM states in TCGA
dfEMT_tcga_hmm<-merge(dfEMTscores[grep(dfEMTscores[,1],pattern="TCGA"),],knn_df_tcga[,c(5,4)],by.x="ID",by.y="samples2")
colnames(dfEMT_tcga_hmm)[4]<-"hmm_states"
  
dfEMT_tcga_hmm$hmm_states<-as.character(dfEMT_tcga_hmm$hmm_states)

dfEMT_tcga_hmm$hmm_states<-gsub(dfEMT_tcga_hmm$hmm_states,pattern="1",replacement="hEMT")
dfEMT_tcga_hmm$hmm_states<-gsub(dfEMT_tcga_hmm$hmm_states,pattern="2",replacement="epi")
dfEMT_tcga_hmm$hmm_states<-gsub(dfEMT_tcga_hmm$hmm_states,pattern="3",replacement="mes")

dfEMT_tcga_hmm$ID<- unlist(lapply(strsplit(as.character(dfEMT_tcga_hmm$ID),split="\\."),FUN=function(X){paste(X[2:4],collapse="-")}))

#merge with the HMM states in POG
load("/home/guidantoniomt/pseudospace/POG570/Pseudotime_POG_on_MCF.RData")
dfEMTscores[,1]<-gsub(dfEMTscores[,1],pattern="^X",replacement="")
             
dfEMT_pog_hmm<-merge(dfEMTscores[grep(dfEMTscores[,1],pattern="TCGA",invert=T),],POG_pseudotime_EMT_metpot_HMM[,c(1,4)],by.x="ID",by.y="POG_ID")

dfEMT_pog_hmm$hmm_states<-gsub(dfEMT_pog_hmm$hmm_states,pattern = "1",replacement="epi")
dfEMT_pog_hmm$hmm_states<-gsub(dfEMT_pog_hmm$hmm_states,pattern = "2",replacement="mes")
dfEMT_pog_hmm$hmm_states<-gsub(dfEMT_pog_hmm$hmm_states,pattern = "3",replacement="hEMT")
# dfEMT_pog_hmm<-dfEMT_pog_hmm[-which(dfEMT_pog_hmm[,1]%in%samples_to_remove),]

dfEMT_global<-rbind(dfEMT_tcga_hmm,dfEMT_pog_hmm)

dfEMT_global$status<-factor(dfEMT_global$status,levels=c("before_treatment","after_treatment"))

setwd("/home/guidantoniomt/pseudospace/POG570/")

pdf("boxplot_EMT_score_TCGA_POG_before_after_treatment.pdf")
p2_sig<-ggplot(dfEMT_global,aes(x=status, y=EMT_score))+
  geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
  theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")
p3<-p2_sig+stat_compare_means(comparisons = list(c("before_treatment","after_treatment")),method="wilcox.test")
print(p3)
dev.off()


# check which are the common treatments between POG and TCGA


therapies_tcga<-unique(therapy_table_tcga[,c("patient_barcode","drug_name")])
colnames(therapies_tcga)<-c("ID","drugs")

therapies_pog<-unique(therapy_table[,c(1,2)])
colnames(therapies_pog)<-c("ID","drugs")

combined_datasets_drugs<-rbind(therapies_tcga,therapies_pog)
common_drugs<-intersect(therapies_tcga[,2],therapies_pog[,2])

combined_datasets_drugs2<-combined_datasets_drugs[which(combined_datasets_drugs[,2]%in%common_drugs),]

dfEMT_global_drugs<-merge(dfEMT_global,combined_datasets_drugs2,by.x="ID",by.y="ID")

dfEMT_global_drugs$status<-factor(dfEMT_global_drugs$status,levels=c("before_treatment","after_treatment"))

pdf("boxplot_EMT_score_TCGA_POG_before_after_treatment_v2.pdf")

stat.test <- dfEMT_global_drugs %>%
  group_by(hmm_states) %>%
  wilcox_test(EMT_score ~ status)

stat.test <- stat.test %>% 
  add_xy_position(x = "hmm_states", dodge = 0.8)

p2_sig<-ggplot(dfEMT_global,aes(x=hmm_states, y=EMT_score,color=status))+
  geom_boxplot()+ stat_pvalue_manual(stat.test,  label = "p", tip.length = 0)+
  theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")

stat.test <- dfEMT_global_drugs %>%
  anova_test(EMT_score ~ hmm_states)

p3<-p2_sig+ggtitle(paste("p-value:",paste(stat.test[,5],stat.test[,6])))

print(p3)

dev.off()

#
# Select the drugs for downstream analysis
#

#table of drugs for status (before and after)
freq_drugs_between_treatment<-table(dfEMT_global_drugs$drugs,dfEMT_global_drugs$status)
#take the drugs in more than 10 patients between the two groups
drugs_to_consider<-names(which(rowSums(freq_drugs_between_treatment)>10))
#each group at minimum must have 10 patients
idx1<-names(which(freq_drugs_between_treatment[drugs_to_consider,1]>=10))
idx2<-names(which(freq_drugs_between_treatment[drugs_to_consider,2]>=10))
final_drug<-intersect(idx1,idx2)

pdf("boxplot_EMT_score_TCGA_POG_before_after_treatment_drug_by_drug.pdf",width=10,height=5)

dfEMT_global_drugs<-dfEMT_global_drugs[which(dfEMT_global_drugs$drugs%in%final_drug),]

stat.test <- dfEMT_global_drugs %>%
  group_by(drugs) %>%
  wilcox_test(EMT_score ~ status)


stat.test<-data.frame(stat.test)
list_drugs<-stat.test[stat.test$p<=0.05,1]

for(pqd in list_drugs){
  
  p2_sig<-ggplot(dfEMT_global_drugs[dfEMT_global_drugs$drugs%in%pqd,],aes(x=status, y=EMT_score))+
    geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
    theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")
  
  p3<-p2_sig+stat_compare_means(comparisons = list(c("before_treatment","after_treatment")),method="wilcox.test")
  p3<-p3+ geom_hline(yintercept=0, linetype="dashed", color = "red")+ggtitle(pqd)

  # i want to see how much samples there are in a specific hmm state and status (Before, after)
  tab_plot<-tableGrob(as.matrix(table(dfEMT_global_drugs[dfEMT_global_drugs$drugs%in%pqd,c("status","hmm_states")])))
  
  print(grid.arrange(p3,tab_plot,ncol=2))
  
}

dev.off()

#
# With only MES and hEMT states
#

pdf("boxplot_EMT_score_TCGA_POG_before_after_treatment_drug_by_drug.ONLY_MES_hEMT.pdf",width=10,height=5)

dfEMT_global_drugs<-dfEMT_global_drugs[which(dfEMT_global_drugs$drugs%in%final_drug),]
dfEMT_global_drugs_transformed<-dfEMT_global_drugs[dfEMT_global_drugs$hmm_states%in%c("hEMT","mes"),]
  
stat.test <- dfEMT_global_drugs_transformed %>%
  group_by(drugs) %>%
  wilcox_test(EMT_score ~ status)

stat.test<-data.frame(stat.test)
list_drugs<-stat.test[stat.test$p<=0.05,1]

for(pqd in list_drugs){
  
  p2_sig<-ggplot(dfEMT_global_drugs_transformed[dfEMT_global_drugs_transformed$drugs%in%pqd,],aes(x=status, y=EMT_score))+
    geom_boxplot()+geom_jitter(width=0.1, aes(color = EMT_score))+scale_colour_gradient2(low = "orange2", mid = "grey", high = "darkmagenta", midpoint = 0)+
    theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")
  
  p3<-p2_sig+stat_compare_means(comparisons = list(c("before_treatment","after_treatment")),method="wilcox.test")
  p3<-p3+ geom_hline(yintercept=0, linetype="dashed", color = "red")+ggtitle(pqd)
  
  # i want to see how much samples there are in a specific hmm state and status (Before, after)
  tab_plot<-tableGrob(as.matrix(table(dfEMT_global_drugs_transformed[dfEMT_global_drugs_transformed$drugs%in%pqd,c("status","hmm_states")])))
  
  print(grid.arrange(p3,tab_plot,ncol=2))
  
}

dev.off()

#
# Anova nested
#
library(rstatix)

list_drugs<-unique(dfEMT_global_drugs$drugs)

pval_aov_all<-NULL
pval_aov_all2<-NULL
out_drugs<-NULL

for(ipa in list_drugs){

input_aov<-dfEMT_global_drugs[which(dfEMT_global_drugs$drugs%in%ipa),]

if(length(unique(input_aov$hmm_states))>=2 & length(unique(input_aov$status))>=2 & nrow(input_aov)>10){
  
  print(ipa)
  
  input_aov$hmm_states<-as.factor(input_aov$hmm_states)
  input_aov$status<-as.factor(input_aov$status)
  
  formula_to_use<-as.formula(paste("EMT_score ~ hmm_states /","factor(status)"))
  
  res.aov <- aov(input_aov$EMT_score ~ input_aov$hmm_states / factor(input_aov$status), data = input_aov)
  
  pval<-summary(res.aov)[[1]][[1,"Pr(>F)"]]
  pval2<-summary(res.aov)[[1]][[2,"Pr(>F)"]]
  
  pval_aov_all<-c(pval_aov_all,pval)
  pval_aov_all2<-c(pval_aov_all2,pval2)
  
  out_drugs<-c(out_drugs,unique(input_aov$drugs))

}

}

df_drugs<-data.frame(drugs=out_drugs,pvalue=pval_aov_all,pvalue2=pval_aov_all2)
df_drugs2<-df_drugs[which(df_drugs$pvalue<=0.05 & df_drugs$pvalue2<=0.05),]
# df_drugs$BH_pval<-p.adjust(df_drugs[,2],"BH")
# df_drugs$BH_pval2<-p.adjust(df_drugs[,3],"BH")

list_drugs2<-unique(df_drugs2[,1])

# SELECT MANUALLY THE SIGNIFICANT GENES
# list_drugs2<-c("FLUOROURACIL","BEVACIZUMAB","EPIRUBICIN","DOCETAXEL","TAMOXIFEN","CYCLOPHOSPHAMIDE","CAPECITABINE","CISPLATIN","PACLITAXEL","LETROZOLE","DOXORUBICIN","GEMCITABINE","ETOPOSIDE","CARBOPLATIN","VINCRISTINE")

pdf("boxplot_EMT_score_TCGA_POG_before_after_treatment_drug_by_drug.ANOVA.pdf")

for(pqd2 in list_drugs2){
  
  pval_for_main<-df_drugs2[df_drugs2[,1]%in%pqd2,2]
  
  stat.test <- dfEMT_global_drugs %>%
    group_by(hmm_states) %>%
    t_test(EMT_score ~ status)

  stat.test <- stat.test %>% 
    add_xy_position(x = "hmm_states", dodge = 0.8)
  
  p2_sig_drug<-ggplot(dfEMT_global_drugs[dfEMT_global_drugs$drugs%in%pqd2,],aes(x=hmm_states, y=EMT_score,color=status))+
    geom_boxplot()+ stat_pvalue_manual(stat.test,  label = "p", tip.length = 0)+
    theme(text = element_text(size=8),legend.position="top")+theme_bw()+ ylab("EMT score")
  
  p3_drug<-p2_sig_drug+stat_compare_means(comparisons = list(c("before_treatment","after_treatment")),method="wilcox.test")
  p4_drug<-p3_drug+ geom_hline(yintercept=0, linetype="dashed", color = "red")+ggtitle(paste(pqd2,":",pval_for_main))
  
  tab_plot<-tableGrob(table(dfEMT_global_drugs[dfEMT_global_drugs$drugs%in%pqd2,c("status","hmm_states")]))
  
  print(grid.arrange(p4_drug,tab_plot,nrow=2))
  
}

dev.off()


