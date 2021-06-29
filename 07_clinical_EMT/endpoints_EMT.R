library(TCGAbiolinks)
library("survival")
library("survminer")
source("ggforest2.R")
library("readxl")
library(plyr)

setwd("/home/guidantoniomt/pseudospace/HMM")

tab_hmm<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)

list_tissue<-paste("TCGA-",unique(sapply(strsplit(tab_hmm[,1],split="\\."),"[[",1)),sep="")
tab_hmm[,1]<-unlist(lapply(strsplit(tab_hmm[,1],split="\\."),FUN=function(x){paste(x[2:4],collapse="-")}))

# 
# Donwload  general clinical information
# 

TCGA_global_tcga<-vector(mode="list",length(list_tissue))

for(tissue in 1:length(list_tissue)){
  
  print(tissue)
  
  clin.tissue<-GDCquery_clinic(list_tissue[tissue], "clinical")
  
  col_to_use<-c("submitter_id","vital_status",grep(colnames(clin.tissue),pattern="days",value=T))
  
  # clin.tissue2<-clin.tissue[,which(colnames(clin.tissue)%in%col_to_use)]
  
  clin.tissue2<-clin.tissue[,which(colnames(clin.tissue)%in%c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up","tumor_stage","synchronous_malignancy","age_at_index","gender","cigarettes_per_day","alcohol_history","bmi"))]
  
  # days_to_follow_up: Number of days between the date used for index and the date of the patient's last follow-up appointment or contact.
  # 
  
  TCGA_global_tcga[[tissue]]<-clin.tissue2
  
}


TCGA_general_clinical<-rbind.fill(TCGA_global_tcga)


# 
# Donwload  the information about the drug treatment
# 

TCGA_global_tcga<-vector(mode="list",length(list_tissue))

for(tissue in 1:length(list_tissue)){
  
  print(tissue)
  
  clin.tissue<-GDCquery(list_tissue[tissue], 
                        "Clinical",
                        data.type = "Clinical Supplement",
                        data.format = "BCR Biotab")
  GDCdownload(clin.tissue)
  
  clinical.BCRtab<- GDCprepare(clin.tissue)
  
  columns_to_use<-paste("clinical_drug",tolower(gsub(list_tissue[tissue],pattern="TCGA-",replacement="")),sep="_")
  
  clinical.drug <- data.frame(clinical.BCRtab[[columns_to_use]])
  colnames(clinical.drug)<-clinical.drug[1,]
  clinical.drug<-clinical.drug[-c(1:2),]

  TCGA_global_tcga[[tissue]]<-clinical.drug
  
}

TCGA_global_drugs<-rbind.fill(TCGA_global_tcga)

setwd("/home/guidantoniomt/pseudospace/survival_analysis")

library(readxl)
library(survival)
library(survminer)

my_data <- read_excel("1-s2.0-S0092867418302290-mmc1.xlsx",sheet="TCGA-CDR")
my_data2 <- my_data[,-1]
TCGA_surv<-merge(my_data2,tab_hmm[,c(1,5)],by.x="bcr_patient_barcode",by.y="samples")
TCGA_surv$biological_states<-as.factor(TCGA_surv$biological_states)
TCGA_surv$biological_states2<-as.character(TCGA_surv$biological_states)
TCGA_surv[which(TCGA_surv$biological_states2%in%"mes"),"biological_states2"]<-"all_mes"
TCGA_surv[which(TCGA_surv$biological_states2%in%"mix"),"biological_states2"]<-"all_mes"
TCGA_surv_save<-TCGA_surv

#TCGA_surv<-merge(TCGA_surv,TCGA_global_drugs,by.x="bcr_patient_barcode",by.y="bcr_patient_barcode",all.y=F)

#
#TCGA_surv contains already the column gender and vital status
#
TCGA_surv<-merge(TCGA_surv,TCGA_general_clinical[,-which(colnames(TCGA_general_clinical)%in%c("gender","vital_status"))],by.x="bcr_patient_barcode",by.y="submitter_id",all.y=F)
dim(TCGA_surv)

# Parse the data for coxph

TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="i/ii nos",value=F),"tumor_stage"]<-"StageI"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="not reported",value=F),"tumor_stage"]<-"not_reported"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="stage 0",value=F),"tumor_stage"]<-"Stage0"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage i$",value=F),"tumor_stage"]<-"StageI"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage ia$",value=F),"tumor_stage"]<-"StageI"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage ib$",value=F),"tumor_stage"]<-"StageI"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage ii$",value=F),"tumor_stage"]<-"StageII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iia$",value=F),"tumor_stage"]<-"StageII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iib$",value=F),"tumor_stage"]<-"StageII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iic$",value=F),"tumor_stage"]<-"StageII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iii$",value=F),"tumor_stage"]<-"StageIII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iiia$",value=F),"tumor_stage"]<-"StageIII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iiib$",value=F),"tumor_stage"]<-"StageIII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iiic$",value=F),"tumor_stage"]<-"StageIII"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iv$",value=F),"tumor_stage"]<-"StageIV"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage iva$",value=F),"tumor_stage"]<-"StageIV"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage ivb$",value=F),"tumor_stage"]<-"StageIV"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage ivc$",value=F),"tumor_stage"]<-"StageIV"
TCGA_surv[grep(TCGA_surv$tumor_stage,pattern="^stage x$",value=F),"tumor_stage"]<-"not_reported" # of these samples we don't know nothing they should be undetermined, i have called them not_reported

# TCGA_surv<-TCGA_surv[-which(TCGA_surv$tumor_stage%in%"not_reported"),]
TCGA_surv$biological_states<-TCGA_surv$biological_states
TCGA_surv$biological_states2<-TCGA_surv$biological_states2

# TCGA_surv$tumor_stage<-ifelse(TCGA_surv$tumor_stage=="High",1,0)
# synchronous_malignancy: No (7689), Not Reported (1043), Yes (42)
# tumor_stage: High 2441, Low 4151, not_reported 2170, undetermined (12)
# days_to_last_follow_up
# cigarettes_per_day: numeric (1,2,3,4,5,6)
# alcohol_history: yes (547), no (269), not reported (7958)
# bmi: numeric 
# gender: female, male
# vital_status: female (2), male (1)
# age_at_index: numeric
# days_to_death: female, male
# biological_states: epi, mes...

# 
# MES, EPI, MIX
# 
library(finalfit)
library(knitr)
library(kableExtra)

sfit <- survfit(Surv(as.numeric(TCGA_surv$OS.time)/365, event = TCGA_surv$OS) ~ TCGA_surv$biological_states,data=TCGA_surv)

setwd("/home/guidantoniomt/pseudospace/survival_analysis")

# 
# Overall survival: OS
# 

pdf("OS_CELL_Liu_hmm_states.May.pdf")
p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               legend.title="biological_states",
               palette=c("dodgerblue2", "red2","orange2"),
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p3)

p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p2)
dev.off()    

TCGA_surv_coxph<-TCGA_surv[-which(TCGA_surv$tumor_stage%in%"not_reported"),]
  
res.cox <- coxph(Surv(OS.time, OS) ~ strata(biological_states)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv_coxph)
waldp<-summary(res.cox)$waldtest[3]

explanatory = c("biological_states","age_at_index", "gender", "bmi", "tumor_stage")
dependent = "Surv(OS.time, OS)"
TCGA_surv_coxph %>%
  finalfit.coxph(dependent, explanatory) -> t6
knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph OS EPI, pEMT, MES",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"OS_CELL_Liu_hmm_states.coxph.May.pdf")

# 
# Disease free interval: DFI
# 

sfit <- survfit(Surv(as.numeric(TCGA_surv$DFI.time)/365, event = TCGA_surv$DFI) ~ TCGA_surv$biological_states,data=TCGA_surv)

pdf("DFI_CELL_Liu_hmm_states.May.pdf")
p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
               palette=c("dodgerblue2", "red2","orange2"),
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p3)
p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p2)
dev.off()  

TCGA_surv_coxph<-TCGA_surv[-which(TCGA_surv$tumor_stage%in%"not_reported"),]

res.cox <- coxph(Surv(DFI.time, DFI) ~ strata(biological_states)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv_coxph)
waldp<-summary(res.cox)$waldtest[3]

explanatory = c("biological_states","age_at_index", "gender", "bmi", "tumor_stage")
dependent = "Surv(DFI.time, DFI)"
TCGA_surv_coxph %>%
  finalfit.coxph(dependent, explanatory) -> t6
knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph DFI EPI, pEMT, MES",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"DFI_CELL_Liu_hmm_states.coxph.May.pdf")

sfit <- survfit(Surv(as.numeric(TCGA_surv$DSS.time)/365, event = TCGA_surv$DSS) ~ TCGA_surv$biological_states,data=TCGA_surv)

# 
# DSS
# 

pdf("DSS_CELL_Liu_hmm_states.May.pdf")
p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
               palette=c("dodgerblue2", "red2","orange2"),
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p3)
p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p2)
dev.off()    

TCGA_surv_coxph<-TCGA_surv[-which(TCGA_surv$tumor_stage%in%"not_reported"),]

res.cox <- coxph(Surv(DSS.time, DSS) ~ strata(biological_states)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv_coxph)
waldp<-summary(res.cox)$waldtest[3]

explanatory = c("biological_states","age_at_index", "gender", "bmi", "tumor_stage")
dependent = "Surv(DSS.time, DSS)"
TCGA_surv_coxph %>%
  finalfit.coxph(dependent, explanatory) -> t6
knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph DSS EPI, pEMT, MES",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"DSS_CELL_Liu_hmm_states.coxph.May.pdf")

sfit <- survfit(Surv(as.numeric(TCGA_surv$PFI.time)/365, event = TCGA_surv$PFI) ~ TCGA_surv$biological_states,data=TCGA_surv)

# 
# Progression free interval
# 

pdf("PFI_CELL_Liu_hmm_states.May.pdf")
p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
               palette=c("dodgerblue2", "red2","orange2"),
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p3)
p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
               title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
               risk.table.height=0.3)
print(p2)
dev.off()    

res.cox <- coxph(Surv(PFI.time, PFI) ~ strata(biological_states)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv_coxph)
waldp<-summary(res.cox)$waldtest[3]

TCGA_surv_coxph<-TCGA_surv[-which(TCGA_surv$tumor_stage%in%"not_reported"),]

explanatory = c("biological_states","age_at_index", "gender", "bmi", "tumor_stage")
dependent = "Surv(PFI.time, PFI)"
TCGA_surv_coxph %>%
  finalfit.coxph(dependent, explanatory) -> t6
knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph PFI EPI, pEMT, MES",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"PFI_CELL_Liu_hmm_states.coxph.May.pdf")

# 
# ALL MES vs EPI
# # 
# 
# sfit <- survfit(Surv(as.numeric(TCGA_surv_save$OS.time)/365, event = TCGA_surv_save$OS) ~ TCGA_surv_save$biological_states2,data=TCGA_surv_save)
# 
# setwd("/home/guidantoniomt/pseudospace/survival_analysis")
# 
# pdf("OS_CELL_Liu_hmm_states_ALLMES_vs_EPI.May.pdf")
# p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                legend.labs=c("all_mes","epi"), legend.title="biological_states",
#                palette=c( "red2","dodgerblue2"),
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p3)
# p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p2)
# dev.off()    
# 
# res.cox <- coxph(Surv(OS.time, OS) ~ strata(biological_states2)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv)
# waldp<-summary(res.cox)$waldtest[3]
# 
# explanatory = c("biological_states2","age_at_index", "gender", "bmi", "tumor_stage")
# dependent = "Surv(OS.time, OS)"
# TCGA_surv %>%
#   finalfit.coxph(dependent, explanatory) -> t6
# knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph OS ALL MES vs EPI",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"OS_CELL_Liu_hmm_states_ALLMES_vs_EPI.coxph.May.pdf")
# 
# 
# sfit <- survfit(Surv(as.numeric(TCGA_surv_save$DFI.time)/365, event = TCGA_surv_save$DFI) ~ TCGA_surv_save$biological_states2,data=TCGA_surv_save)
# 
# pdf("DFI_CELL_Liu_hmm_states_ALLMES_vs_EPI.May.pdf")
# p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                legend.labs=c("all_mes","epi"), legend.title="biological_states",
#                palette=c( "red2","dodgerblue2"),
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p3)
# p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p2)
# dev.off()  
# 
# res.cox <- coxph(Surv(DFI.time, DFI) ~ strata(biological_states2)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv)
# waldp<-summary(res.cox)$waldtest[3]
# 
# explanatory = c("biological_states2","age_at_index", "gender", "bmi", "tumor_stage")
# dependent = "Surv(DFI.time, DFI)"
# TCGA_surv %>%
#   finalfit.coxph(dependent, explanatory) -> t6
# knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph DFI ALL MES vs EPI",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"DFI_CELL_Liu_hmm_states_ALLMES_vs_EPI.coxph.May.pdf")
# 
# sfit <- survfit(Surv(as.numeric(TCGA_surv_save$DSS.time)/365, event = TCGA_surv_save$DSS) ~ TCGA_surv_save$biological_states2,data=TCGA_surv_save)
# 
# pdf("DSS_CELL_Liu_hmm_states_ALLMES_vs_EPI.May.pdf")
# p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                legend.labs=c("all_mes","epi"), legend.title="biological_states",
#                palette=c( "red2","dodgerblue2"),
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p3)
# p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p2)
# dev.off()    
# 
# res.cox <- coxph(Surv(DSS.time, DSS) ~ strata(biological_states2)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv)
# waldp<-summary(res.cox)$waldtest[3]
# 
# explanatory = c("biological_states2","age_at_index", "gender", "bmi", "tumor_stage")
# dependent = "Surv(DSS.time, DSS)"
# TCGA_surv %>%
#   finalfit.coxph(dependent, explanatory) -> t6
# knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph DFF ALL MES vs EPI",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"DSS_CELL_Liu_hmm_states_ALLMES_vs_EPI.coxph.May.pdf")
# 
# 
# sfit <- survfit(Surv(as.numeric(TCGA_surv_save$PFI.time)/365, event = TCGA_surv_save$PFI) ~ TCGA_surv_save$biological_states2,data=TCGA_surv_save)
# 
# pdf("PFI_CELL_Liu_hmm_states_ALLMES_vs_EPI.May.pdf")
# p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                legend.labs=c("all_mes","epi"), legend.title="biological_states",
#                palette=c( "red2","dodgerblue2"),
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p3)
# p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
#                risk.table.height=0.3)
# print(p2)
# dev.off()   
# 
# res.cox <- coxph(Surv(PFI.time, PFI) ~ strata(biological_states2)+tumor_stage+age_at_index+gender+bmi, data =  TCGA_surv)
# waldp<-summary(res.cox)$waldtest[3]
# 
# explanatory = c("biological_states2","age_at_index", "gender", "bmi", "tumor_stage")
# dependent = "Surv(PFI.time, PFI)"
# TCGA_surv %>%
#   finalfit.coxph(dependent, explanatory) -> t6
# knitr::kable(t6, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"),caption = paste("Coxph PFI ALL MES vs EPI",paste("pvalue:",waldp)))%>%kable_classic(full_width = F, html_font = "Cambria") %>% kableExtra::save_kable(bs_theme=,"PFI_CELL_Liu_hmm_states_ALLMES_vs_EPI.coxph.May.pdf")
# 
# # 
# # look at the Disease Free Interval tumor by tumor
# # 
# 
# # 
# # Step1: check which are the drugs with problem of drugs nomenclaature
# # 
# drugs_tab<-data.frame(table(TCGA_surv$drug_name))
# 
# drugs_tab<-drugs_tab[order(drugs_tab[,2],decreasing=T),]
# 
# # 
# # Correct fluo uracil 
# # 
# 
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="^Fluorouracil$")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5- FU")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-FU")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-fu")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-Flourouracil")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5 FU")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5FU")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5 fluorouracil")]<-"5-Fluorouracil"
# 
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="Fluorouracil (5-Fluorouracil)")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-5-Fluorouracil")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5 5-Fluorouracil")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5 fluorouracil")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-Fluorouracil?")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-Fluoruoracil")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="5-flurouracil")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="Fluorouracilum")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="^Fluorouracil (5-Fluorouracil)$")]<-"5-Fluorouracil"
# TCGA_surv$drug_name[grep(TCGA_surv$drug_name,pattern="^5-fluorouracil$")]<-"5-Fluorouracil"
# 
# 
# # 
# # make uniform thetop ~50 drugs
# # 
# 
# TCGA_surv[grep(TCGA_surv$drug_name,pattern="carboplatin",value=F),"drug_name"]<-"carboplatinum"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Cisplatin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:3,5:8)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"cisplatin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Taxol",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:2)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Taxol"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Paclitaxel",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:2,4,5)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Paclitaxel"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Gemcitabine",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:6,8,9)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Gemcitabine"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Cytoxan",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:3)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Cytoxan"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Oxaliplatin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:2,5,7)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Oxaliplatin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Doxorubicin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:12,18,20,21)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Doxorubicin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Taxotere",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:4)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Taxotere"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Cyclophosphamide",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:3)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Cyclophosphamide"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Leucovorin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:5,7:8)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Leucovorin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Adriamycin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:3)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Adriamycin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Arimidex",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Arimidex"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Docetaxel",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:4)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Docetaxel"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Doxil",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:2)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Doxil"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="gemzar",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Gemcitabine"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Topotecan",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[2:5]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Topotecan"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Bevacizumab",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1,5,6,9,10)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Bevacizumab"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Anastrozole",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Anastrozole"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Epirubicin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Epirubicin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Avastin",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1,5)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Avastin"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Xeloda",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Xeloda"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Folinic acid",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:2,4)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Folinic acid"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Capecitabine",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Capecitabine"
# 
# drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Irinotecan",value=T,fixed=F,ignore.case=T))
# drug_name_to_change<-drug_name[c(1:4,6,7)]
# TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Irinotecan"
# 
# drugs_tab3<-data.frame(table(TCGA_surv$drug_name))
# 
# drugs_tab3<-drugs_tab3[order(drugs_tab3[,2],decreasing=T),]
# 
# drugs_to_use<-grep(drugs_tab3[drugs_tab3[,2]>50,1],pattern="Available",invert=T,value=T)
# 
# library(cowplot)
# 
# # 
# # DFI Cell paper drugs
# # 
# 
# pdf("DFI_drugs.May.pdf",width=15)
# 
# for(i in drugs_to_use){
#   
#   print(i)
#   
#   sub_surv<-TCGA_surv[TCGA_surv$drug_name%in%i,]
#   sfit <- survfit(Surv(as.numeric(sub_surv$DFI.time)/365, event = sub_surv$DFI) ~ sub_surv$biological_states,data=sub_surv)
#   
#   idx<-which(table(sub_surv$biological_states)<=101)
#   
#   if(length(idx)==0){
#     
#   p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                  legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
#                  palette=c("dodgerblue2", "red2","orange2"),
#                  risk.table.height=0.3, ylab="DFI")
#   }else{
#     
#     p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                    legend.title="biological_states",
#                    risk.table.height=0.3, ylab="DFI")
#     
#   }
#   
#   sfit <- survfit(Surv(as.numeric(sub_surv$DFI.time)/365, event = sub_surv$DFI) ~ sub_surv$biological_states2,data=sub_surv)
#   
#   p4<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                  legend.labs=c("all_mes","epi"), legend.title="biological_states",
#                  palette=c( "red2","dodgerblue2"),
#                  risk.table.height=0.3, ylab="DFI")
# 
#   splots <- list()
#   splots[[1]]<-p3
#   splots[[2]]<-p4
#     
# 
#   arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 1, risk.table.height = 0.4,title=paste("DFI:",i))
# }
# 
# dev.off()
# 
# # 
# # DFI Cell paper drugs check
# #  
# 
# pdf("DFI_drugs_check.May.pdf",width=15)
# 
# for(i in drugs_to_use){
#   
#   print(i)
#   
#   sub_surv<-TCGA_surv[TCGA_surv$drug_name%in%i,]
#   sfit <- survfit(Surv(as.numeric(sub_surv$DFI.time)/365, event = sub_surv$DFI) ~ sub_surv$biological_states,data=sub_surv)
#   
#   idx<-which(table(sub_surv$biological_states)<=101)
#   
#   if(length(idx)==0){
#     
#     p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                    legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
#                    palette=c("dodgerblue2", "red2","orange2"),
#                    risk.table.height=0.3,
#                    ylab="DFI")
#   }else{
#     
#     p3<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
#                    legend.title="biological_states",
#                    risk.table.height=0.3)
#     
#   }
#   
#   sfit <- survfit(Surv(as.numeric(sub_surv$DFI.time)/365, event = sub_surv$DFI) ~ sub_surv$biological_states2,data=sub_surv)
#   
#   p4<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,legend.title="biological_states",
#                  risk.table.height=0.3)
#   
#   splots <- list()
#   splots[[1]]<-p3
#   splots[[2]]<-p4
#   
#   
#   arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 1, risk.table.height = 0.4,title=paste("DFI:",i))
# }
# 
# dev.off()
# 
# # 
# # DFI Custom
# # TIME: The event time is the shortest time from initial diagnosis date to the date of an event
# #
# 
# customDFItime<-function(X=TCGA_surv,var_to_use="Complete Remission/Response"){
#                 
#               #for the event is the new_tumour_event_dx_days_to
#   
#               X$time_custom_DFI<-rep("NA",nrow(X))
#               
#               X[X$treatment_outcome_first_course%in%var_to_use,"time_custom_DFI"]<-X[X$treatment_outcome_first_course%in%var_to_use,"new_tumor_event_dx_days_to"]
#               
#               #Treat censored data 
#               X[X$treatment_outcome_first_course%in%"[Not Available]","time_custom_DFI"]<-X[X$treatment_outcome_first_course%in%"[Not Available]","last_contact_days_to"]
#               
#               X[X$treatment_outcome_first_course%in%"[Discrepancy]","time_custom_DFI"]<-X[X$treatment_outcome_first_course%in%"[Discrepancy]","last_contact_days_to"]
#               
#               X[X$treatment_outcome_first_course%in%"[Not Applicable]","time_custom_DFI"]<-X[X$treatment_outcome_first_course%in%"[Not Applicable]","last_contact_days_to"]
#               
#               X[X$treatment_outcome_first_course%in%"[Not Evaluated]","time_custom_DFI"]<-X[X$treatment_outcome_first_course%in%"[Not Evaluated]","last_contact_days_to"]
#               
#               X[X$treatment_outcome_first_course%in%"[Unknown]","time_custom_DFI"]<-X[X$treatment_outcome_first_course%in%"[Unknown]","last_contact_days_to"]
#               
#               X$DFI_status_complete<-rep(0,nrow(X))
#               X[X$treatment_outcome_first_course%in%"Complete Remission/Response","DFI_status_complete"]<-1
#               
#               X$DFI_status_partial<-rep(0,nrow(X))
#               X[X$treatment_outcome_first_course%in%"Partial Remission/Response","DFI_status_partial"]<-1
#               
#               X$DFI_status_persistent<-rep(0,nrow(X))
#               X[X$treatment_outcome_first_course%in%"Persistent Disease","DFI_status_persistent"]<-1
#               
#               X$DFI_status_progressive<-rep(0,nrow(X))
#               X[X$treatment_outcome_first_course%in%"Progressive Disease","DFI_status_progressive"]<-1
#               
#               X$DFI_status_stable<-rep(0,nrow(X))
#               X[X$treatment_outcome_first_course%in%"Stable Disease","DFI_status_stable"]<-1
#               
#               return(X)
# }
# 
# df_complete_remission<-customDFItime(X=TCGA_surv,var_to_use="Complete Remission/Response")
# sfit <- survfit(Surv(as.numeric(df_complete_remission$time_custom_DFI)/365, event = df_complete_remission$DFI_status_complete) ~ df_complete_remission$biological_states,data=df_complete_remission)
# pcomplete<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,risk.table.height=0.3,title="Complete remission",ylim=c(0.4,1))
# 
# df_partial<-customDFItime(X=TCGA_surv,var_to_use="Partial Remission/Response")
# sfit <- survfit(Surv(as.numeric(df_partial$time_custom_DFI)/365, event = df_partial$DFI_status_partial) ~ df_partial$biological_states,data=df_partial)
# ppartial<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,risk.table.height=0.3,title="Partial remission",ylim=c(0.9,1))
# 
# df_persistent<-customDFItime(X=TCGA_surv,var_to_use="Persistent Disease")
# sfit <- survfit(Surv(as.numeric(df_persistent$time_custom_DFI)/365, event = df_persistent$DFI_status_persistent) ~ df_persistent$biological_states,data=df_persistent)
# ppersistent<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,risk.table.height=0.3, title="Persistent",ylim=c(0.9,1))
# 
# df_progressive<-customDFItime(X=TCGA_surv,var_to_use="Progressive Disease")
# sfit <- survfit(Surv(as.numeric(df_progressive$time_custom_DFI)/365, event = df_progressive$DFI_status_progressive) ~ df_progressive$biological_states,data=df_progressive)
# pprogressive<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,risk.table.height=0.3,title="Progressive",ylim=c(0.6,1))
# 
# df_stable<-customDFItime(X=TCGA_surv,var_to_use="Stable Disease")
# sfit <- survfit(Surv(as.numeric(df_stable$time_custom_DFI)/365, event = df_stable$DFI_status_stable) ~ df_stable$biological_states,data=df_stable)
# pstable<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,risk.table.height=0.3,title="Stable",ylim=c(0.6,1))
# 
# 
# splots <- list()
# splots[[1]]<-pcomplete
# splots[[2]]<-ppartial
# splots[[3]]<-ppersistent
# splots[[4]]<-pprogressive
# splots[[5]]<-pstable
# 
# 
# pdf("custom_DFI.May.pdf",width=15,height=15)
# arrange_ggsurvplots(splots, print = TRUE, ncol = 3, nrow = 2, risk.table.height = 0.4)
# dev.off()
# 
