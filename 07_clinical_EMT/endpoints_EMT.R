library(TCGAbiolinks)
library("survival")
library("survminer")
source("ggforest2.R")
library("readxl")
library(plyr)

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

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

setwd(output_dir)

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
