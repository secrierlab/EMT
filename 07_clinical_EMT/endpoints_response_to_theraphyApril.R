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

TCGA_surv<-merge(TCGA_surv,TCGA_global_drugs,by.x="bcr_patient_barcode",by.y="bcr_patient_barcode",all.y=F)

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

TCGA_surv$biological_states<-factor(TCGA_surv$biological_states,levels=c("epi","mix","mes"))
TCGA_surv$biological_states2<-factor(TCGA_surv$biological_states2,levels=c("epi","all_mes"))

# 
# make uniform thetop ~50 drugs
# 

TCGA_surv[grep(TCGA_surv$drug_name,pattern="carboplatin",value=F),"drug_name"]<-"carboplatinum"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Cisplatin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:3,5:8)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"cisplatin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Taxol",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:2)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Taxol"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Paclitaxel",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:2,4,5)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Paclitaxel"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Gemcitabine",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:6,8,9)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Gemcitabine"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Cytoxan",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:3)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Cytoxan"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Oxaliplatin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:2,5,7)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Oxaliplatin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Doxorubicin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:12,18,20,21)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Doxorubicin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Taxotere",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:4)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Taxotere"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Cyclophosphamide",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:3)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Cyclophosphamide"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Leucovorin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:5,7:8)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Leucovorin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Adriamycin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:3)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Adriamycin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Arimidex",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Arimidex"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Docetaxel",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:4)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Docetaxel"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Doxil",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:2)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Doxil"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="gemzar",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Gemcitabine"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Topotecan",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[2:5]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Topotecan"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Bevacizumab",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1,5,6,9,10)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Bevacizumab"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Anastrozole",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Anastrozole"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Epirubicin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Epirubicin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Avastin",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1,5)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Avastin"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Xeloda",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Xeloda"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Folinic acid",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:2,4)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Folinic acid"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Capecitabine",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Capecitabine"

drug_name<-unique(grep(TCGA_surv$drug_name,pattern="Irinotecan",value=T,fixed=F,ignore.case=T))
drug_name_to_change<-drug_name[c(1:4,6,7)]
TCGA_surv[TCGA_surv$drug_name%in%drug_name_to_change,"drug_name"]<-"Irinotecan"

drugs_tab3<-data.frame(table(TCGA_surv$drug_name))

drugs_tab3<-drugs_tab3[order(drugs_tab3[,2],decreasing=T),]

drugs_to_use<-grep(drugs_tab3[drugs_tab3[,2]>50,1],pattern="Available",invert=T,value=T)

TCGA_surv2<-TCGA_surv[TCGA_surv$drug_name%in%drugs_to_use,]

TCGA_surv2[TCGA_surv2$treatment_outcome_first_course %in% "Complete Remission/Response","treatment_outcome_first_course"]<- "response"
TCGA_surv2[TCGA_surv2$treatment_outcome_first_course %in% "Partial Remission/Response","treatment_outcome_first_course"]<- "response"
TCGA_surv2[TCGA_surv2$treatment_outcome_first_course %in% "Persistent Disease","treatment_outcome_first_course"]<- "not_response"
TCGA_surv2[TCGA_surv2$treatment_outcome_first_course %in% "Progressive Disease","treatment_outcome_first_course"]<- "not_response"
TCGA_surv2[TCGA_surv2$treatment_outcome_first_course %in% "Stable Disease","treatment_outcome_first_course"]<- "not_response"

TCGA_surv3<-TCGA_surv2[TCGA_surv2$treatment_outcome_first_course %in% c("response","not_response"),]

list_comparisons<-vector(mode="list",2)

list_comparisons[[1]]<-c("mes","epi")
list_comparisons[[2]]<-c("mix","epi")
list_comparisons[[3]]<-c("all_mes","epi")
#
# Find the significant drugs related with the EMT and response
#

drugs<-NULL
pvalue<-NULL
comparison<-NULL
OR<-NULL

for(i in 1:length(list_comparisons)){
  
  print(i)
  
  current_comparison<-list_comparisons[[i]]
  
  print(current_comparison)
  
  if(i!=3){

  TCGA_surv_classes<-TCGA_surv3[TCGA_surv3$biological_states%in%current_comparison,]
  
  TCGA_surv_classes$biological_states<-as.character(TCGA_surv_classes$biological_states)

  }else{
     
  TCGA_surv_classes<-TCGA_surv3[TCGA_surv3$biological_states2%in%current_comparison,]
  
  TCGA_surv_classes$biological_states2<-as.character(TCGA_surv_classes$biological_states2)
	
  }

  # start to iterate for each drug
  list_drugs<-unique(TCGA_surv_classes$drug_name)
  
  for(d in list_drugs){
    
  if(i!=3){	  
   
  TCGA_surv_drug<-TCGA_surv_classes[TCGA_surv_classes$drug_name%in%d,colnames(TCGA_surv_classes)%in%c("biological_states","treatment_outcome_first_course")]
  
  }else{
  
  TCGA_surv_drug<-TCGA_surv_classes[TCGA_surv_classes$drug_name%in%d,colnames(TCGA_surv_classes)%in%c("biological_states2","treatment_outcome_first_course")]

  }

  class1<-current_comparison[1]
  class2<-current_comparison[2]
  
  contingency_drug<-t(table(TCGA_surv_drug))
  
  if(nrow(contingency_drug)==2 & ncol(contingency_drug)==2){
    
  contingency_drug_reorder<-contingency_drug[c(class1,class2),c("response","not_response")]
  
  fishres<-fisher.test(contingency_drug_reorder)

  # contrib <- 100*chisqres$residuals^2/chisqres$statistic
  # 
  # contrib <- round(contrib, 3)
  # contrib2<- as.numeric(contrib)
  
  comparison<-c(comparison,paste(class1,class2,sep="_"))
  drugs<-c(drugs,d)
  pvalue<-c(pvalue,fishres$p.value)
  OR<-c(OR,fishres$estimate)
  
  }
  
  }
  
}

resDrugs<-data.frame(comparison,drugs,pvalue,OR)

 
