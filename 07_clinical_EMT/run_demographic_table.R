library(TCGAbiolinks)
library(plyr)

setwd("/home/guidantoniomt/pseudospace/HMM")

tab_hmm<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)

list_tissue<-paste("TCGA-",unique(sapply(strsplit(tab_hmm[,1],split="\\."),"[[",1)),sep="")
tab_hmm[,1]<-unlist(lapply(strsplit(tab_hmm[,1],split="\\."),FUN=function(x){paste(x[2:4],collapse="-")}))


TCGA_global_tcga<-vector(mode="list",length(list_tissue))

for(tissue in 1:length(list_tissue)){

  print(tissue)

  clin.tissue<-GDCquery_clinic(list_tissue[tissue], "clinical")

  clin.tissue2<-clin.tissue[,which(colnames(clin.tissue)%in%c("submitter_id","ajcc_pathologic_m", "tumor_stage","disease"))]

  # days_to_follow_up: Number of days between the date used for index and the date of the patient's last follow-up appointment or contact.
  # 

  TCGA_global_tcga[[tissue]]<-clin.tissue2

}


TCGA_global_tcga2<-do.call("rbind.fill",TCGA_global_tcga)
TCGA_surv<-merge(TCGA_global_tcga2,tab_hmm[,c(1,5)],by.x="submitter_id",by.y="samples")
TCGA_surv$custom_stage<-TCGA_surv$tumor_stage
#TCGA_surv$custom_stage[TCGA_surv$ajcc_pathologic_m=="M1"]<-"Metastatic"
TCGA_surv$custom_stage2<-TCGA_surv$custom_stage

df_stages<-table(TCGA_surv$tumor_stage,TCGA_surv$biological_states)

write.table(df_stages,file="stages_hmm.txt",sep="\t",quote=F)

low_stage_cancers<-apply(df_stages[c(1,3:10),],2,sum)
high_stage_cancers<-apply(df_stages[c(11:nrow(df_stages)),],2,sum)
not_reported<-df_stages[2,]



df_stage<-data.frame(high_stage=high_stage_cancers,
	   low_stage=low_stage_cancers,
	   not_reported=not_reported)

df_stage<-df_stage[c(2,3,1),]

write.table(df_stage,file="stages_hmm_collapsed2.txt",sep="\t",quote=F)


#chisq.test using epi, mes, mix, and high,low, not reported
stat_mes_epi_mix_all_stages<-chisq.test(df_stage)

setwd("/home/guidantoniomt/pseudospace/survival_analysis")

library(corrplot)
pdf("chisq_all_groups_and_stages.pdf")
corrplot(stat_mes_epi_mix_all_stages$residuals, is.cor = FALSE)
contrib <- 100*stat_mes_epi_mix_all_stages$residuals^2/stat_mes_epi_mix_all_stages$statistic
corrplot(contrib, is.cor = FALSE)
dev.off()

#chisq.test using 1)EPI 2)MES 3)pEMT and 1) low stage 2) high stage
chisq2<-chisq.test(df_stage[,-3])

pdf("chisq_all_groups_and_low_high.pdf")
corrplot(chisq2$residuals, is.cor = FALSE)
contrib <- 100*chisq2$residuals^2/chisq2$statistic
corrplot(contrib, is.cor = FALSE)
dev.off()


setwd("/home/guidantoniomt/pseudospace/survival_analysis")
write.table(df_stage,file="TCGA_EMT_STATES_clinical_stage.txt",sep="\t",row.names=T,quote=F)
