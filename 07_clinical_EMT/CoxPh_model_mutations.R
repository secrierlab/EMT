library(TCGAbiolinks)
library(VennDetail)
library(data.table)
library(survival)
library(survminer)
library(readxl)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

tab_hmm<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)
list_tissue<-paste("TCGA-",unique(sapply(strsplit(tab_hmm[,1],split="\\."),"[[",1)),sep="")
tab_hmm[,1]<-unlist(lapply(strsplit(tab_hmm[,1],split="\\."),FUN=function(x){paste(x[2:4],collapse="-")}))

#
# Get end-points
#

my_data <- read_excel("1-s2.0-S0092867418302290-mmc1.xlsx",sheet="TCGA-CDR")
my_data2 <- my_data[,-1]
TCGA_surv<-merge(my_data2,tab_hmm[,c(1,5)],by.x="bcr_patient_barcode",by.y="samples")

#
# Get clinical data from TCGA
#
TCGA_global_tcga<-vector(mode="list",length(list_tissue))

for(tissue in 1:length(list_tissue)){
  
  print(tissue)
  
  clin.tissue<-GDCquery_clinic(list_tissue[tissue], "clinical")
  
  clin.tissue2<-clin.tissue[,which(colnames(clin.tissue)%in%c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up","tumor_stage","synchronous_malignancy","age_at_index","gender","cigarettes_per_day","alcohol_history","bmi"))]
  
  # days_to_follow_up: Number of days between the date used for index and the date of the patient's last follow-up appointment or contact.
  # 
  
  TCGA_global_tcga[[tissue]]<-data.frame(clin.tissue2,tissue=list_tissue[tissue])
  
}


TCGA_global_tcga2<-do.call(rbind,TCGA_global_tcga)

TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="i/ii nos",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="not reported",value=F),"tumor_stage"]<-"not_reported"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="stage 0",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage i$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage ia$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage ib$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage ii$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iia$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iib$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iic$",value=F),"tumor_stage"]<-"Low"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iii$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iiia$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iiib$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iiic$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iv$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage iva$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage ivb$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage ivc$",value=F),"tumor_stage"]<-"High"
TCGA_global_tcga2[grep(TCGA_global_tcga2$tumor_stage,pattern="^stage x$",value=F),"tumor_stage"]<-"not_reported" # of these samples we don't know nothing they should be undetermined, i have called them not_reported

TCGA_global_tcga2<-TCGA_global_tcga2[-which(TCGA_global_tcga2$tumor_stage%in%"not_reported"),]
TCGA_global_tcga2$tumor_stage<-ifelse(TCGA_global_tcga2$tumor_stage=="High",1,0)

list_rdata_lasso<-c("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData","HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData","HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

list_class1<-c("mes","mes","mix")
list_class2<-c("epi","mix","epi")

markers_lasso<-vector(mode="list",3)

for(i in 1:length(list_rdata_lasso)){

        print(i)

        class1<-list_class1[i]
        class2<-list_class2[i]

        load(list_rdata_lasso[i])
        # get the most frequent (80%) genes from lasso
        df_coef2<-do.call(rbind,res_lasso)

        df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

        features<-data.frame(table(df_coef2[,1]))

        features_freq_df<-features[order(features[,2],decreasing=T),]

        lasso_features_80<-features_freq_df[features_freq_df[,2]>=800,1]
	
	markers_lasso[[i]]<-lasso_features_80
}

names(markers_lasso)<-paste(list_class1,list_class2,sep="_vs_")

setwd(input_dir)

ven <- venndetail(list(mes_vs_epi=markers_lasso[["mes_vs_epi"]],
mes_vs_mix=markers_lasso[["mes_vs_mix"]],
mix_vs_epi=markers_lasso[["mix_vs_epi"]]))

mes_vs_epi_spec<-ven@wide[ven@wide$mes_vs_epi==1 & ven@wide$mix_vs_epi==0,1]
mes_vs_mix_spec<-ven@wide[ven@wide$mes_vs_epi==0 & ven@wide$mix_vs_epi==0,1]
mix_vs_epi_spec<-ven@wide[ven@wide$mes_vs_epi==0 & ven@wide$mix_vs_epi==1,1]

input_ml<-fread("input_for_ml_hmm_states_3_mock_as_gd.txt",data.table=F)
input_ml[,1]<-unlist(lapply(strsplit(input_ml[,1],split="\\-"),FUN=function(x){paste(x[1:3],collapse="-")}))

list_comparisons<-list(mes_vs_epi=grep(mes_vs_epi_spec,pattern="dndscv",value=T),
  mes_vs_mix=grep(mes_vs_mix_spec,pattern="dndscv",value=T),
  mix_vs_epi=grep(mix_vs_epi_spec,pattern="dndscv",value=T))

getFeatures<-function(list_markers,string){

	stat2<-NULL

	for(X in list_markers){

	print(X)

	x2<-as.character(X)

	idx1<-which(colnames(input_ml)%in%"biological_states")
	idx2<-which(colnames(input_ml)%in%x2)

	res<-table(input_ml[,c(idx1,idx2)])

	stat2<-c(stat2,ifelse(res[rownames(res)%in%string,2]>0,1,0))

	}
	names(stat2)<-list_markers

	return(stat2)
}

get_mes<-getFeatures(list_comparisons[[1]],"mes")
markers_mes<-get_mes[get_mes%in%1]

get_mes2<-getFeatures(list_comparisons[[2]],"mes")
markers_mes2<-get_mes2[get_mes2%in%1]

get_mix<-getFeatures(list_comparisons[[3]],"mix")
markers_mix<-get_mix[get_mix%in%1]


list_comparisons2<-list(mes_vs_epi=markers_mes,
                       mes_vs_mix=markers_mes2,
                       mix_vs_epi=markers_mix)


list_comparisons2<-list(mes_vs_epi=markers_mes,
                        mix_vs_epi=markers_mix)

list_endpoints<-c("OS","PFI")

for(end_points in list_endpoints){
  
statistic_coxph2<-data.frame()

# Create new groups
for(am in 1:length(list_comparisons2)){
  
  
  markers2=list_comparisons2[[am]]
  output_comparison=names(list_comparisons2)[[am]]
  
  print(output_comparison)
  
  class1<-sapply(strsplit(output_comparison,split="_vs_"),"[[",1)
  class2<-sapply(strsplit(output_comparison,split="_vs_"),"[[",2)
  
  statistic_coxph<-data.frame()
  
  setwd(output_dir)

  for(cl in 1:length(markers2)){
    
    current_marker<-names(markers2)[cl]
    
    sub_input<-input_ml[,which(colnames(input_ml)%in%c("patients","biological_states",current_marker))]
    
    sub_input2<-sub_input[sub_input$biological_state%in%c(class1,class2),]
    
    TCGA_surv2<-merge(TCGA_surv,sub_input2[,-2],by.x="bcr_patient_barcode",by.y="patients",all.x=F,all.y=F)
    TCGA_surv2_temp<-merge(TCGA_surv2[,-which(colnames(TCGA_surv2)%in%"gender")],TCGA_global_tcga2,by.x="bcr_patient_barcode",by.y="submitter_id")
    
    variable_time<-paste(end_points,"time",sep=".")
    variable_event<-end_points
    
    TCGA_surv3<-TCGA_surv2_temp[,colnames(TCGA_surv2_temp)%in%c(variable_event,variable_time,"biological_states",current_marker,"tumor_stage","age_at_index","gender")]
    TCGA_surv3$gender<-ifelse(TCGA_surv3$gender=="male",1,0)
    TCGA_surv3<-TCGA_surv3[TCGA_surv3$biological_states%in%c(class1,class2),]
    
    colnames(TCGA_surv3)[4]<-"marker"
    
    TCGA_surv3$marker<-ifelse(TCGA_surv3$marker,"yes","no")
    
    if(length(table(TCGA_surv3$marker))!=1){
    
    variable_time<-paste(end_points,"time",sep=".")
    variable_event<-end_points
      
    xsurv <- coxph(Surv(as.numeric(TCGA_surv3[,variable_time])/365, event = TCGA_surv3[,variable_event]) ~ TCGA_surv3$biological_states+marker+tumor_stage+age_at_index+gender, data = TCGA_surv3)
    x <- summary(xsurv)
    
    p.value<-signif(x$wald["pvalue"], digits=2)
    
    # wald.test<-signif(x$wald["test"], digits=2)
    # beta<-signif(x$coef[1], digits=2);#coeficient beta
    # For each treatment group, exp(coef) gives the hazard against whatever you're using as the control or baseline group.

    HR <-signif(x$coef["markeryes",2], digits=2);#exp(beta)
    HR.confint.lower <- signif(x$conf.int["markeryes","lower .95"],2)
    HR.confint.upper <- signif(x$conf.int["markeryes","upper .95"],2)
    
    temp_coxph<-data.frame(comparison=paste(class1,class2,sep="_"),
                           markers=current_marker,
                           p.value=p.value,
                           HR_marker=HR,
                           HR.lower_marker=HR.confint.lower,
                           HR.upper_marker= HR.confint.upper)
    statistic_coxph<-rbind(statistic_coxph,temp_coxph)
    
    }
    
  }
  statistic_coxph<-cbind(statistic_coxph,padjust=p.adjust(statistic_coxph$p.value,"BH"))
  
  statistic_coxph2<-rbind(statistic_coxph2,statistic_coxph)
  
  }
  
    setwd(output_dir)
    library(forcats)
    library(ggforce)
    library(wesanderson)
    library("gtools")
    
    statistic_coxph3<-statistic_coxph2[which(statistic_coxph2$padjust<=0.05),]
    
    statistic_coxph3$Sig = NA
    statistic_coxph3$Sig<-ifelse(statistic_coxph3$HR_marker>1,"bad","good")
    statistic_coxph3$stars<-stars.pval(statistic_coxph3$padjust)
    statistic_coxph3$padjust<- -log(statistic_coxph3$padjust,10)
    
    
    pdf(paste(end_points,"copxh_dotplot_mutations.pdf",sep="."))
    
    for(group in unique(statistic_coxph2$comparison)){
      
      print(group)
      
      stat_group2<-statistic_coxph3[statistic_coxph3$comparison%in%group,]
      stat_group2_bad<-stat_group2[stat_group2$HR_marker > 1.5,]
      
      stat_group2_good<-stat_group2[stat_group2$HR_marker <1,]
      stat_group2_good2<-stat_group2_good[stat_group2_good$HR_marker > 0.7,]
      
      stats_to_use<-rbind(stat_group2_bad,stat_group2_good2)
      
      # stats_to_use$HR_marker<-log(stats_to_use$HR_marker,10)
      # stats_to_use$HR.lower_marker<-log(stats_to_use$HR.lower_marker,10)
      # stats_to_use$HR.upper_marker<-log(stats_to_use$HR.upper_marker,10)
 
      write.table(stats_to_use,file=paste(end_points,group,"copxh_dotplot_mutations.txt",sep="."),sep="\t",row.names=F,quote=F)

      P6 = ggplot(data=stats_to_use, aes(x=fct_reorder(markers,HR_marker,.desc =FALSE), y=HR_marker, ymin=HR.lower_marker, ymax=HR.upper_marker)) +
        geom_errorbar(width=0.2,size=0.1)+geom_point(aes(colour=factor(Sig),size=padjust)) +   scale_colour_manual( values = c("bad" = "#F21A00","good" = "#3B9AB2")) +
        geom_hline(yintercept=1, lty=2,color="blue",size=0.1)+   # add a dotted line at x=1 after flip
        xlab("Markers") + ylab(paste("log10 Hazard Ratio", end_points,"(95% CI)")) +
        theme_classic()+theme(axis.text.x=element_text(angle= 50,size=8,vjust=1,hjust=1))+ ggtitle(group) 
      
      print(P6)
      
    }
    dev.off()
}
