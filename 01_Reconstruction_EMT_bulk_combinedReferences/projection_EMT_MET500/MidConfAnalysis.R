####################################
#Mid confidence pseudotime analysis

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500")
input_cp_final <- read.table("AveragedPseudotime_midconf.txt", header = TRUE, sep = "\t")

# 
# Compute the EMT scores of the samples along the MOCK psuedotime
# 

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
load('input_ge.RData')
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
epi_genes<-intersect(input_ge[,1], markers_genes_read[markers_genes_read$status=="Epithelial_marker",1])
mes_genes<-intersect(input_ge[,1], markers_genes_read[markers_genes_read$status=="Mesenchymal_marker",1])

# Compute an EMT score of the patients
# VIM + CDH2 + FOXC2 + SNAI1 + SNAI2 + TWIST1 + FN1 + ITGB6 + MMP2 + MMP3 + MMP9 + SOX10 + GCS − CDH1 − DSP − OCLN
# https://www.nature.com/articles/s41598-018-21061-1

input_ge2<-t(input_ge[input_ge[,1]%in%c(epi_genes,mes_genes),which(colnames(input_ge)%in%input_cp_final$patients)])

zscore<-function(x){(x-mean(x))/sd(x)}
input_ge3<-t(apply(input_ge2,2,zscore))

all_scores_emt<-NULL

for(iemt in 1:ncol(input_ge3)){
  
  mean_epi<-mean(input_ge3[epi_genes,iemt])
  mean_mes<-mean(input_ge3[mes_genes,iemt])
  
  score_emt<-mean_mes-mean_epi
  
  all_scores_emt<-c(all_scores_emt,score_emt)
  
}

df_emt_samples<-data.frame(ID=colnames(input_ge3),EMT_scores=all_scores_emt)
df_emt_samples$ID<-as.character(df_emt_samples$ID)
input_cp_final$patients<-as.character(input_cp_final$patients)

input_cp_final2<-merge(x=input_cp_final[,-c(4:5)],y=df_emt_samples,by.x="patients",by.y="ID")

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/output_dir")
write.table(input_cp_final2,file='MidConf_pseudotime_withEMT.txt',row.names=F,quote=F,sep='\t')

input_cp_final2$new_ann<-rep("no",nrow(input_cp_final2))
input_cp_final2$new_ann[grep(input_cp_final2$patients,pattern="TCGA")]<-"TCGA"
input_cp_final2$new_ann[grep(input_cp_final2$patients,pattern="TCGA",invert=T)]<-"MET500"
input_cp_final2$new_ann<-as.factor(input_cp_final2$new_ann)


input_cp_final2$emt_status<-ifelse(input_cp_final2$EMT_scores>0,"High_EMT","Low_EMT")
summary(input_cp_final2$mock_pseudospace)
input_cp_final2$pseudospace_status<-ifelse(input_cp_final2$mock_pseudospace<=50,"Late_Pseudotime","Early_Pseudotime")

library(ggplot2)
library(ggridges)
library(corrplot)

pdf(paste("TCGA_MET500_pseudotime_vs_emt_scores_mock_mid_conf.pdf",sep="."),width=12,pointsize=12)
p<-ggplot(input_cp_final2, 
          aes(x=mock_pseudospace, y=EMT_scores,color=new_ann)) +
  geom_point(shape=18,size=2,alpha=0.5)+scale_x_reverse()+ 
  scale_color_manual(values=c('#ee0c0c','#B5B7FF'))+
  geom_smooth(aes(group=new_ann),method = "lm", formula = y ~ poly(x, 10),se=TRUE, linetype="dashed",color="blue")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+theme_classic()
print(p)

input_cp_final3<-input_cp_final2[input_cp_final2$new_ann%in%"MET500",]
chisq <- chisq.test(table(input_cp_final3$pseudospace_status,input_cp_final3$emt_status)[2:1,])
contrib <- 100*chisq$residuals^2/chisq$statistic
corrplot(chisq$residuals, is.cor = FALSE)
corrplot(contrib, is.cor = FALSE)
dev.off()


# 
# Compare the EMT scores of the MET500 samples and the samples defined using HMM
# 

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/02_HMM_macrostases_EMT/output_dir_mid_withzeros")
tab_hmm<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)
samples_with_HMM_and_MET500<-merge(tab_hmm,input_cp_final2,by.x="samples",by.y="patients",all.x=T,all.y=T)
samples_with_HMM_and_MET500$HMM_states[which(is.na(samples_with_HMM_and_MET500$HMM_states))]<-"MET500"
samples_with_HMM_and_MET500$HMM_states<-gsub(samples_with_HMM_and_MET500$HMM_states,pattern="1",replacement="epi")
samples_with_HMM_and_MET500$HMM_states<-gsub(samples_with_HMM_and_MET500$HMM_states,pattern="2",replacement="pEMT")
samples_with_HMM_and_MET500$HMM_states<-gsub(samples_with_HMM_and_MET500$HMM_states,pattern="3",replacement="mes")
samples_with_HMM_and_MET500$HMM_states<-factor(samples_with_HMM_and_MET500$HMM_states,levels=c("epi","pEMT","mes","MET500"))

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/output_dir")

library(ggpubr)

pdf("Compare_EMT_scores_primary_HMMstates_with_MET500_midconf_full.pdf")
p<-ggviolin(samples_with_HMM_and_MET500,
            x = "HMM_states",
            y = "EMT_scores",
            fill = "HMM_states",
            palette  = c("#56b4e9","#E69f00","#EE0C0C","#A10070"),
            add = "boxplot",
            add.params = list(fill = "white"))+
  stat_compare_means()+geom_hline(yintercept=0)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

dev.off()

# 
# Compare the pseudotime values for each category HMM
# 

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/02_HMM_macrostases_EMT/output_dir_mid_withzeros")
tab_hmm_test<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)

#tab_hmm_test$pseudospace<-max(tab_hmm_test$pseudospace)-tab_hmm_test$pseudospace
tab_hmm_test$pseudospace<-1-(tab_hmm_test$pseudospace/max(tab_hmm_test$pseudospace))

# 
# Compare the hypoxia scores of the MET500 samples with the hypoxia scores of the TCGA samples
# 

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/output_dir")

pdf("Pseudotime_values_for_HMMstate_midconf.pdf")
p<-ggviolin(tab_hmm_test,
            x = "biological_states",
            y = "pseudospace",
            fill = "biological_states",
            palette  = c("#56b4e9","#E69f00","#EE0C0C","#A10070"),
            add = c("boxplot"),
            add.params = list(fill = "white"))+
  stat_compare_means()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

