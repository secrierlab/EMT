#
# Upload results pseudotime
#
require(sjPlot)

setwd('/data/pseudospace/HMM')
input_file<-c("HMM_results_nstates_3.txt")
knn_df_tcga<-read.delim(file=input_file)

scores_emt_file<-"HMM_results_nstates_tumors_for_states3.withEMT.txt"
scores_EMT<-read.delim(file=scores_emt_file)

code_for_met<-sapply(strsplit(as.character(scores_EMT[,2]),split="\\."),"[[",5)
tum_for_met<-sapply(strsplit(as.character(scores_EMT[,2]),split="\\."),"[[",1)

# see which are the samples with at least a metastatic samples
table(code_for_met,tum_for_met)

#
# Read mutational signatures data SigProfiler
#

setwd("/data/pseudospace/mutational_signature/sbs_solutions_deconstructSigs_Sigprofiler")

sampleinfo<-c("sigs.defaultnounknown.ACC.RData", "sigs.defaultnounknown.BLCA.RData","sigs.defaultnounknown.BRCA.metastatic.RData",
             "sigs.defaultnounknown.CESC.metastatic.RData","sigs.defaultnounknown.CHOL.RData","sigs.defaultnounknown.COAD.metastatic.RData",
             "sigs.defaultnounknown.ESCA.metastatic.RData","sigs.defaultnounknown.GBM.RData","sigs.defaultnounknown.HNSC.metastatic.RData",
             "sigs.defaultnounknown.KICH.RData","sigs.defaultnounknown.KIRC.RData","sigs.defaultnounknown.KIRP.RData",
             "sigs.defaultnounknown.LGG.RData","sigs.defaultnounknown.LIHC.RData","sigs.defaultnounknown.LUAD.RData",
             "sigs.defaultnounknown.LUSC.RData","sigs.defaultnounknown.MESO.RData","sigs.defaultnounknown.OV.RData",
             "sigs.defaultnounknown.PAAD.metastatic.RData","sigs.defaultnounknown.PRAD.metastatic.RData","sigs.defaultnounknown.READ.RData",
             "sigs.defaultnounknown.SARC.RData","sigs.defaultnounknown.SKCM.metastatic.RData","sigs.defaultnounknown.STAD.RData",
             "sigs.defaultnounknown.TGCT.RData","sigs.defaultnounknown.THCA.RData","sigs.defaultnounknown.THYM.RData",
             "sigs.defaultnounknown.UCEC.RData","sigs.defaultnounknown.UCS.RData")

store_matrices<-vector(mode='list',length=length(sampleinfo))

list_tumours<-NULL

for(i in 1:length(sampleinfo)){

file_to_read<-get(load(sampleinfo[i]))

tumour_name<-gsub(gsub(sampleinfo[i],pattern="sigs.defaultnounknown.",replacement=""),pattern=".RData",replacement="")

tab_input<-data.frame(tumor=as.character(sampleinfo[i]),samples=rownames(file_to_read),file_to_read)

store_matrices[[i]]<-tab_input
list_tumours<-c(list_tumours,tumour_name)

}

names(store_matrices)<-list_tumours
  
#
# ggballoon chart
#
library(reshape2)
library(ggpubr)

sigs_ballon<-melt(store_matrices)[,-1]
colnames(sigs_ballon)[ncol(sigs_ballon)]<-"tum"
sigs_ballon$tum<-as.factor(sigs_ballon$tum)

setwd("/data/pseudospace/mutational_signature/")

avg_contribution<-aggregate(.~tum+variable, sigs_ballon, mean)

sigs_artefacts<-c('SBS27','SBS43','SBS45','SBS47','SBS48','SBS49','SBS50','SBS51','SBS52','SBS53','SBS54','SBS55','SBS56','SBS57','SBS58','SBS59','SBS60')

avg_contribution2<-avg_contribution[-which(avg_contribution[,2]%in%sigs_artefacts),]
avg_contribution2[,1]<- factor(avg_contribution2[,1],rev(levels(avg_contribution2[,1])))
avg_contribution2<-avg_contribution2[avg_contribution2$value>0,]

pdf("sigs_for_cancer_balloonchart.pdf")
p<-ggballoonplot(avg_contribution2,x="variable",y="tum",fill = "value",size="value",size.range=c(0,5))+scale_fill_viridis_c(option = "C")+theme(axis.text.x = element_text(size = 5,angle=90, hjust=1),axis.text.y = element_text(size = 5))
print(p)
dev.off()

library(plyr)
matrix_mutsigs<-do.call(rbind.fill,store_matrices)

#remove artifacts

sigs_artefacts<-c('SBS27','SBS43','SBS45','SBS47','SBS48','SBS49','SBS50','SBS51','SBS52','SBS53','SBS54','SBS55','SBS56','SBS57','SBS58','SBS59','SBS60')

if(length(which(colnames(matrix_mutsigs)%in%sigs_artefacts)) != 0){

matrix_mutsigs<-matrix_mutsigs[,-c(which(colnames(matrix_mutsigs)%in%sigs_artefacts))]

}


#
# Merge pseudotime with mutational signatures
#

knn_df_tcga[,1]<- sapply(strsplit(as.character(knn_df_tcga[,1]),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")})

matrix_mutsigs[,2]<-sapply(strsplit(as.character(matrix_mutsigs[,2]),split="\\-"),FUN=function(X){paste(X[1:4],collapse="-")})

matrix_mutsigs2<-matrix_mutsigs[matrix_mutsigs[,2]%in%knn_df_tcga[,1],]

scores_EMT[,2]<-sapply(strsplit(as.character(scores_EMT[,2]),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")})

pseudospace_with_score<-merge(knn_df_tcga,scores_EMT,by="samples")

mutsig_with_pseudotime<-merge(pseudospace_with_score,matrix_mutsigs2,by.x="samples",by.y="samples")


#
# Create ballonchart for emt group
#
mutsig_with_pseudotime$tumor<-gsub(mutsig_with_pseudotime$tumor,pattern="sigs.defaultnounknown.",replacement="")
mutsig_with_pseudotime$tumor<-gsub(mutsig_with_pseudotime$tumor,pattern=".RData",replacement="")
mutsig_with_pseudotime$states<-factor(mutsig_with_pseudotime$states,levels=c(1,2,3))

mutsig_with_pseudotime_ballonchart<-melt(mutsig_with_pseudotime[,c(8:ncol(mutsig_with_pseudotime))])

avg_contribution2<-aggregate(.~tumor+variable+states, mutsig_with_pseudotime_ballonchart, mean)

avg_contribution2<-avg_contribution2[avg_contribution2$value>0,]

pdf("sigs_for_cancer_EMT_groups_balloonchart.pdf")
p<-ggballoonplot(avg_contribution2,x="variable",y="tumor",fill = "value",size="value",size.range=c(0,5))+scale_fill_viridis_c(option = "C")+theme(axis.text.x = element_text(size = 5,angle=90, hjust=1),axis.text.y = element_text(size = 5))+facet_wrap(.~states)
print(p)
dev.off()

#
# Measure correlation between EMT and MUT sig
#
list_states<-unique(mutsig_with_pseudotime$states)

all_cor<-data.frame()

for(lss in list_states){
  
  est_cor<-cor(mutsig_with_pseudotime[mutsig_with_pseudotime$states%in%lss,c(7,10:ncol(mutsig_with_pseudotime))])[1,] 
  est_cor2<-data.frame(lss,rownames(data.frame(est_cor[-1])),data.frame(est_cor[-1]))
  colnames(est_cor2)<-c("states","sigs","cor_emt_sigs")
  
  all_cor<-rbind(all_cor,est_cor2)
}

pdf("cor_emt_for_sigs_bySTATES.pdf")
all_cor$states<-factor(all_cor$states,levels=c("2","1","3"))
p<-ggballoonplot(all_cor,x="states",y="sigs",fill = "cor_emt_sigs",size="cor_emt_sigs",size.range=c(0,5))+gradient_fill(c("blue", "white", "red"))+theme(axis.text.x = element_text(size = 5,angle=90, hjust=1),axis.text.y = element_text(size = 5))
print(p)
dev.off()


#
# Make a regression between the pseudotime and the mutational signatures
#

mutsig_with_pseudotime2<-mutsig_with_pseudotime[order(mutsig_with_pseudotime$pseudospace,decreasing=T),]

mutsig_with_pseudotime2[,-c(1:9)]<-mutsig_with_pseudotime2[,-c(1:9)]
mutsig_with_pseudotime2[is.na(mutsig_with_pseudotime2)]<-0

count_zero<-function(X){
length(which(X==0))
}


#
# Linear Regression between the mutational signatures and the pseudo time correcting by tissue
#

library(jtools)
library(ggstance)
library(broom.mixed)
library(MASS)
library(ggfortify)
library(caret)

#
# Get only the columns in which the sum is greater than 0 (signals of mutational signatures)
#

# columns_to_save<-c("score_emt","tumors","states",names(which(apply(mutsig_with_pseudotime2[,10:ncol(mutsig_with_pseudotime2)],2,sum) > 0)))
# mutsig_with_pseudotime3<-mutsig_with_pseudotime2[,colnames(mutsig_with_pseudotime2)%in%columns_to_save]
# 
# mutsig_with_pseudotime3$states<-gsub(mutsig_with_pseudotime3$states,pattern="1",replacement="pEMT")
# mutsig_with_pseudotime3$states<-gsub(mutsig_with_pseudotime3$states,pattern="2",replacement="epi")
# mutsig_with_pseudotime3$states<-gsub(mutsig_with_pseudotime3$states,pattern="3",replacement="mes")
# 
# form_lm<-as.formula(paste(paste("score_emt","~",paste(grep(colnames(mutsig_with_pseudotime3),pattern="SBS",value=T),collapse="+")),"+tumors"))
# 

# input_for_boxplot<-mutsig_with_pseudotime2[,c(1,5,6,7,c(10:ncol(mutsig_with_pseudotime2)))]
# sigs_to_use<-names(which(apply(input_for_boxplot[,-c(1:4)],2,sum)>0))
# 
# col_to_use<-c(colnames(input_for_boxplot)[1:4],sigs_to_use)
# 
# input_for_boxplot<-input_for_boxplot[,colnames(input_for_boxplot) %in%col_to_use]
# 
# input_for_boxplot2<-data.frame()
# 
# for(sbs in 5:ncol(input_for_boxplot)){
# 
#         sbs_part<-input_for_boxplot[,c(1,2,3,4,sbs)]
#         sbs_part<-sbs_part[sbs_part[,5]>0.05,]
#         res<-data.frame(sbs=colnames(input_for_boxplot)[sbs],sbs_part)
#         colnames(res)[6]<-"sbs_freq"
# 
#         input_for_boxplot2<-rbind(input_for_boxplot2,res)
# }
# 
# set.seed(123)
# 
# training.samples <- mutsig_with_pseudotime3$score_emt %>% createDataPartition(p = 0.7, list = FALSE)
# 
# train.data  <- mutsig_with_pseudotime3[training.samples, ]
# test.data <- mutsig_with_pseudotime3[-training.samples, ]
# 
# res_lm2<-lm(form_lm,train.data,na.omit=T)
# 
# res_sm<-summary(res_lm2)$coefficients[,4]
# update_sigs<-grep(names(res_sm[res_sm<=0.05]),pattern="SBS",value=T)
# 
# require(ggplot2)

setwd("/data/pseudospace/mutational_signature")
# 
# pdf("coefficients_mutsigs_with_pseudotime_tissuecorrect.pdf")
# p1<-plot_summs(res_lm2, scale = TRUE, inner_ci_level = .9,coefs=update_sigs)
# print(p1)
# p2<-plot_model(res_lm2, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,terms=update_sigs)
# print(p2)
# 
# 
# input_for_boxplot3<-input_for_boxplot2[input_for_boxplot2$sbs%in%update_sigs,]
# 
# p3<-ggplot(input_for_boxplot3, aes(x = sbs_freq)) +
#   stat_ecdf(aes(color = biological_states,linetype = biological_states),
#               geom = "step", size = 0.5) +
#   scale_color_manual(values = c("#56b4e9","#EE0C0C","#E69f00"))+
#   labs(y = "f(contr. mut.sig.)")+facet_wrap(.~sbs)+theme_bw()
# 
# print(p3)
# dev.off()
# 
# form_lm<-as.formula(paste(paste("score_emt","~",paste(grep(colnames(mutsig_with_pseudotime3),pattern="SBS",value=T),collapse="+"))))
# 
# set.seed(123)
# 
# res_lm2_notissue<-lm(form_lm,train.data,na.omit=T)
# 
# res_sm<-summary(res_lm2_notissue)$coefficients[,4]
# update_sigs<-grep(names(res_sm[res_sm<=0.05]),pattern="SBS",value=T)
# 
# pdf("coefficients_mutsigs_with_pseudotime_notissuecorrect.pdf")
# p1<-plot_summs(res_lm2, scale = TRUE, inner_ci_level = .9,coefs=update_sigs)
# print(p1)
# p2<-plot_model(res_lm2, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,terms=update_sigs)
# print(p2)
# input_for_boxplot3<-input_for_boxplot2[input_for_boxplot2$sbs%in%update_sigs,]
# p3<-ggplot(input_for_boxplot3, aes(x = sbs_freq)) +
#   stat_ecdf(aes(color = biological_states,linetype = biological_states),
#               geom = "step", size = 0.5) +
#   scale_color_manual(values =  c("#56b4e9","#EE0C0C","#E69f00"))+
#   labs(y = "f(contr. mut.sig.)")+facet_wrap(.~sbs)+theme_bw()
# 
# print(p3)
# dev.off()


# #
# # glm EPI-MES
# #
# 
# set.seed(123)
# 
# mat_temp<-mutsig_with_pseudotime3[,-2]
# mat_temp<-mat_temp[mat_temp$states%in%c("mes","epi"),]
# mat_temp$states<-as.factor(ifelse(mat_temp$states=="mes",1,0))
# mat_temp$tumors<-as.factor(mat_temp$tumors)
#   
# training.samples <- mat_temp$states %>% createDataPartition(p = 0.7, list = FALSE)
# 
# train.data  <- mat_temp[training.samples, ]
# test.data <- mat_temp[-training.samples, ]
# 
# # tissue effect
# form_glm<-as.formula(paste(paste("states","~",paste(grep(colnames(train.data),pattern="SBS",value=T),collapse="+")),"+tumors"))
# 
# model_mes_epi<-glm(form_glm,train.data,family = "binomial")
# step.model <- model_mes_epi %>% stepAIC(trace = TRUE)
# 
# update_sigs<-grep(names(which(summary(step.model)$coefficients[,4]<=0.05)),pattern="SBS",value=T)
# 
# update_sigs<-grep(names(summary(step.model)$coefficients[,4]),pattern="SBS",value=T)
# 
# form_glm<-as.formula(paste(paste("states","~",paste(update_sigs,collapse="+")),"+tumors"))
# 
# model_mes_epi<-glm(form_glm,train.data,family = "binomial")
# 
# require(sjPlot)
# 
# pdf("glm_sbs_mes_epi.pdf")
# p1<-plot_summs(step.model, scale = TRUE, inner_ci_level = .9,coefs=update_sigs)
# print(p1)
# p2<-plot_model(model_mes_epi, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,terms=update_sigs)
# print(p2)
# dev.off()
# 
# # no tissue effect
# 
# form_glm<-as.formula(paste(paste("states","~",paste(grep(colnames(train.data),pattern="SBS",value=T),collapse="+"))))
# 
# model_mes_epi<-glm(form_glm,train.data,family = "binomial")
# step.model <- model_mes_epi %>% stepAIC(trace = TRUE)
# 
# update_sigs<-grep(names(which(summary(step.model)$coefficients[,4]<=0.05)),pattern="SBS",value=T)
# 
# form_glm<-as.formula(paste(paste("states","~",paste(update_sigs,collapse="+"))))
# 
# model_mes_epi<-glm(form_glm,train.data,family = "binomial")
# 
# require(sjPlot)
# 
# pdf("glm_sbs_mes_epi_notissue.pdf")
# # plot_summs(step.model, scale = TRUE, inner_ci_level = .9,coefs=update_sigs)
# p<-plot_model(model_mes_epi, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,terms=update_sigs)
# print(p)
# dev.off()

# no tissue effect

#
# glm EPI-pEMT
#
# 
# set.seed(123)
# 
# mat_temp<-mutsig_with_pseudotime3[,-2]
# mat_temp<-mat_temp[mat_temp$states%in%c("pEMT","epi"),]
# mat_temp$states<-as.factor(ifelse(mat_temp$states=="pEMT",1,0))
# mat_temp$tumors<-as.factor(mat_temp$tumors)
# 
# training.samples <- mat_temp$states %>% createDataPartition(p = 0.7, list = FALSE)
# 
# train.data  <- mat_temp[training.samples, ]
# test.data <- mat_temp[-training.samples, ]
# 
# # tissue effect
# form_glm<-as.formula(paste(paste("states","~",paste(grep(colnames(train.data),pattern="SBS",value=T),collapse="+")),"+tumors"))
# 
# model_pemt_epi<-glm(form_glm,train.data,family = "binomial")
# step.model <- model_pemt_epi %>% stepAIC(trace = TRUE)
# 
# update_sigs<-grep(names(summary(step.model)$coefficients[,4]),pattern="SBS",value=T)
# 
# form_glm<-as.formula(paste(paste("states","~",paste(update_sigs,collapse="+")),"+tumors"))
# 
# model_pemt_epi<-glm(form_glm,train.data,family = "binomial")
# 
# pdf("glm_sbs_pemt_epi.pdf")
# p1<-plot_summs(step.model, scale = TRUE, inner_ci_level = .9,coefs=update_sigs)
# print(p1)
# p2<-plot_model(model_pemt_epi, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,terms=update_sigs)
# print(p2)
# dev.off()
# 
# # no tissue effect
# form_glm<-as.formula(paste(paste("states","~",paste(grep(colnames(train.data),pattern="SBS",value=T),collapse="+"))))
# 
# model_pemt_epi<-glm(form_glm,train.data,family = "binomial")
# step.model <- model_pemt_epi %>% stepAIC(trace = TRUE)
# 
# update_sigs<-grep(names(which(summary(step.model)$coefficients[,4]<=0.05)),pattern="SBS",value=T)
# 
# form_glm<-as.formula(paste(paste("states","~",paste(update_sigs,collapse="+"))))
# 
# model_pemt_epi<-glm(form_glm,train.data,family = "binomial")
# 
# pdf("glm_sbs_pemt_epi_notissue.pdf")
# p1<-plot_summs(step.model, scale = TRUE, inner_ci_level = .9,coefs=update_sigs)
# print(p1)
# p2<-plot_model(model_pemt_epi, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,terms=update_sigs)
# print(p2)
# dev.off()

#
# Boxplot EMT scores vs mut sigs
#

#
# Filter2: Get only the columns in which the sum is greater than 0 (signals of mutational signatures)
#

input_for_boxplot<-mutsig_with_pseudotime2[,c(1,5,6,7,c(10:ncol(mutsig_with_pseudotime2)))]
sigs_to_use<-names(which(apply(input_for_boxplot[,-c(1:4)],2,sum)>0))

col_to_use<-c(colnames(input_for_boxplot)[1:4],sigs_to_use)

input_for_boxplot<-input_for_boxplot[,colnames(input_for_boxplot) %in%col_to_use]

input_for_boxplot2<-data.frame()

for(sbs in 5:ncol(input_for_boxplot)){
	
	sbs_part<-input_for_boxplot[,c(1,2,3,4,sbs)]
	# sbs_part<-sbs_part[sbs_part[,5]>0.05,]
	res<-data.frame(sbs=colnames(input_for_boxplot)[sbs],sbs_part)
	colnames(res)[6]<-"sbs_freq"

	input_for_boxplot2<-rbind(input_for_boxplot2,res)
}


pdf("emt_vs_sbs_new.pdf",width=15)
p<-ggviolin(input_for_boxplot2, x.text.angle = 90 , "sbs", "score_emt", fill = "sbs",add = "boxplot", add.params = list(fill = "white"))
print(p)
dev.off()

table(input_for_boxplot2$sbs,input_for_boxplot2$biological_states)

sigs_in_patients<-names(which(rowSums(table(input_for_boxplot2$sbs,input_for_boxplot2$biological_states))>=50))

pdf("emt_vs_sbs_new_atleast50patients.pdf",width=15)
input_for_boxplot_subsetpts<-input_for_boxplot2[input_for_boxplot2$sbs%in%sigs_in_patients,]
p<-ggviolin(input_for_boxplot_subsetpts, x.text.angle = 90 , "sbs", "score_emt", fill = "sbs",add = "boxplot", add.params = list(fill = "white"))
print(p)
dev.off()


library(rstatix)

list_sbs<-intersect(as.character(unique(input_for_boxplot2$sbs)),sigs_in_patients)

all_pvalue<-NULL

for(lss in list_sbs){

	if(length(unique(input_for_boxplot2[input_for_boxplot2$sbs%in%lss,"biological_states"]))>=2){
	res.aov2 <- aov(sbs_freq ~ biological_states, data = input_for_boxplot2[input_for_boxplot2$sbs%in%lss,])
	sum_test <- unlist(summary(res.aov2))
	paov<-sum_test["Pr(>F)1"]
	all_pvalue<-c(all_pvalue,paov)}else{
	all_pvalue<-c(all_pvalue,NA)
	}
}
res_anova<-data.frame(list_sbs,all_pvalue)

pdf("emt_vs_sbs_for_group_anova.pdf")

significant_sigs<-na.omit(res_anova[res_anova[,2]<=0.05,1])
input_for_boxplot3<-input_for_boxplot2[input_for_boxplot2$sbs %in% significant_sigs,]

input_for_boxplot3$biological_states<-factor(input_for_boxplot3$biological_states,levels=c("epi","mix","mes"))
#p<-ggviolin(input_for_boxplot3, x.text.angle = 90 , "sbs", "score_emt", fill = "sbs",add = "boxplot", add.params = list(fill = "white"))+coord_flip()+facet_wrap(.~biological_states)+theme(legend.position="right")
#print(p)

p<-ggviolin(input_for_boxplot3, x.text.angle = 90 , "biological_states", "sbs_freq", fill = "biological_states",add = "boxplot", add.params = list(fill = "white"),palette=c("#56b4e9","#E69f00","#EE0C0C"))+stat_summary(fun=mean, geom="line", aes(group=1))+ stat_summary(fun=mean, geom="point")+coord_flip()+facet_wrap(.~sbs)+theme(legend.position="right")
print(p)

dev.off()

pdf("emt_vs_sbs_for_group_anova_heatmap.pdf",height=2)
p1 <- ggplot(input_for_boxplot3,aes(x=samples,y=sbs,fill=sbs_freq))+
   geom_tile()+scale_fill_distiller(palette = "Greys", direction = 1)+facet_wrap(.~biological_states)+theme_bw()
 
 print(p1)
dev.off()

pdf("emt_vs_sbs_for_group_anova_hist.pdf")

input_for_boxplot3<-input_for_boxplot2[which(input_for_boxplot2$sbs %in% significant_sigs),]

# p1 <- ggplot(input_for_boxplot3,aes(x=samples,y=sbs,fill=sbs_freq))+
#   geom_tile()+scale_fill_distiller(palette = "RdYlGn", direction = -1)+facet_wrap(.~biological_states)+theme_bw()
# 
# print(p1)

p2<-gghistogram(input_for_boxplot3, x = "score_emt",
   add = "mean", rug = TRUE,
   color = "biological_states", fill = "biological_states",palette=c("#56b4e9","#EE0C0C","#E69f00"))+facet_wrap(.~sbs)+theme(legend.position="right")
print(p2)

p3<-ggplot(input_for_boxplot3, aes(x = sbs_freq)) +
  stat_ecdf(aes(color = biological_states,linetype = biological_states),
              geom = "step", size = 0.5) +
  scale_color_manual(values = c("#56b4e9","#EE0C0C","#E69f00")) + labs(y = "f(contr. mut.sig.)")+facet_wrap(.~sbs)+theme_bw()
print(p3)

p4<-ggplot(input_for_boxplot3, aes(x = score_emt)) +
  stat_ecdf(aes(color = biological_states,linetype = biological_states),
              geom = "step", size = 0.5) +
  scale_color_manual(values = c("#56b4e9","#EE0C0C","#E69f00")) + labs(y = "f(score emt)")+facet_wrap(.~sbs)+theme_bw()
print(p4)

dev.off()

# 
# Test fixed effects 
# 


# 
# Correct model by tissue 
# 

library(lme4)
library(data.table)
library(MuMIn)

list_sbs<-intersect(as.character(unique(input_for_boxplot2$sbs)),sigs_in_patients)

table(bin_sign$sbs,bin_sign$sbs_freq)

# 
# The signature must be present in more than 100 patients 
# 
bin_sign<-input_for_boxplot2
bin_sign$sbs_freq<-ifelse(bin_sign$sbs_freq >0.05,1,0)

test<-table(bin_sign$sbs,bin_sign$sbs_freq)
xsbs<-rownames(test[test[,2]>50,])

all_pvalue<-NULL
aggregated_sbs_for_lme<-data.frame()
all_rmarginal<-NULL
all_rconditional<-NULL

for(lss in intersect(list_sbs,xsbs)){
  
    input_lmer<-setDT(unique(input_for_boxplot2[input_for_boxplot2$sbs%in%lss,-3]))

    input_lmer2<-input_lmer[ , .(mean(score_emt), mean(sbs_freq),unique(tumors)), by = .(sbs, samples)]  
    colnames(input_lmer2)<-c("sbs","samples","score_emt","sbs_freq","tumors")

    print(length(unique(input_lmer2$tumors)))
    
    full.lmer <- lmer(score_emt ~ sbs_freq + (1|tumors),input_lmer2, REML = FALSE)
    rmarginal<-MuMIn::r.squaredGLMM(full.lmer)[1]
    rconditional<-MuMIn::r.squaredGLMM(full.lmer)[2]
    
    
    full.lmer <- lmer(score_emt ~ sbs_freq + (1|tumors),input_lmer2, REML = FALSE)
    
    reduced.lmer <- lmer(score_emt ~ 1 + (1|tumors),input_lmer2, REML = FALSE)
    
    res.aov2 <-anova(reduced.lmer, full.lmer)
    # sum_test <- unlist(summary(res.aov2))
    pvalue_lme<- res.aov2[2,"Pr(>Chisq)"]
    
    all_pvalue<-c(all_pvalue,pvalue_lme)
    all_rmarginal<-c(all_rmarginal,rmarginal)
    all_rconditional<-c(all_rconditional,rconditional)
    
    aggregated_sbs_for_lme<-rbind(aggregated_sbs_for_lme,input_lmer2)
    
}

res_lme4<-data.frame(intersect(list_sbs,xsbs),all_pvalue,all_rmarginal,all_rconditional)

res_lme4$padjust<-p.adjust(res_lme4[,2],"fdr")

significant_sigs<-na.omit(res_lme4[res_lme4$padjust<=0.05,1])
input_for_boxplot3<-aggregated_sbs_for_lme[aggregated_sbs_for_lme$sbs %in% significant_sigs,]

input_for_boxplot3<-unique(input_for_boxplot3)
input_for_boxplot3_sorted_emt<-input_for_boxplot3[order(input_for_boxplot3$score_emt),]
input_for_boxplot3_sorted_emt$samples<-factor(input_for_boxplot3_sorted_emt$samples,levels=unique(input_for_boxplot3_sorted_emt$samples))

library(cowplot)

pdf("emt_vs_sbs_for_group_heatmap_lmer4.pdf",height=3,width=4)

p1<-ggplot(input_for_boxplot3_sorted_emt, aes(samples, score_emt))+  theme_bw()+geom_segment(aes(x = samples, y = 0, xend = samples, yend = score_emt), data = input_for_boxplot3_sorted_emt)+theme(
                                                                                                  
                                                                                    legend.position="bottom",
                                                                                    axis.title.x=element_blank(),
                                                                                    axis.text.x=element_blank(),
                                                                                    axis.ticks.x=element_blank())+ geom_hline(yintercept=0, linetype="dashed", color = "red")

input_for_boxplot3_sorted_emt$sbs<-gsub(input_for_boxplot3_sorted_emt$sbs,pattern="SBS",replacement="")

#https://stackoverflow.com/questions/23572660/ggplot2-continuous-colors-for-discrete-scale-and-delete-a-legend

p2 <- ggplot(input_for_boxplot3_sorted_emt,aes(reorder(samples, score_emt, max),y=sbs,fill=sbs_freq))+
  geom_tile()+scale_fill_gradient2(low = "white", mid = "brown", high = "red", midpoint = 0.5)+theme_bw() + theme(legend.position="bottom",axis.title.x=element_blank(),
                                                                                        axis.text.x=element_blank(),
                                                                                        axis.ticks.x=element_blank())
p3<-plot_grid(p1,p2,nrow=2,rel_heights=c(0.3,0.7))
print(p3)

dev.off()


pdf("boxplot_results_lmer4.pdf",height=3)

boxplot_lmer<-unique(input_for_boxplot2[input_for_boxplot2$sbs%in%significant_sigs,c(1,2,3,6)])

boxplot_lmer$biological_states<-factor(boxplot_lmer$biological_states,levels=c("epi","mix","mes"))

boxplot_lmer_kappa1<-boxplot_lmer[boxplot_lmer$sbs%in%c("SBS1","SBS13"),]

my_comparisons <- list(c("epi", "mix"),c("epi","mes"),c("mix","mes"))

p<-ggboxplot(boxplot_lmer_kappa1, x.text.angle = 90 , "biological_states", "sbs_freq", fill = "biological_states", outlier.shape = NA,palette=c("#56b4e9","#E69f00","#EE0C0C"))+stat_summary(fun=median, geom="line", aes(group=1))+ stat_summary(fun=median, geom="point")+facet_wrap(.~sbs)+theme(legend.position="right")+ylim(c(0,0.65))
print(p+stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif",label.y=c(0.6,0.63,0.55)))

boxplot_lmer_kappa2<-boxplot_lmer[boxplot_lmer$sbs%in%c("SBS16"),]

boxplot_lmer_kappa2<-boxplot_lmer_kappa2[boxplot_lmer_kappa2$sbs_freq>0,]

boxplot_lmer_kappa2$sbs_freq<-log(boxplot_lmer_kappa2$sbs_freq+1,2)

p2<-ggboxplot(boxplot_lmer_kappa2, x.text.angle = 90 , "biological_states", "sbs_freq",trim=T, outlier.shape = NA, fill = "biological_states",palette=c("#56b4e9","#E69f00","#EE0C0C"))+stat_summary(fun=median, geom="line", aes(group=1))+ stat_summary(fun=median, geom="point")+facet_wrap(.~sbs)+theme(legend.position="right")+ylim(c(0,0.2))
print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif",label.y=c(0.30,0.33,0.40)))

dev.off()


# 
# uncorrected model by tissue lineare modelling
# 

bin_sign<-input_for_boxplot2
bin_sign$sbs_freq<-ifelse(bin_sign$sbs_freq >0.05,1,0)

test<-table(bin_sign$sbs,bin_sign$sbs_freq)
xsbs<-rownames(test[test[,2]>50,])

library(data.table)
# input_for_lm<-input_for_boxplot2[input_for_boxplot2$sbs_freq>0,]
input_for_lm<-input_for_boxplot2[input_for_boxplot2$sbs%in%xsbs,]

list_sbs<-intersect(as.character(unique(input_for_lm$sbs)),sigs_in_patients)

all_pvalue_lm<-NULL
aggregated_sbs_for_lm<-data.frame()

for(lss in list_sbs){
  
  print(lss)
  
  input_lmer<-setDT(unique(input_for_lm[input_for_lm$sbs%in%lss,-3]))
  
  input_lmer2<-input_lmer[ , .(mean(score_emt), mean(sbs_freq),unique(tumors)), by = .(sbs, samples)]  
  colnames(input_lmer2)<-c("sbs","samples","score_emt","sbs_freq","tumors")
  
  full.lmer <- lm(score_emt ~ sbs_freq,input_lmer2)
  pvalue_lme<-coef(summary(full.lmer))[,"Pr(>|t|)"][2]
  
  all_pvalue_lm<-c(all_pvalue_lm,pvalue_lme)
  
  aggregated_sbs_for_lm<-rbind(aggregated_sbs_for_lm,input_lmer2)
  
}

res_lm4<-data.frame(list_sbs,all_pvalue_lm)
res_lm4$padjust<-p.adjust(all_pvalue_lm,"fdr")

significant_sigs<-na.omit(res_lm4[res_lm4[,3]<=0.01,1])
input_for_boxplot3_lm<-aggregated_sbs_for_lm[aggregated_sbs_for_lm$sbs %in% significant_sigs,]

input_for_boxplot3_lm<-unique(input_for_boxplot3_lm)
input_for_boxplot3_lm_sorted_emt<-input_for_boxplot3_lm[order(input_for_boxplot3_lm$score_emt),]
input_for_boxplot3_lm_sorted_emt$samples<-factor(input_for_boxplot3_lm_sorted_emt$samples,levels=unique(input_for_boxplot3_lm_sorted_emt$samples))

library(cowplot)

pdf("emt_vs_sbs_for_group_heatmap_linear_model.pdf",height=5)

p1<-ggplot(input_for_boxplot3_lm_sorted_emt, aes(samples, score_emt))+geom_segment(aes(x = samples, y = 0, xend = samples, yend = score_emt), data = input_for_boxplot3_lm_sorted_emt)+theme_bw()+ theme(legend.position="bottom",axis.title.x=element_blank(),
                                                                                                                                                                                                         axis.text.x=element_blank(),
                                                                                                                                                                                                         axis.ticks.x=element_blank())+ geom_hline(yintercept=0, linetype="dashed", color = "red")

input_for_boxplot3_lm_sorted_emt$sbs<-gsub(input_for_boxplot3_lm_sorted_emt$sbs,pattern="SBS",replacement="")

p2 <- ggplot(input_for_boxplot3_lm_sorted_emt,aes(reorder(samples, score_emt, max),y=sbs,fill=sbs_freq))+
  geom_tile()+scale_fill_gradient2(low = "white", mid = "brown", high = "red", midpoint = 0.5)+theme_bw() + theme(legend.position="bottom",axis.title.x=element_blank(),
                                                                                        axis.text.x=element_blank(),
                                                                                        axis.ticks.x=element_blank())

p3<-plot_grid(p1,p2,nrow=2,rel_heights=c(0.2,0.8))
print(p3)

dev.off()

#
# Multiple linear regression
#
bin_sign<-melt(input_for_boxplot[,-c(1:4)])
colnames(bin_sign)<-c("sbs","sbs_freq")
bin_sign$sbs_freq<-ifelse(bin_sign$sbs_freq >0.05,1,0)

test<-table(bin_sign$sbs,bin_sign$sbs_freq)
xsbs<-rownames(test[test[,2]>50,])

input_for_boxplot_mullm<-data.frame(input_for_boxplot[,c(1:4)],input_for_boxplot[,which(colnames(input_for_boxplot)%in%xsbs)])
  
form_lm<-as.formula(paste(paste("score_emt","~",paste(grep(colnames(input_for_boxplot_mullm),pattern="SBS",value=T),collapse="+"))))
res_multiple_lm<-data.frame(coef(summary(lm(form_lm,input_for_boxplot_mullm))))
res_multiple_lm$padjust<-p.adjust(res_multiple_lm[,4],"fdr")
mut_sig_to_use<-rownames(res_multiple_lm[res_multiple_lm$padjust<=0.01,])

multiple_lm_melt<-data.frame()

for(sbs in 5:ncol(input_for_boxplot_mullm)){
  
  sbs_part<-input_for_boxplot_mullm[,c(1,2,3,4,sbs)]
  # sbs_part<-sbs_part[sbs_part[,5]>0.05,]
  res<-data.frame(sbs=colnames(input_for_boxplot_mullm)[sbs],sbs_part)
  colnames(res)[6]<-"sbs_freq"
  
  multiple_lm_melt<-rbind(multiple_lm_melt,res)
}

multiple_lm_melt2<-setDT(unique(multiple_lm_melt[multiple_lm_melt$sbs%in%mut_sig_to_use,]))

multiple_lm_melt3<-multiple_lm_melt2[ , .(mean(score_emt), mean(sbs_freq),unique(tumors)), by = .(sbs, samples)]  
colnames(multiple_lm_melt3)<-c("sbs","samples","score_emt","sbs_freq","tumors")

multiple_lm_melt3<-unique(multiple_lm_melt3)
multiple_lm_melt3_sorted_emt<-multiple_lm_melt3[order(multiple_lm_melt3$score_emt),]
multiple_lm_melt3_sorted_emt$samples<-factor(multiple_lm_melt3_sorted_emt$samples,levels=unique(multiple_lm_melt3_sorted_emt$samples))

library(cowplot)

pdf("emt_vs_sbs_for_group_heatmap_multiple_linear_model.pdf",height=5)

p1<-ggplot(multiple_lm_melt3_sorted_emt, aes(samples, score_emt))+geom_segment(aes(x = samples, y = 0, xend = samples, yend = score_emt), data = input_for_boxplot3_lm_sorted_emt)+theme_bw()+ theme(legend.position="bottom",axis.title.x=element_blank(),
                                                                                                                                                                                                         axis.text.x=element_blank(),
                                                                                                                                                                                                         axis.ticks.x=element_blank())+ geom_hline(yintercept=0, linetype="dashed", color = "red")

multiple_lm_melt3_sorted_emt$sbs<-gsub(multiple_lm_melt3_sorted_emt$sbs,pattern="SBS",replacement="")

p2 <- ggplot(multiple_lm_melt3_sorted_emt,aes(reorder(samples, score_emt, max),y=sbs,fill=sbs_freq))+
  geom_tile()+ scale_fill_gradient2(low = "white", mid = "brown", high = "red", midpoint = 0.5)+theme_bw() + theme(legend.position="bottom",axis.title.x=element_blank(),
                                                                                        axis.text.x=element_blank(),
                                                                                        axis.ticks.x=element_blank(),
                                                                                        panel.background = element_rect(fill = 'white'))

p3<-plot_grid(p1,p2,nrow=2,rel_heights=c(0.2,0.8))
print(p3)

dev.off()


#
# Create a boxplot in which discretize the SBS values
#
multiple_lm_melt3_sorted_emt_disc<-merge(multiple_lm_melt3_sorted_emt,knn_df_tcga,by.x="samples",by.y="samples")

multiple_lm_melt3_sorted_emt_disc$status_sigs<-ifelse(multiple_lm_melt3_sorted_emt_disc$sbs_freq>0.20,1,0)

pdf("boxplot_multiple_linear_model.pdf",height=5)

multiple_lm_melt3_sorted_emt_disc$biological_states<-factor(multiple_lm_melt3_sorted_emt_disc$biological_states,levels=c("epi","mix","mes"))
           
                    
for(iss in unique(multiple_lm_melt3_sorted_emt_disc$sbs)){
  
xmat<-multiple_lm_melt3_sorted_emt_disc[multiple_lm_melt3_sorted_emt_disc$sbs%in%iss,]
  
p<-ggboxplot(xmat, 
             x.text.angle = 90 , 
             "biological_states", 
             "sbs_freq", 
             outlier.shape = NA,
             fill = "biological_states", 
             palette=c("#56b4e9","#E69f00","#EE0C0C"))+stat_summary(fun=median, geom="line", aes(group=1))+stat_summary(fun=median, geom="point")+facet_wrap(~sbs, scales = "free")+theme(legend.position="right")+ylim(c(0,max(xmat$sbs_freq)-0.20))

print(p+stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif",label.y=c(max(xmat$sbs_freq)-0.20,max(xmat$sbs_freq)-0.25,max(xmat$sbs_freq)-0.30)))

}

dev.off()



# pdf("emt_vs_sbs_v2.pdf",width=15)
# 
# p1 <- input_for_boxplot2 %>% ggplot(aes(x = sbs, y = score_emt, fill = sbs)) +  
#   geom_boxplot() +
#   geom_point(position = position_jitterdodge(jitter.width = 0),
#              aes(color = factor(biological_states)), show.legend = TRUE)
# p1
# dev.off()


# #
# # Is the expression of the EPI/MES markers related with the mutational signatures?
# #
# setwd("/data/pseudospace/input_pseudospace")
# load("TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")
# 
# setwd("/data/pseudospace/HMM")
# markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
# markers_genes<-markers_genes_read[,1]
# 
# TCGA_GEXP_ALL_rid<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,]
# colnames(TCGA_GEXP_ALL_rid)<-sapply(strsplit(as.character(colnames(TCGA_GEXP_ALL_rid)),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")})
# 
# TCGA_GEXP_ALL_rid2<-t(TCGA_GEXP_ALL_rid[,-1])
# colnames(TCGA_GEXP_ALL_rid2)<-TCGA_GEXP_ALL_rid[,1]
# TCGA_GEXP_ALL_rid3<-data.frame(samples=rownames(TCGA_GEXP_ALL_rid2),log(TCGA_GEXP_ALL_rid2+1,2))
# 
# all_pvalues_genes_per_sbs<-data.frame()
# 
# setwd("/data/pseudospace/mutational_signature")
# 
# all_pvalues_genes_per_sbs<-data.frame()
# 
# for(g in colnames(TCGA_GEXP_ALL_rid3)[-1]){
# 
# 	current_genes_df<-TCGA_GEXP_ALL_rid3[,c("samples",g)]
# 	current_genes_df_sigs<-merge(current_genes_df,mutsig_with_pseudotime2)
# 	current_genes_df_sigs<-current_genes_df_sigs[order(current_genes_df_sigs$pseudospace,decreasing=T),]
# 	form_lm<-as.formula(paste(paste(g,"~",paste(grep(colnames(mutsig_with_pseudotime2),pattern="SBS",value=T),collapse="+")),"+tumors"))
# 	res_bygenes<-lm(form_lm,current_genes_df_sigs)
# 	res_df_part<-data.frame(t(summary(res_bygenes)$coefficients[,1]))[,-1]
# 	
# 	all_pvalues_genes_per_sbs<-rbind(all_pvalues_genes_per_sbs,res_df_part)
# }
# 
# all_pvalues_genes_per_sbs2<-all_pvalues_genes_per_sbs[,grep(colnames(all_pvalues_genes_per_sbs),pattern="SBS")]
# 
# rownames(all_pvalues_genes_per_sbs2)<- colnames(TCGA_GEXP_ALL_rid3)[-1]
# 
# all_pvalues_genes_per_sbs2_log<- all_pvalues_genes_per_sbs2
# 
# 
# library(pheatmap)
# 
# setwd("/data/pseudospace/mutational_signature")
# 
# pdf(paste("All_groups","genes_signatures_pvalue.pdf",sep="."),width=15)
# 
# all_pvalues_genes_per_sbs2_log[is.na(all_pvalues_genes_per_sbs2_log)]<-0
# 
# ann_row<-markers_genes_read[match(rownames(all_pvalues_genes_per_sbs2_log),markers_genes_read[,1]),]
# 
# rownames(all_pvalues_genes_per_sbs2_log)<-ann_row[,1]
# ann_row2<-data.frame(GeneType=as.factor(ann_row[,2]))
# rownames(ann_row2)<-ann_row[,1]
# 
# pheatmap(all_pvalues_genes_per_sbs2_log,annotation_row=ann_row2)
# 
# dev.off()
# 
