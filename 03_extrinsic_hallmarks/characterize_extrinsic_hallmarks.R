library(tidyverse)
library(ggpubr)
library(GSVA)
library(parallel)
library(data.table)
library(readxl)
library(VennDetail)
library(jtools)

#
# Compute hypoxia scores in my data
#

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)
load("TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")
tcga<-TCGA_GEXP_ALL
tcga[,-1]<-log(tcga[,-1]+1,2)
tcga<-setDT(tcga)

tcga2 <- data.frame(tcga[, lapply(.SD, mean), by = gene_symbol])
rownames(tcga2)<-tcga2[,1]

hypoxia_catalog<-as.data.frame(read_excel("hypoxia_genes.xlsx",col_names=F))
hypoxia_catalog_list <- split(hypoxia_catalog[,-c(1)], seq(nrow(hypoxia_catalog)))         # Convert rows to list
hypoxia_catalog_list<-lapply(hypoxia_catalog_list,FUN=function(X){as.character(na.omit(as.character(X)))})
names(hypoxia_catalog_list)<-hypoxia_catalog[,1]

list_cancers<-unique(sapply(strsplit(colnames(tcga2)[-1],split="\\."),"[[",1))
hypoxia_list_results<-vector(mode="list",length=length(list_cancers))
names(hypoxia_list_results)<-list_cancers

knn_df_tcga<-read.delim(file="HMM_results_nstates_tumors_for_states3.withEMT.txt")

input_for_hypoxia_score<-as.matrix(tcga2[,-1])
rownames(input_for_hypoxia_score)<-tcga2[,1]

# get only the expression data for epithelial tumors
input_for_hypoxia_score<-input_for_hypoxia_score[,colnames(input_for_hypoxia_score)%in%knn_df_tcga$samples]  

#create an empty list in which stores the hypoxia scores
hypoxia_list_results<-vector(mode="list",length(hypoxia_catalog_list))
  
#for each signature
  for(hypox in 1:length(hypoxia_list_results)){
   
# select the gene-expression data of the signature
  input_exp_current_hypoxia_gs<-input_for_hypoxia_score[rownames(input_for_hypoxia_score)%in%hypoxia_catalog_list[[hypox]],]

# create a temporary matrix
  scores_by_gs<-vector(mode="list",nrow(input_exp_current_hypoxia_gs))
  
# for each gene in the hypoxia signature
  for(i in 1:nrow(input_exp_current_hypoxia_gs)){
    
# select expression values for a gene
  temp_matrix<-input_exp_current_hypoxia_gs[i,]
  
# compute quantile 
  quantile_gene_hypoxia<-quantile(input_exp_current_hypoxia_gs[i,])
  
# discretize values scores https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(14)71021-6/fulltext
  temp_matrix2<-ifelse(temp_matrix>quantile_gene_hypoxia[3],1,-1)
  
  scores_by_gs[[i]]<-temp_matrix2
  
  }
  
  scores_current_gs<-apply(do.call(rbind,scores_by_gs),2,sum)
  
  hypoxia_list_results[[hypox]]<-scores_current_gs
  
  }
  
  hypoxia_results_df<-do.call(cbind,hypoxia_list_results)
  hypoxia_results_df2<-data.frame(rownames(hypoxia_results_df),hypoxia_results_df)
  colnames(hypoxia_results_df2)<-c("samples",names(hypoxia_catalog_list))
  
  hypoxia_results_df2[,1]<-as.character(hypoxia_results_df2[,1])

setwd(output_dir)
save(hypoxia_results_df2,file="PancancerHypoxia.RData")

colnames(hypoxia_results_df2)[1]<-"samples"
hypoxia_results_df2[,1]<-as.factor(hypoxia_results_df2[,1])
hypoxia_results_df2[is.na(hypoxia_results_df2)]<-0
hypoxia_results_df2<-data.frame(hypoxia_results_df2 %>% dplyr::group_by(samples) %>% dplyr::summarize(Buffa_hypoxia_genes=mean(Buffa_hypoxia_genes),
                                                                                                      Ragnum_hypoxia_genes=mean(Ragnum_hypoxia_genes),
                                                                                                      Winter_hypoxia_genes=mean(Winter_hypoxia_genes),
                                                                                                      Elvidge_hypoxia_genes=mean(Elvidge_hypoxia_genes),
                                                                                                      Eustace_hypoxia_genes=mean(Eustace_hypoxia_genes),
                                                                                                      Sorensen_hypoxia_genes=mean(Sorensen_hypoxia_genes),
                                                                                                      Hu_hypoxia_genes=mean(Hu_hypoxia_genes),
                                                                                                      Seigneuric_hypoxia_genes_early=mean(Seigneuric_hypoxia_genes_early),
                                                                                                      Seigneuric_hypoxia_genes_late=mean(Seigneuric_hypoxia_genes_late)))

#
# Integrate the genomic data
#
setwd(input_dir)
knn_df_tcga<-read.delim(file="HMM_results_nstates_tumors_for_states3.withEMT.txt")
knn_df_tcga$samples<- unlist(lapply(strsplit(as.character(knn_df_tcga$samples),split="\\."),FUN=function(X){paste(X[2:4],collapse="-")}))
colnames(knn_df_tcga)[2]<-"samples"
knn_df_tcga[,2]<-as.factor(knn_df_tcga[,2])
knn_df_tcga<-data.frame(knn_df_tcga %>% dplyr::group_by(samples) %>% dplyr::summarize(tumors=unique(tumors),
                                                                                      score_emt=mean(score_emt),
                                                                                      states=unique(states)))

# Load Aneuploidy data
setwd(input_dir)
aneuploidy_tab<-read.delim(file="aneuploidy_score_taylor_et_al2018.legacy.txt")[,c(1:13)]
aneuploidy_tab$Sample<-unlist(lapply(strsplit(as.character(aneuploidy_tab$Sample),split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
aneuploidy_tab<-aneuploidy_tab[,c(1,3,6,12,13)]
colnames(aneuploidy_tab)[1]<-"samples"
aneuploidy_tab[,1]<-as.factor(aneuploidy_tab[,1])
aneuploidy_tab[is.na(aneuploidy_tab)]<-0
aneuploidy_tab<-data.frame(aneuploidy_tab %>% dplyr::group_by(samples) %>% dplyr::summarize(AneuploidyScore.AS.=mean(AneuploidyScore.AS.),
                                                                                            Genome_doublings=mean(Genome_doublings),
                                                                                            SilentMutationspeMb=mean(SilentMutationspeMb),
                                                                                            Non.silentMutationsperMb=mean(Non.silentMutationsperMb)))



# Load Stemness data
stemness_tab<-read.delim(file="StemnessScores_RNAexp.txt")
stemness_tab$TCGAlong.id<-unlist(lapply(strsplit(as.character(stemness_tab$TCGAlong.id),split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
colnames(stemness_tab)[1]<-"samples"
stemness_tab[,1]<-as.factor(stemness_tab[,1])
stemness_tab<-data.frame(stemness_tab %>% dplyr::group_by(samples) %>% dplyr::summarize(mRNAsi=mean(mRNAsi)))

# Load CA20 data
ca20_tab<-read.delim(file="journal.pcbi.1006832.s018.txt")
ca20_tab$Sample.ID<-gsub(ca20_tab$Sample.ID,pattern="\\.",replacement="-")
ca20_tab$Sample.ID<-unlist(lapply(strsplit(as.character(ca20_tab$Sample.ID),split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
colnames(ca20_tab)[1]<-"samples"
ca20_tab<-ca20_tab[,c(1,5)]
ca20_tab[,1]<-as.factor(ca20_tab[,1])
ca20_tab<-data.frame(ca20_tab %>% dplyr::group_by(samples) %>% dplyr::summarize(mean(CA20)))
ca20_tab[,1]<-as.character(ca20_tab[,1])
colnames(ca20_tab)[2]<-"CA20"


# Load mut burden
load("PancancerMB_HMM_emtstates.gsva.RData")
MUT_BURDEN_ALL$sample<-unlist(lapply(strsplit(as.character(MUT_BURDEN_ALL$sample),split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
colnames(MUT_BURDEN_ALL)[2]<-"samples"

setwd(output_dir)
hypoxia_results_df2[,1]<-unlist(lapply(strsplit(as.character(hypoxia_results_df2$samples),split="\\."),FUN=function(X){paste(X[2:4],collapse="-")}))

resOverlapSamples<-venndetail(list(samples_pseudospace = knn_df_tcga[,1], 
                samples_aneuploidy = aneuploidy_tab[,1],
                samples_stemness = stemness_tab[,1],
                samples_ca20 = ca20_tab[,1],
                samples_hypoxia=hypoxia_results_df2[,1],
                samples_mutburden=MUT_BURDEN_ALL[,2]))

resOverlapSamples_wide<-getSet(resOverlapSamples,wide=TRUE)
commonSamples<-resOverlapSamples_wide[resOverlapSamples_wide$SharedSets==6,1]

genomics_alterations<-list(unique(knn_df_tcga[knn_df_tcga[,1]%in%commonSamples,]),
                           unique(aneuploidy_tab[aneuploidy_tab[,1]%in%commonSamples,]),
                           unique(stemness_tab[stemness_tab[,1]%in%commonSamples,]), 
                           unique(ca20_tab[ca20_tab[,1]%in%commonSamples,]),
                           unique(hypoxia_results_df2[hypoxia_results_df2[,1]%in%commonSamples,]),
                           MUT_BURDEN_ALL[MUT_BURDEN_ALL[,2]%in%commonSamples,]) %>% reduce(inner_join, by = "samples")

save(genomics_alterations,file="aneuploidy_stem_hypo_other_features_TCGA.RData")

col_to_consider<-c("states","score_emt","AneuploidyScore.AS.",
                   "Genome_doublings","SilentMutationspeMb","Non.silentMutationsperMb",
                   "mRNAsi","CA20","Buffa_hypoxia_genes",
                   "Winter_hypoxia_genes",
                   "Ragnum_hypoxia_genes",
                   "West_hypoxia_score_pan_cancer",
                   "Sorensen_hypoxia_genes",
                   "Elvidge_hypoxia_genes",
                   "Hu_hypoxia_genes",
                   "Seigneuric_hypoxia_genes_early",
                   "Seigneuric_hypoxia_genes_late",
                   "Synonymous","Non.Synonymous")

genomics_alterations_clean<-genomics_alterations[,which(colnames(genomics_alterations)%in%col_to_consider)]

genomics_alterations_clean$states<-gsub(genomics_alterations_clean$states,pattern="1",replacement="pEMT")
genomics_alterations_clean$states<-gsub(genomics_alterations_clean$states,pattern="2",replacement="epi")
genomics_alterations_clean$states<-gsub(genomics_alterations_clean$states,pattern="3",replacement="mes")

genomics_alterations_clean$states<- factor(genomics_alterations_clean$states, levels = c("epi", "pEMT", "mes"))


setwd(output_dir)

#
# Boxplot Aneuploidy
#
my_comparisons<-list(c("epi","mes"),c("pEMT","mes"),c("epi","pEMT"))

pdf("boxplot_aneuploidy.pdf")

ggviolin(genomics_alterations_clean, 
         x = "states", 
         y = "AneuploidyScore.AS.", 
         fill = "states",
         palette  = c("blue2","orange2", "red2"),
         add = "boxplot", 
         add.params = list(fill = "white"))+
         stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")


dev.off()
#
# Boxplot CA20
#

pdf("boxplot_CA20.pdf")

ggviolin(genomics_alterations_clean, 
         x = "states", 
         y = "CA20", 
         fill = "states",
         palette  = c("blue2","orange2", "red2"),
         add = "boxplot", 
         add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")
dev.off()

#
# Hypoxia
#
list_hypoxia<-grep(colnames(genomics_alterations_clean),pattern="hypoxia",value=T)

pdf("boxplot_hypoxia_score_pan_cancer.pdf")

for(i in list_hypoxia){
  
p<-ggviolin(genomics_alterations_clean, 
         x = "states", 
         y = i, 
         fill = "states",
         palette  = c("blue2","orange2", "red2"),
         add = "boxplot", 
         add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+ ggtitle(i)

print(p)

}

dev.off()


#
# Boxplot MB_Synonymous
#

pdf("boxplot_MB_Synonymous.pdf")

tempmatrix<-genomics_alterations_clean

tempmatrix$Synonymous<-log(tempmatrix$Synonymous+1,2)

ggviolin(tempmatrix, 
         x = "states", 
         y = "Synonymous", 
         fill = "states",
         palette  = c("blue2","orange2", "red2"),
         add = "boxplot", 
         add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")


dev.off()
#
# Boxplot MB_Synonymous
#

pdf("boxplot_MB_Non.Synonymous.pdf")

tempmatrix<-genomics_alterations_clean

tempmatrix$Non.Synonymous<-log(tempmatrix$Non.Synonymous+1,2)

ggviolin(tempmatrix, 
         x = "states", 
         y = "Non.Synonymous", 
         fill = "states",
         palette  = c("blue2","orange2", "red2"),
         add = "boxplot", 
         add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")
dev.off()


#
# Hypoxia vs emt score
#

pdf("scatterplot_emt_vs_hypoxia.pdf")

for(i in list_hypoxia){
  
  p<-ggscatter(genomics_alterations_clean, 
              x = "score_emt", 
              y = i, 
              color = "states",
              add = "reg.line",
              conf.int = TRUE,
              add.params = list(color = "black", fill = "lightgray"),
              palette  = c("blue2","orange2", "red2"))+geom_hline(yintercept=0, linetype="dashed")+geom_vline(xintercept=0, linetype="dashed")+stat_cor(method = "pearson", label.x = 2, label.y = 20)
  border()                                         
  
  xplot <- ggdensity(genomics_alterations_clean, "score_emt", fill = "states",
                     palette  = c("blue2","orange2", "red2"))
  
  yplot <- ggdensity(genomics_alterations_clean, i, fill = "states", 
                     palette  = c("blue2","orange2", "red2"))+rotate()
    
  sp <- p + rremove("legend")
  
  yplot <- yplot + clean_theme() + rremove("legend")
  xplot <- xplot + clean_theme() + rremove("legend")
  # Arranging the plot using cowplot
  library(cowplot)
  p2<-plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", rel_widths = c(2, 1), rel_heights = c(1, 2))
  print(p2)
}

dev.off()


#
# Correlation between the variables
#

input_for_cor<-genomics_alterations_clean

library(corrplot)
library(RColorBrewer)

pdf("cor_between_genomicmarkers.pdf",width=15,height=15)
M <-cor(input_for_cor[,-2],use="complete.obs")
corrplot(M, method="number")
corrplot(M, method="color")
dev.off()

#
# LM to measure EMT
#

input_lm<-genomics_alterations_clean
colnames(input_lm)[3]<-"AneuploidyScore"

lm.fit2=lm(score_emt ~ (AneuploidyScore+Genome_doublings+mRNAsi+CA20+Buffa_hypoxia_genes+Synonymous+Non.Synonymous)^2,data=input_lm)

library(dotwhisker)
library(interactions)

pdf("lm_emt_vs_genomics_features_dotwisher.pdf")
dwplot(lm.fit2)+ xlab("Coefficient Estimate") + ylab("")+geom_vline(xintercept=0,lty=2)
plot_summs(lm.fit2)
dev.off()

model1 <- lm(score_emt~Buffa_hypoxia_genes*Non.Synonymous, input_lm)

pdf("interaction_chart_Y_EMT_Buffa_hypoxia_vs_non.synonymous.pdf")
sim_slopes(lm.fit2, pred = Buffa_hypoxia_genes, modx = Non.Synonymous, jnplot = TRUE)
interact_plot(lm.fit2, pred = Buffa_hypoxia_genes, modx = Non.Synonymous)
interact_plot(lm.fit2, pred = Buffa_hypoxia_genes, modx = Non.Synonymous, plot.points = TRUE)
dev.off()

#
# GLM binary epi vs mes
#
col_to_consider<-c("states","score_emt","AneuploidyScore",
                   "Genome_doublings","SilentMutationspeMb","Non.silentMutationsperMb",
                   "mRNAsi","CA20","Buffa_hypoxia_genes",
                   "Synonymous","Non.Synonymous")

colnames(genomics_alterations_clean)[3]<-"AneuploidyScore"

input_lm_epi_mes<-genomics_alterations_clean[input_lm$states%in%c("epi","mes"),colnames(genomics_alterations_clean)%in%col_to_consider]
mylogit_epi_mes <- glm(states~score_emt+AneuploidyScore+SilentMutationspeMb+Non.silentMutationsperMb+mRNAsi+CA20+Buffa_hypoxia_genes+Synonymous+Non.Synonymous, data = input_lm_epi_mes, family = "binomial",maxit=200)

input_lm_epi_pEMT<-genomics_alterations_clean[input_lm$states%in%c("epi","pEMT"),colnames(genomics_alterations_clean)%in%col_to_consider]
mylogit_epi_pEMT <- glm(states~score_emt+AneuploidyScore+SilentMutationspeMb+mRNAsi+CA20+Buffa_hypoxia_genes+Synonymous, data = input_lm_epi_pEMT, family = "binomial",maxit=100)

pdf("test.pdf")
plot_summs(mylogit_epi_mes, mylogit_epi_pEMT)
dev.off()

#
# What is the relation between the TME and hypoxia?
#
setwd(input_dir)
load("PancancerTME.RData")
all_signatures<-unique(unlist(lapply(list_results,rownames)))
all_samples<-sum(unlist(lapply(list_results,ncol)))

#
# Create the matrix with the values of the TME
#
library(plyr)

input_tme<-data.frame(matrix(0,nrow=1,ncol=length(all_signatures)))
colnames(input_tme)<-all_signatures

for(i in 1:length(list_results)){
  
  sub_tme<-data.frame(t(list_results[[i]]))
  input_tme<-rbind.fill(input_tme,sub_tme)
  
}

input_tme<-data.frame(samplesID=as.character(unlist(lapply(list_results,colnames))),input_tme[-1,])
input_tme[,1]<-unlist(lapply(strsplit(as.character(input_tme$samplesID),split="\\."),FUN=function(X){paste(X[2:4],collapse="-")}))

col_to_consider<-c("samples","score_emt","Buffa_hypoxia_genes",
                   "Winter_hypoxia_genes",
                   "Ragnum_hypoxia_genes",
                   "West_hypoxia_score_pan_cancer",
                   "Sorensen_hypoxia_genes",
                   "Elvidge_hypoxia_genes",
                   "Hu_hypoxia_genes",
                   "Seigneuric_hypoxia_genes_early",
                   "Seigneuric_hypoxia_genes_late")
                  

input_hypoxia<-genomics_alterations[,which(colnames(genomics_alterations)%in%col_to_consider)]

input_cor_tme_hypoxia<-merge(input_tme,input_hypoxia,by.x="samplesID",by.y="samples")

library(corrplot)
library(RColorBrewer)

setwd(output_dir)

pdf("cor_TME_hypoxia.pdf",width=15,height=15)
M <-cor(input_cor_tme_hypoxia[,-1],use="complete.obs")
corrplot(M, method="number")
corrplot(M, method="color")
dev.off()

# There is a sort of relation between hypoxia and TME?

relation_hypoxia_tme<-input_cor_tme_hypoxia
relation_hypoxia_tme$Buffa_hypoxia_genes_thr<-factor(ifelse(relation_hypoxia_tme$Buffa_hypoxia_genes>0,"High_hypoxia","Low_Hypoxia"), levels = c("Low_Hypoxia", "High_hypoxia"))
relation_hypoxia_tme$score_emt_thr<-factor(ifelse(relation_hypoxia_tme$score_emt>0,"High_EMT","Low_EMT"), levels = c("Low_EMT", "High_EMT"))
relation_hypoxia_tme$combined_emt_hypoxia<-paste(relation_hypoxia_tme$score_emt_thr,relation_hypoxia_tme$Buffa_hypoxia_genes_thr,sep="_vs_")
  
columns_to_consider<-c(                     
  "B_cells","Cytotoxic_cells",                
  "Dendritic_cells","Eosinophils",                    
  "Macrophages","Mast_cells",                     
  "NK_cells", "Neutrophils",                    
  "T_cells_CD4" ,"T_cells_CD8",                    
  "T_cells_gamma_delta","T_regulatory_cells",             
  "Macrophages_M1","Macrophages_M2",                 
  "Endothelial","Fibroblasts",                    
  "Monocytes","Plasma_cells",                   
  "Immune_Score","combined_emt_hypoxia")          

relation_hypoxia_tme2<-melt(relation_hypoxia_tme[,colnames(relation_hypoxia_tme)%in%columns_to_consider])

pdf("boxplot_TME_vs_Hypoxia_with_EMT_High_Low.pdf")
ggviolin(relation_hypoxia_tme2, x = "combined_emt_hypoxia", y = "value", fill = "combined_emt_hypoxia",
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("High_hypoxia","Low_Hypoxia")),method="wilcox" ,label = "p.signif")+facet_wrap(~variable)+
  stat_summary(fun.y=median, geom="line", aes(group=1))
dev.off()

#
#  Do a model to explore the relation between multiple genomics features, TME and fibroblasts
# 

setwd(input_dir)

load("PancancerFibroblastTME.ssgsea.RData")

all_signatures<-unique(unlist(lapply(fibroblast_list_results,rownames)))
all_samples<-sum(unlist(lapply(fibroblast_list_results,ncol)))

#
# Create the matrix with the values of the TME
#

input_fibroblast<-data.frame(matrix(0,nrow=1,ncol=length(all_signatures)))
colnames(input_fibroblast)<-all_signatures

for(i in 1:length(fibroblast_list_results)){
  
  sub_fibroblast<-data.frame(t(fibroblast_list_results[[i]]))
  input_fibroblast<-rbind.fill(input_fibroblast,sub_fibroblast)
  
}

input_fibroblast<-data.frame(samplesID=as.character(unlist(lapply(fibroblast_list_results,colnames))),input_fibroblast[-1,])
input_fibroblast[,1]<-unlist(lapply(strsplit(as.character(input_fibroblast$samplesID),split="\\."),FUN=function(X){paste(X[2:4],collapse="-")}))
input_fibroblast<-input_fibroblast[,-c(3,45:47)]
remove_fibroblasts<-c("tumors","samples","samples2","LUNG","CRC","Breast_pCAF",
                      "Breast_NMF","Breast_sCAF","General_FS1",
                      "General_FS2","Qian_FS2","C7_RGS5","Nurmik_Fibroblast_markers",
                      "Wu_CAF_markers","Wu_dPVL_cells_TNBC","Wu_iCAFs_TNBC","Wu_imPVL_cells_TNBC","Han_CAF_markers",
                      "PROSTATE_CAF_markers")
col_fibroblast<-setdiff(colnames(input_fibroblast),remove_fibroblasts)

tme_fibroblast<-merge(input_fibroblast,input_tme,by="samplesID",all.x=F,all.y=F)

col_to_consider<-c("samples","states","score_emt","AneuploidyScore.AS.",
                   "Genome_doublings","SilentMutationspeMb","Non.silentMutationsperMb",
                   "mRNAsi","CA20","Buffa_hypoxia_genes",
                   "Synonymous","Non.Synonymous")

genomics_alterations_clean<-unique(genomics_alterations[,which(colnames(genomics_alterations)%in%col_to_consider)])

samples_genomics<-genomics_alterations_clean$samples
samples_tme<-tme_fibroblast$samplesID

common_samples<-intersect(samples_genomics,samples_tme)

mat_genomics<-unique(genomics_alterations_clean[genomics_alterations_clean[,1]%in%common_samples,])
colnames(mat_genomics)[1]<-"samplesID"

mat_tme<-unique(tme_fibroblast[tme_fibroblast[,1]%in%common_samples,])
col_tme<-colnames(mat_tme)

input_gf_tme_fb<-unique(merge(mat_genomics,
                       mat_tme,
                       by="samplesID",
                       all.x=F,
                       all.y=F))

setwd(output_dir)

save(input_gf_tme_fb,file="aneuploidy_stem_hypo_TME_Fibroblasts_TCGA.RData")

col_to_consider<-setdiff(c("states","score_emt","AneuploidyScore.AS.",
                   "Genome_doublings","SilentMutationspeMb","Non.silentMutationsperMb",
                   "mRNAsi","CA20","Buffa_hypoxia_genes",
                   "Winter_hypoxia_genes",
                   "Ragnum_hypoxia_genes",
                   "West_hypoxia_score_pan_cancer",
                   "Sorensen_hypoxia_genes",
                   "Elvidge_hypoxia_genes",
                   "Hu_hypoxia_genes",
                   "Seigneuric_hypoxia_genes_early",
                   "Seigneuric_hypoxia_genes_late",
                   "Synonymous","Non.Synonymous",
                   col_fibroblast,col_tme),remove_fibroblasts)

#
# Run GLM model
#

input_gf_tme_fb2<-input_gf_tme_fb[,colnames(input_gf_tme_fb)%in%col_to_consider]
input_gf_tme_fb2$states<-gsub(input_gf_tme_fb2$states,pattern="1",replacement="pEMT")
input_gf_tme_fb2$states<-gsub(input_gf_tme_fb2$states,pattern="2",replacement="epi")
input_gf_tme_fb2$states<-gsub(input_gf_tme_fb2$states,pattern="3",replacement="mes")

input_gf_tme_fb2$states<- factor(input_gf_tme_fb2$states, levels = c("epi", "pEMT", "mes"))

pdf("cor_TME_Fibroblast_GenomicFeatures.pdf",width=15,height=15)
M <-cor(input_gf_tme_fb2[,-c(1,3)],use="complete.obs")
corrplot(M, method="color",col=colorRampPalette(c("blue","white","red"))(200),type="upper")
M2<-M
M2[M2<0.2]<-0
corrplot(M2, method="color",col=colorRampPalette(c("blue","white","red"))(200),type="upper")
dev.off()

#
# GLM between MES vs EPI 
#

input_glm_epi_mes_extended<-input_gf_tme_fb2[input_gf_tme_fb2$states%in%c("epi","mes"),]
input_glm_epi_mes_extended$states<-as.character(input_glm_epi_mes_extended$states)
input_glm_epi_mes_extended$states<-gsub(input_glm_epi_mes_extended$states,patter="epi",replacement="0")
input_glm_epi_mes_extended$states<-gsub(input_glm_epi_mes_extended$states,patter="mes",replacement="1")
input_glm_epi_mes_extended$states<-as.factor(input_glm_epi_mes_extended$states)
input_glm_epi_mes_extended$states<-relevel(input_glm_epi_mes_extended$states, ref = "0")

formula_no_interactions<-paste("states~",paste(colnames(input_glm_epi_mes_extended[,-which(colnames(input_glm_epi_mes_extended)%in%c("samplesID","score_emt","states"))]),collapse="+"),sep="")

require(sjPlot)

#Run a GLM model to select the most relevant features
mylogit_epi_mes <- glm(formula_no_interactions, data = input_glm_epi_mes_extended, family = "binomial")
features_epi_mes<-summary(mylogit_epi_mes)$coeff[-1,4]
significant_terms_epi_mes<-names(features_epi_mes[which(features_epi_mes<0.01)])

#create new formulas
formula_no_interactions_sig<-paste("states~",paste(significant_terms_epi_mes,collapse="+"),sep="")
formula_interactions_sig<-paste("states~","(",paste(paste(significant_terms_epi_mes,collapse="+"),")^2"),sep="")

mylogit_epi_mes_reduced <- glm(formula_no_interactions_sig, data = input_glm_epi_mes_extended, family = "binomial")

pdf("mes_vs_epi_TME_Fibroblast_GenomicFeatures_reduced_model_nointeractions_pvalue0.05.pdf")
plot_model(mylogit_epi_mes_reduced, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,value.offset = .3)
dev.off()

mylogit_epi_mes_interactions <- glm(formula_interactions_sig, data = input_glm_epi_mes_extended, family = "binomial")
coef_int_mes_epi<-exp(cbind(mes_vs_epi=coef(mylogit_epi_mes_interactions), confint.default(mylogit_epi_mes_interactions)))
stats_int_mes_epi<-cbind(coef_int_mes_epi, p_value=summary(mylogit_epi_mes_interactions)$coeff[,4])
stats_int_mes_epi<-stats_int_mes_epi[-1,]
stats_int_mes_epi<-stats_int_mes_epi[stats_int_mes_epi[,"p_value"]<=0.05,]
write.table(stats_int_mes_epi,file="mes_vs_epi_TME_Fibroblast_GenomicFeatures_reduced_model_interactions_pvalue0.05.txt",sep="\t",row.names=T,quote=F)

#
# GLM between pEMT vs EPI 
#

input_glm_epi_pEMT_extended<-input_gf_tme_fb2[input_gf_tme_fb2$states%in%c("epi","pEMT"),]
input_glm_epi_pEMT_extended$states<-as.character(input_glm_epi_pEMT_extended$states)
input_glm_epi_pEMT_extended$states<-gsub(input_glm_epi_pEMT_extended$states,patter="epi",replacement="0")
input_glm_epi_pEMT_extended$states<-gsub(input_glm_epi_pEMT_extended$states,patter="pEMT",replacement="1")
input_glm_epi_pEMT_extended$states<-as.factor(input_glm_epi_pEMT_extended$states)
input_glm_epi_pEMT_extended$states<-relevel(input_glm_epi_pEMT_extended$states, ref = "0")

formula_no_interactions<-paste("states~",paste(colnames(input_glm_epi_pEMT_extended[,-which(colnames(input_glm_epi_pEMT_extended)%in%c("samplesID","score_emt","states"))]),collapse="+"),sep="")

require(sjPlot)

#Run a GLM model to select the most relevant features
mylogit_epi_pEMT <- glm(formula_no_interactions, data = input_glm_epi_pEMT_extended, family = "binomial")
features_epi_pEMT<-summary(mylogit_epi_pEMT)$coeff[-1,4]
significant_terms_epi_pEMT<-names(features_epi_pEMT[which(features_epi_pEMT<0.05)])

#create new formulas
formula_no_interactions_sig<-paste("states~",paste(significant_terms_epi_pEMT,collapse="+"),sep="")
formula_interactions_sig<-paste("states~","(",paste(paste(significant_terms_epi_pEMT,collapse="+"),")^2"),sep="")

mylogit_epi_pEMT_reduced <- glm(formula_no_interactions_sig, data = input_glm_epi_pEMT_extended, family = "binomial")

pdf("pEMT_vs_epi_TME_Fibroblast_GenomicFeatures_reduced_model_nointeractions_pvalue0.05.pdf")
plot_model(mylogit_epi_pEMT_reduced, vline.color = "red",transform = NULL,show.p = T,show.values = TRUE,value.offset = .3)
dev.off()

mylogit_epi_pEMT_interactions <- glm(formula_interactions_sig, data = input_glm_epi_pEMT_extended, family = "binomial")
coef_int_pEMT_epi<-exp(cbind(pEMT_vs_epi=coef(mylogit_epi_pEMT_interactions), confint.default(mylogit_epi_pEMT_interactions)))
stats_int_pEMT_epi<-cbind(coef_int_pEMT_epi, p_value=summary(mylogit_epi_pEMT_interactions)$coeff[,4])
stats_int_pEMT_epi<-stats_int_pEMT_epi[-1,]
stats_int_pEMT_epi<-stats_int_pEMT_epi[stats_int_pEMT_epi[,"p_value"]<=0.05,]
terms_int_pEMT_epi<-rownames(stats_int_pEMT_epi)
write.table(stats_int_pEMT_epi,file="pEMT_vs_epi_TME_Fibroblast_GenomicFeatures_reduced_model_interactions_pvalue0.05.txt",sep="\t",row.names=T,quote=F)

#
# scatter plot hypoxia and fibroblasts for which there is an interaction
#

pdf("mes_vs_epi_hypoxia_vs_fibroblasts_glm_interactions.pdf",width=5,height=5)

list_fibroblasts=c("Esophageal_CAF_markers_carcinogenesis","BLADDER_CAF_markers")

for(i in list_fibroblasts){
  
part_input<-input_glm_epi_mes_extended[,colnames(input_glm_epi_mes_extended)%in%c("states","Buffa_hypoxia_genes",i)]

p<-ggscatter(part_input, 
             x = "Buffa_hypoxia_genes", 
             y = i, 
             color = "states",
             add = "reg.line",
             conf.int = TRUE,
             add.params = list(color = "black", fill = "lightgray"),
             palette  = c("blue2","red2"))+geom_hline(yintercept=0, linetype="dashed")+geom_vline(xintercept=0, linetype="dashed")+stat_cor(method = "pearson")
border()                                         

xplot <- ggdensity(part_input, "Buffa_hypoxia_genes", fill = "states",palette  = c("blue2","red2"))

yplot <- ggdensity(part_input, i, fill = "states", palette  = c("blue2","red2"))+rotate()

sp <- p + rremove("legend")

yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
p2<-cowplot::plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", rel_widths = c(2, 0.5), rel_heights = c(0.5, 2))
print(p2)

}

dev.off()

#
# scatter plot hypoxia and fibroblasts for which there is an interaction
#

pdf("pEMT_vs_epi_hypoxia_vs_fibroblasts_glm_interactions.pdf",width=5,height=5)

list_fibroblasts=c("THYM_CAF_markers","X3_HCC_Early_CAF_markers","Plasma_cells","Esophageal_CAF_markers_proliferation")

for(i in list_fibroblasts){
  
  part_input<-input_glm_epi_pEMT_extended[,colnames(input_glm_epi_pEMT_extended)%in%c("states","Buffa_hypoxia_genes",i)]
  
  p<-ggscatter(part_input, 
               x = "Buffa_hypoxia_genes", 
               y = i, 
               color = "states",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "black", fill = "lightgray"),
               palette  = c("blue2","orange2"))+geom_hline(yintercept=0, linetype="dashed")+geom_vline(xintercept=0, linetype="dashed")+stat_cor(method = "pearson")
  border()                                         
  
  xplot <- ggdensity(part_input, "Buffa_hypoxia_genes", fill = "states",palette  = c("blue2","orange2"))
  
  yplot <- ggdensity(part_input, i, fill = "states", palette  = c("blue2","orange2"))+rotate()
  
  sp <- p + rremove("legend")
  
  yplot <- yplot + clean_theme() + rremove("legend")
  xplot <- xplot + clean_theme() + rremove("legend")
  # Arranging the plot using cowplot
  library(cowplot)
  p2<-cowplot::plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", rel_widths = c(2, 0.5), rel_heights = c(0.5, 2))
  print(p2)
  
}

dev.off()
