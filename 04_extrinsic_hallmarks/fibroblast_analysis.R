library(plyr)
library(circlize)
library(ComplexHeatmap)
library(ggpubr)
library(readxl)
library(data.table)
library(GSVA)
library(parallel)

setwd("/home/guidantoniomt/pseudospace/input_pseudospace")
load("TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")
tcga<-TCGA_GEXP_ALL
tcga[,-1]<-log(tcga[,-1]+1,2)
tcga<-setDT(tcga)

tcga2 <- data.frame(tcga[, lapply(.SD, mean), by = gene_symbol])
rownames(tcga2)<-tcga2[,1]


setwd("/home/guidantoniomt/pseudospace/pathway_characterization/fibroblast_analysis")
fibro_catalog<-as.data.frame(read_excel("Fibroblast_catalog.xlsx",col_names=F))
fibro_catalog_list <- split(fibro_catalog[,-c(1,2)], seq(nrow(fibro_catalog)))         # Convert rows to list
fibro_catalog_list<-lapply(fibro_catalog_list,FUN=function(X){as.character(na.omit(as.character(X)))})
names(fibro_catalog_list)<-fibro_catalog[,2]

purity_data<-read.delim(file="/home/guidantoniomt/datasets/TCGA_supp/purity/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",stringsAsFactors=F)
purity_data$samples<- unlist(lapply(strsplit(as.character(purity_data$sample),split="-"),FUN=function(X){paste(X[1:4],collapse="-")}))
purity_data<-purity_data[-which(is.na(purity_data$purity)),]

knn_df_tcga<-read.delim(file="/home/guidantoniomt/pseudospace/HMM/HMM_results_nstates_tumors_for_states3.withEMT.txt",stringsAsFactors=F)
knn_df_tcga$samples2<-knn_df_tcga$samples
knn_df_tcga$samples<- unlist(lapply(strsplit(as.character(knn_df_tcga$samples),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")}))

emt_with_purity<-merge(knn_df_tcga,purity_data,by="samples")

setwd("/home/guidantoniomt/pseudospace/pathway_characterization")

emt_with_purity$states<-gsub(emt_with_purity$states,pattern=1,replacement="pEMT")
emt_with_purity$states<-gsub(emt_with_purity$states,pattern=2,replacement="epi")
emt_with_purity$states<-gsub(emt_with_purity$states,pattern=3,replacement="mes")

emt_with_purity$states<-as.factor(emt_with_purity$states)
emt_with_purity$states<- factor(emt_with_purity$states, levels = c("epi", "pEMT", "mes"))

list_cancers<-unique(sapply(strsplit(colnames(tcga2)[-1],split="\\."),"[[",1))
fibroblast_list_results<-vector(mode="list",length=length(list_cancers))
names(fibroblast_list_results)<-list_cancers

for(i in 1:length(list_cancers)){
  
  tum=list_cancers[i]
  
  idx_tum<-grep(unlist(sapply(strsplit(colnames(tcga2),split="\\."),"[[",1)),pattern=tum)
  
  matrix_tum<-tcga2[,idx_tum]
  matrix_tum$var<-apply(matrix_tum,1,var)
  matrix_tum2<-matrix_tum[order(matrix_tum$var,decreasing=T),-c(ncol(matrix_tum))][1:round(nrow(matrix_tum)*30/100),]
  
  res<-gsva(as.matrix(matrix_tum2),fibro_catalog_list,method="ssgsea",parallel.sz=60)
  
  fibroblast_list_results[[i]]<-res
}

setwd("/home/guidantoniomt/pseudospace/pathway_characterization/fibroblast_analysis")
save(fibroblast_list_results,file="PancancerFibroblastTME.ssgsea.RData")

#
# Analysis Fibroblast
#

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
input_fibroblast[,1]<-unlist(lapply(strsplit(as.character(input_fibroblast$samplesID),split="\\."),FUN=function(X){paste(X[2:5],collapse="-")}))

matrix_pseudotime_fibroblast<-merge(emt_with_purity,input_fibroblast,by.x="samples",by.y="samplesID")
matrix_pseudotime_fibroblast<-matrix_pseudotime_fibroblast[,-c(17,59:61)]

quantile_purity<-quantile(matrix_pseudotime_fibroblast$purity)
matrix_pseudotime_fibroblast$purity_thresholded<-rep("0",nrow(matrix_pseudotime_fibroblast))
matrix_pseudotime_fibroblast$purity_thresholded[quantile_purity<=quantile_purity[2]]<-"low_purity"
matrix_pseudotime_fibroblast$purity_thresholded[quantile_purity>=quantile_purity[4]]<-"high_purity"
matrix_pseudotime_fibroblast$purity_thresholded<-gsub(matrix_pseudotime_fibroblast$purity_thresholded,pattern="0",replacement="intermediate_purity")

#
# Heatmap with the fibroblasts values
#
annotation_colors=list(tumors=c("ACC"="#996b2f","BLCA"="#724fb3","BRCA"="#ba449c",
                                "CESC"="#adba3d","CHOL"="#8f2e43","COAD"="#191970",
                                "ESCA"="#5060ac","HNSC"="#65b06f","KICH"="#B266FF",
                                "KIRC"="#969841","KIRP"="#bf6d29","LIHC"="#d1972c",
                                "LUAD"="#36dee6","LUSC"="#cc4b3e","OV"="#d97eb8",
                                "PAAD"="#c43d7e","PRAD"="#d077d8","READ"="#FFB266","SKCM"="#452c7c",
                                "STAD"="#43c29e","THYM"="#de6991","THCA"="#757ee8",
                                "UCS"="#88321a","UCEC"="#75b550","UVM"="#8f2931"),
                       states=c("epi"="blue2","pEMT"="orange2","mes"="red2"),
                       purity_thresholded=c("low_purity"="green","intermediate_purity"="orange","high_purity"="red"))

col_fun = colorRamp2(c(min(matrix_pseudotime_fibroblast$score_emt), 0, max(matrix_pseudotime_fibroblast$score_emt)), c("deepskyblue3", "gainsboro", "darkorange2"))

col_ann<-HeatmapAnnotation(tumors=matrix_pseudotime_fibroblast$tumors,
                           score_emt=matrix_pseudotime_fibroblast$score_emt,
                           states=matrix_pseudotime_fibroblast$states,
                           purity_thresholded=matrix_pseudotime_fibroblast$purity_thresholded,
                           col=list(
                             tumors=annotation_colors[[1]],
                             score_emt=col_fun,
                             states=annotation_colors[[2]],
                             purity_thresholded=annotation_colors[[3]]))

pdf(paste("Fibroblasts_heatmap_with_emt_nstates_forced_purity",3,"ssgsea","pdf",sep="."),width=12)
matrix_pseudotime_fibroblast[is.na(matrix_pseudotime_fibroblast)]<-0

mat_fibroblast<-t(matrix_pseudotime_fibroblast[,-c(1:15,ncol(matrix_pseudotime_fibroblast))])

p<-Heatmap(mat_fibroblast,
           show_column_names=F,
           top_annotation = col_ann,
           column_split=as.factor(matrix_pseudotime_fibroblast$purity_thresholded),
           row_km=3)
print(p)
dev.off()

pdf(paste("Fibroblasts_heatmap_with_emt_nstates_forced_EMTstates",3,"ssgsea","pdf",sep="."),width=12)
matrix_pseudotime_fibroblast[is.na(matrix_pseudotime_fibroblast)]<-0

mat_fibroblast<-t(matrix_pseudotime_fibroblast[,-c(1:15,ncol(matrix_pseudotime_fibroblast))])

p<-Heatmap(mat_fibroblast,
           show_column_names=F,
           top_annotation = col_ann,
           column_split=as.factor(matrix_pseudotime_fibroblast$states),
           row_km=3)
print(p)
dev.off()


#
# Heatmap of only C1-C11 genes
#
gene_sets_ctype<-c("Pancreas","CRC","General_FS2","C1_KCNN3","C6_CALB2","C7_MYH11","C9_CFD","C10_COMP","C11_SERPINE1","Nurmik_CAF_markers","Gastric_CAF_markers","Esophageal_CAF_markers_premalignant","Esophageal_CAF_markers_carcinogenesis","Esophageal_CAF_markers_invasion","BCC_CAF_markers","HNSCC_CAF_markers","HCC_CAF_markers","X12_HCC_CAF_markers","8_HCC_Advanced_CAF_markers")

pdf(paste("Selected_Fibroblasts_heatmap_with_emt_nstates_forced_purity",3,"ssgsea","pdf",sep="."),width=12)
matrix_pseudotime_fibroblast[is.na(matrix_pseudotime_fibroblast)]<-0

mat_fibroblast<-t(matrix_pseudotime_fibroblast[,-c(1:15,ncol(matrix_pseudotime_fibroblast))])
mat_fibroblast<-mat_fibroblast[rownames(mat_fibroblast)%in%gene_sets_ctype,]
matrix_pseudotime_fibroblast$states<-factor(matrix_pseudotime_fibroblast$states,levels=c("2","1","3"))

p<-Heatmap(mat_fibroblast,
           show_column_names=F,
           top_annotation = col_ann,
           column_split=as.factor(matrix_pseudotime_fibroblast$states),
           row_km=3)
print(p)
dev.off()


COL11A1_FS_markers_fibroblasts<-c("COL11A1","COL1A1","COL1A2","COL3A1","FAP","COL5A1",
  "CTHRC1","SULF1","VCAN","FN1","SFRP2","OLFML2B","COL6A3",
  "INHBA","ASPN","ADAMTS12","LUM","NTM","SPARC","AEBP1","CDH11",
  "GREM1","ITGBL1","FNDC1","PRRX1","MMP11","FBN1","CTSK","COL8A1",
  "MMP2","WISP1","COMP","SFRP4","HTRA3","EPYC","PCOLCE","PODNL1","ANTXR1",
  "ZNF469","CORIN","THY1","MFAP5","GLT8D2","ISLR","COL6A2","PDGFRB","LOXL2","
  OMD","P4HA3","LOX","DCN","TIMP2","WNT2","NOX4","TNFSF4","EMILIN1","NID2",
  "COL6A1","COL8A2","TNFAIP6","RCN3","FAM26E","COL5A3","DACT1","GFPT2",
  "MXRA8","CCDC80","ANGPTL2","MMP14","KIF26B","FIBIN","MRC2","KCND2","ALPK2",
  "ADAMTS14","TWIST1","ACTA2","LAMA4","TAGLN","CERCAM","FBLN2","CALD1","ECM2",
  "MEDAG","TMEM119","TGFB3","EDNRA","FSTL1","C1QTNF6","SCARF2","RGS4","FAM180A"
  ,"OLFML1","MSRB3","CRISPLD2","GPC6","SOX11","GPR1","GAS1","SERPINF1","PODN","CD248",
  "PLXDC1","MMP13","COL16A1","CPZ","VGLL3","COL15A1","SGIP1","CMTM3","LOXL1","NALCN",
  "COL24A1","PLAU","MFAP2","LTBP2","MATN3","THBS1","CTGF","DPYSL3","ARSI","ADAMTS16","GRP",
  "RAB31","HTRA1","CLMP","IGFL2","BICC1","PDGFRL","FAM101A","SGCD","CSMD2","LAMP5","SERPINH1","PDGFRA",
  "BNC2","COL4A1", "IBSP","ITGA5","MMP16","ST6GALNAC5","ADAM19","MEIS3","APCDD1L",
  "XIRP1","ADAMTS4","PDPN","COLEC12","MAGEL2","C1S","DDR2","BMP1","TENM3","DACT3",
  "DIO2","SPOCK1","GGT5","PALLD","NTNG2","LGALS1","CPXM1","ADAMTS6","HMCN1","NNMT","PMP22","TGFBI",
  "CLEC11A","CCL11","ZFPM2","GUCA1A","SPOCD1","TGFB1I1","KERA","NUAK1","ALDH1L2",
  "RUNX2","CLEC5A","CCIN","MSC","PRR16","LZTS1","FILIP1L","ZNF521","PMEPA1","TMEM158","CHN1","TSHZ3",
  "MICAL2","LMCD1","SNAI2","GALNT15","CSGALNACT2","BCAT1","AXLP")

exp_markers_fibroblasts<-tcga2[tcga2[,1]%in%COL11A1_FS_markers_fibroblasts,]
rownames(exp_markers_fibroblasts)<-exp_markers_fibroblasts[,1]
exp_markers_fibroblasts2<-exp_markers_fibroblasts[,colnames(exp_markers_fibroblasts)%in%knn_df_tcga$samples2]
exp_markers_fibroblasts2<-exp_markers_fibroblasts2[,match(knn_df_tcga$samples2,colnames(exp_markers_fibroblasts2))]
exp_markers_fibroblasts2<-apply(exp_markers_fibroblasts2,1,FUN=function(X){(X-mean(X))/sd(X)})

annotation_colors=list(tumors=c("ACC"="#996b2f","BLCA"="#724fb3","BRCA"="#ba449c",
                                "CESC"="#adba3d","CHOL"="#8f2e43","COAD"="#191970",
                                "ESCA"="#5060ac","HNSC"="#65b06f","KICH"="#B266FF",
                                "KIRC"="#969841","KIRP"="#bf6d29","LIHC"="#d1972c",
                                "LUAD"="#36dee6","LUSC"="#cc4b3e","OV"="#d97eb8",
                                "PAAD"="#c43d7e","PRAD"="#d077d8","READ"="#FFB266","SKCM"="#452c7c",
                                "STAD"="#43c29e","THYM"="#de6991","THCA"="#757ee8",
                                "UCS"="#88321a","UCEC"="#75b550","UVM"="#8f2931"),
                       states=c("epi"="steelblue2","pEMT"="orange2","met"="red2"))

library(FSelector)
input_fibro_exp2<-cbind(knn_df_tcga$states,exp_markers_fibroblasts2)
colnames(input_fibro_exp2)[1]<-"states"

res_gain<-information.gain(states~., data.frame(input_fibro_exp2))
thr<-quantile(res_gain[,1])[4]
idx_genes<-which(res_gain>thr)

genes_to_save_FS<-rownames(res_gain)[idx_genes]

setwd("/home/guidantoniomt/pseudospace/pathway_characterization/fibroblast_analysis")

mat_markers_fibroblast<-data.frame(sampleID=rownames(exp_markers_fibroblasts2),exp_markers_fibroblasts2)
mat_markers_fibroblast2<-merge(knn_df_tcga,mat_markers_fibroblast,by.x="samples2",by.y="sampleID")

mat_markers_fibroblast2$states<-gsub(mat_markers_fibroblast2$states,pattern="1",replacement="pEMT")
mat_markers_fibroblast2$states<-gsub(mat_markers_fibroblast2$states,pattern="2",replacement="epi")
mat_markers_fibroblast2$states<-gsub(mat_markers_fibroblast2$states,pattern="3",replacement="met")

mat_markers_fibroblast2$states<-as.factor(mat_markers_fibroblast2$states)

test<-melt(mat_markers_fibroblast2[,c(5:ncol(mat_markers_fibroblast2))])

pdf("test.pdf")
ggviolin(test, x = "states", y = "value", fill = "states",
         palette = c("blue2", "orange2", "red2"),
         add = "boxplot", add.params = list(fill = "white"))
dev.off()

col_fun = colorRamp2(c(min(knn_df_tcga$score_emt), 0, max(knn_df_tcga$score_emt)), c("deepskyblue3", "gainsboro", "darkorange2"))

ann_data<-mat_markers_fibroblast2[,c(1:5)]

col_ann<-HeatmapAnnotation(tumors=ann_data$tumors,
                           score_emt=ann_data$score_emt,
                           states=ann_data$states,
                           col=list(
                             tumors=annotation_colors[[1]],
                             score_emt=col_fun,
                             states=annotation_colors[[2]]))

pdf(paste("Expression_markers_fibroblast_emt_states",3,"ssgsea","pdf",sep="."),width=12)
p<-Heatmap(t(mat_markers_fibroblast2[,colnames(mat_markers_fibroblast2)%in%genes_to_save_FS]),
           show_column_names=F,
           top_annotation = col_ann,
           column_split = factor(mat_markers_fibroblast2$states,levels=(c("epi","pEMT","met"))),
	   row_km=3,
	   cluster_column_slices = FALSE)
print(p)
dev.off()

#
# Boxplot for each fibroblast gene set and EMT state
#

mat_fibroblast_emt<-merge(knn_df_tcga,input_fibroblast[,-c(3,45:47)],by.x="samples",by.y="samplesID")

mat_fibroblast_emt$states<-gsub(mat_fibroblast_emt$states,pattern=1,replacement="pEMT")
mat_fibroblast_emt$states<-gsub(mat_fibroblast_emt$states,pattern=2,replacement="epi")
mat_fibroblast_emt$states<-gsub(mat_fibroblast_emt$states,pattern=3,replacement="mes")

mat_fibroblast_emt$states<-as.factor(mat_fibroblast_emt$states)
mat_fibroblast_emt$states<- factor(mat_fibroblast_emt$states, levels = c("epi", "pEMT", "mes"))

input_for_boxplot_fibroblast<-melt(mat_fibroblast_emt[,c(4,6:ncol(mat_fibroblast_emt))])

#
# Create an heatmap with the activation of the fibroblasts
#
library(rstatix)

res_aov_TME<-input_for_boxplot_fibroblast %>%
  group_by(variable) %>%
  anova_test(value ~ states)

#they are all significants
res_aov_TME_df<-data.frame(res_aov_TME)
res_aov_TME_df$padjust<-p.adjust(res_aov_TME_df[,6],"BH")

aggregated_tme<-aggregate(value~states+variable,input_for_boxplot_fibroblast,mean)

setwd("/home/guidantoniomt/pseudospace/pathway_characterization/fibroblast_analysis")

pdf("heatmap_TME_HMM_3_bytumor.ssgsea.Fibroblasts.pdf")
aggregated_tme[which(aggregated_tme[,3]<0),3]<-0
p <- ggplot(aggregated_tme,aes(x=variable,y=states,fill=value))+
  geom_tile(colour="white",size=0.2)+scale_fill_distiller(palette = "Blues", direction = 1)+theme_bw()
print(p+coord_polar())
dev.off()





#
# Boxplot TME vs STATES
#

my_comparisons<-list(c("mes","epi"),c("epi","pEMT"),c("mes","pEMT"))

pdf("boxplot_Fibroblasts_HMM_3_bytumor.ssgsea.pdf",width=15,height=15)

ggviolin(input_for_boxplot_fibroblast, x = "states", y = "value", fill = "states",
         palette = c("blue2", "orange2", "red2"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+facet_wrap(~variable)+
  stat_summary(fun.y=median, geom="line", aes(group=1))
dev.off()

#
# Multinomial logistic regression 
#

#remove some markers from the multinomial analysis
remove_fibroblasts<-c("tumors","samples","samples2","LUNG","CRC","Breast_pCAF",
"Breast_NMF","Breast_sCAF","General_FS1",
"General_FS2","Qian_FS2","C7_RGS5","Nurmik_Fibroblast_markers",
"Wu_CAF_markers","Wu_dPVL_cells_TNBC","Wu_iCAFs_TNBC","Wu_imPVL_cells_TNBC","Han_CAF_markers",
"PROSTATE_CAF_markers")

library(nnet)

setwd("/home/guidantoniomt/pseudospace/pathway_characterization/fibroblast_analysis")

### Select train and test sets, dividing into 2/3 for train, 1/3 for test:
set.seed(19875)

col_to_use<-setdiff(c("states","score_emt",colnames(mat_fibroblast_emt)),remove_fibroblasts)

input_mnlogit<-mat_fibroblast_emt[,colnames(mat_fibroblast_emt)%in%col_to_use]

randTrainingSet <- floor(0.70 * nrow(input_mnlogit))
train_ind <- sample(seq_len(nrow(input_mnlogit)), size = randTrainingSet)

dat.train <- input_mnlogit[train_ind,]
dat.test <- input_mnlogit[-train_ind,]

# Fit multinomial logistic regression model:
glm.fit=multinom(states~0+., data=dat.train[,-c(1)])
summary(glm.fit)
z <- summary(glm.fit)$coefficients/summary(glm.fit)$standard.error
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

pmultinom<-t(data.frame(p))
RRmultinom<-t(summary(glm.fit)$coefficients)
colnames(RRmultinom)<-paste("RR",colnames(RRmultinom),sep="_")
out_multinom<-cbind(pmultinom,RRmultinom)

write.table(out_multinom,file="Fibroblasts_multinomial_res_hmm_3states.ssgsea.txt",sep="\t",row.names=T,quote=F)

## Predictions:
predicted=predict(glm.fit,dat.test[,-c(1:2)],type="probs")
bpp=cbind(dat.test[,-2], predicted)

bpp2 <- melt(bpp,id.vars = setdiff(colnames(bpp),c("epi","pEMT","mes")))

colnames(bpp2)[32:33] <- c("Prediction","Probability")

all_tme<-setdiff(colnames(bpp2)[-1],c("states","Prediction","Probability"))

library("gridExtra")

input_for_heatmap_Fibroblasts<-data.frame()

pdf("Fibroblast_multinomial_with_emt_ssgsea.pdf",width=10)

for(tmemn in 1:length(all_tme)){
  
  print(tmemn)
  current_tme<-all_tme[tmemn]
  print(current_tme)
  
  bpp3<-bpp2[,colnames(bpp2)%in%c(current_tme,"score_emt","Probability","Prediction")]
  colnames(bpp3)[2]<-"tme"
  
  input_for_heatmap_Fibroblasts<-rbind(input_for_heatmap_Fibroblasts,data.frame(immcell=current_tme,bpp3))

  ptme <- ggplot(bpp3, aes(x = tme, y = Probability, colour = Prediction)) +
    geom_line() + facet_grid(Prediction ~ ., scales="free")+
    ylim(c(0,1))+
    ylab("Probability")+labs(title=all_tme[tmemn]) + geom_smooth(method = "lm",color="black")
  
#  bpp3$tme_status<-rep(NULL,nrow(bpp3))
  
#  bpp3$tme_status<-ifelse(bpp3$tme>0,"High","Low")
#  bpp3$tme_status<-factor(bpp3$tme_status, levels = c("Low", "High"))
  
#  my_comparisons <- list(c("High", "Low"))
  
#  pbox<-ggviolin(bpp3, x = "tme_status", y = "score_emt", fill = "tme_status",
#                 palette = c("blue2","red2"),
#                 add = "boxplot", add.params = list(fill = "white"))+stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")
  
#  bpp3$tme_status_prob<-rep(NULL,nrow(bpp3))
#  bpp3[which(is.na(bpp3$Probability)),"Probability"]<-0
#  bpp3$tme_status_prob<-ifelse(bpp3$Probability>=0.8,"High_Prob","Low_Prob")
#  bpp3$tme_status_prob<-factor(bpp3$tme_status_prob, levels = c("High_Prob", "Low_Prob"))
  
#  bpp3$combined_states<- paste(bpp3$tme_status_prob,bpp3$tme_status,sep="-")
#  bpp3$combined_states<-factor(bpp3$combined_states, levels = c("Low_Prob-Low", "Low_Prob-High","High_Prob-Low","High_Prob-High"))
  
#  if(length(which(is.na(bpp3$tme)))!=0){
#    bpp3<-bpp3[-which(is.na(bpp3$tme)),]
#  }else{
#    bpp3<-bpp3
#  }
  
#  my_comparisons<- list(c("Low_Prob-Low", "High_Prob-High"))
  
#  pbox2<-ggviolin(bpp3, x = "combined_states", y = "score_emt", fill = "combined_states",
#                  add = "boxplot", add.params = list(fill = "white"))+stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+facet_wrap(~Prediction,nrow=3,ncol=1)+stat_summary(fun.y=median, geom="line", aes(group=1))
  
#  my_comparisons <- list(c("epi", "mes"),c("epi", "pEMT"),c("pEMT", "mes"))
  
#  pbox3<-ggviolin(bpp3, x = "Prediction", y = "score_emt", fill = "Prediction",
#                  add = "boxplot", add.params = list(fill = "white"))+stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+facet_wrap(~combined_states,nrow=2)+stat_summary(fun.y=median, geom="line", aes(group=1))
  
#  bpp4<-bpp3[bpp3$combined_states %in% c("High_Prob-High","High_Prob-Low"),]
#  bpp4$combined_states<-as.character(bpp4$combined_states)
#  bpp4$combined_states<-as.factor(bpp4$combined_states)
  
#  bpp4<-bpp4[,c(1,2,3,7)]
  
#  my_comparisons<- list(c("High_Prob-Low", "High_Prob-High"))
  
#  library(ggpubr)
#  library(rstatix)
  
#  idxzero<-which(table(as.character(bpp4$Prediction),bpp4$combined_states)==0)
#  idxcol<-ncol(table(as.character(bpp4$Prediction),bpp4$combined_states))
  
#  if(length(idxzero)>=1 | idxcol == 1){
    
    # melt_check<-melt(table(bpp4$Prediction,bpp4$combined_states))
    # emt_group<-melt_check[idxzero,1]
    # group_prob<-as.character(melt_check[idxzero,2])
    
    # bpp5<-bpp4[-which(bpp4$combined_states%in%group_prob & bpp4$Prediction %in% emt_group ),]
    
#    status_stat<-"NO"
#    bpp5<-bpp4
    
#  } else {
    
#    status_stat<-"YES"
#    bpp5<-bpp4
    
#  }
  
#  if(status_stat!="NO"){
    
#    stat.test <- bpp5 %>% group_by(Prediction) %>% wilcox_test(score_emt ~ combined_states)
    
#    stat.test <- stat.test %>%
#      add_xy_position(x = "Prediction", dodge = 0.8)
    
#    stat.test2 <- bpp4 %>%
#      group_by(combined_states) %>%
#      wilcox_test(score_emt ~ Prediction)
    
#    dfstats<-data.frame(stat.test2)[,c(1,3,4,8)]
    
#   stat.test2 <- stat.test2 %>%
#      add_xy_position(x = "Prediction", dodge = 0.8)
    
#    pbox4<-ggviolin(bpp4, x = "Prediction", y = "score_emt", color = "combined_states",palette = c("#00AFBB", "#E7B800"),add = "boxplot", add.params = list(fill = "white"))+stat_pvalue_manual(stat.test,  label = "p", tip.length = 0)+geom_hline(yintercept=0, linetype="dashed")
    
#  }else{
    
#    pbox4<-ggviolin(bpp4, x = "Prediction", y = "score_emt", color = "combined_states",palette = c("#00AFBB", "#E7B800"),add = "boxplot", add.params = list(fill = "white"))+geom_hline(yintercept=0, linetype="dashed")
    
#  }
  
  #i am interested only in the samples that are High Prob and High and low TME
#  ss <- tableGrob(dfstats)

print(ptme)  
#  ptme_box<-grid.arrange(ptme,pbox2,ncol=2,nrow=1,widths=c(8,3))
  
#  print(ptme_box)
  
#  print(pbox3)
  
#  ptme_box2<-grid.arrange(pbox4,ss,ncol=1,nrow=2)
  
#  print(ptme_box2)
  
}

dev.off()

pemt_tme<-rownames(out_multinom[out_multinom[,"pEMT"]<=0.01 & out_multinom[,"RR_pEMT"]>log(100),])
mes_tme<-rownames(out_multinom[out_multinom[,"mes"]<=0.01 & out_multinom[,"RR_mes"]>log(100),])

terms_to_consider<-unique(c(mes_tme,pemt_tme))

library(VennDetail)

ven <- venndetail(list(pemt = pemt_tme, mes = mes_tme))

write.table(file="Fibroblasts_multinomial_res_hmm_3states.ssgsea.TERMSFORCHARTS.txt",result(ven, wide = TRUE),row.names=T,quote=F,sep="\t")

pdf("heatmap_multnom_fibroblasts.pdf",height=4)
p <- ggplot(input_for_heatmap_Fibroblasts[input_for_heatmap_Fibroblasts$immcell%in%terms_to_consider,],aes(x=tme,y=immcell,fill=Probability))+
  geom_tile()+ geom_vline(xintercept = 0, linetype="dotted",color = "black", size=0)+scale_fill_distiller(palette = "RdYlGn", direction = -1)+facet_grid(.~Prediction,scales="free_x")+theme_bw()
print(p)
dev.off()

pdf("heatmap_multnom_fibroblasts_scatterplotgrid.pdf",height=4)
p2<-ggplot(input_for_heatmap_Fibroblasts[input_for_heatmap_Fibroblasts$immcell%in%terms_to_consider,], aes(x=tme, y=Probability, color=Prediction)) +geom_point(shape=18,size=0.5)+ylim(c(0,1)) +geom_smooth(method = "lm", formula = y ~ poly(x, 10),se=TRUE, linetype="dashed",color="blue") + facet_grid(immcell~Prediction)
print(p2)
dev.off()

library(ggpubr)

input_for_boxplot_tme<-melt(input_mnlogit[,2:ncol(input_mnlogit)])
my_comparisons<-list(c("mes","epi"),c("epi","pEMT"),c("mes","pEMT"))

pdf("boxplot_Fibroblasts_HMM_3_bytumor.ssgsea.selectedTME_multinom.pdf")

p1<-ggviolin(input_for_boxplot_tme[input_for_boxplot_tme$variable %in% terms_to_consider,], x = "states", y = "value", fill = "states",
         palette = c("#56b4e9", "#E69f00", "#EE0C0C"),
         add = "boxplot", add.params = list(fill = "white"),ylim=c(-0.5,0.5))+
         stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+facet_wrap(~variable,nrow=5,ncol=5)+
         stat_summary(fun.y=median, geom="line", aes(group=1))

print(p1)

dev.off()

