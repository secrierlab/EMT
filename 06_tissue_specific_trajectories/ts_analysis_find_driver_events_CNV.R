library("readxl")
library("data.table")
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

setwd("/data/pseudospace/")

getFeaturesSS<-function(res_lasso,ss=800){
  
  df_coef2<-do.call(rbind,res_lasso)
  
  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]
  
  features<-data.frame(table(df_coef2[,1]))
  
  features_freq_df<-features[order(features[,2],decreasing=T),]
  
  features_to_select<-features_freq_df[features_freq_df[,2]>=ss,1]
  
  return(features_to_select)
}


getCoefSS<-function(res_lasso,features_to_select){
  
  df_coef2<-do.call(rbind,res_lasso)
  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]
  
  features<-data.frame(table(df_coef2[,1]))
  
  input_for_boxplot<-df_coef2[which(df_coef2[,1]%in%features_to_select),]
  colnames(input_for_boxplot)[2]<-"values"
  
  return(input_for_boxplot)
  
}

FisherByModule<-function(x,table_nsamples_for_state){
  
  #              M1  REST
  # LASP1 (YES)  A - B
  # LASP1 (NO)   C - D
  #
  all_pvalue<-NULL
  
  for(m in 1:length(x)){
    
    # A) number of samples in a specific module (m) with a copy number alteration 
    nsamples_current_module<-x[m]
    
    # C) number of samples in "m" (e.g.M1) without alterations in LASP1. How many m1 samples are without the alteration?
    nsamples_module_withoutgene<-table_nsamples_for_state[m]-x[m]
    
    # B) the number of samples in which an alteration is not in the current module (e.g. LASP1, M1) (Rest group)
    nsamples_restmodule_withgene<-sum(x)-x[m]
    
    # D) number of samples without an alteration in any module
    
    # 1) Total number of samples  (except that one in the module of interest) and samples with alterations in M1
    # The final number of samples is given by the sum obtained in 1) minus all the samples with an alterations in M1 (excluded the current module m)
    # check sum(mtx_fisher) and sum(table_nsamples_for_state)
    
    
    # nsamples<-table_nsamples_for_state[m]
    nsamples<-sum(table_nsamples_for_state)-table_nsamples_for_state[m]
    
    allsamples_nogene<-nsamples-(sum(x)-x[m])
    
    vector_for_fisher<-c(nsamples_current_module,nsamples_restmodule_withgene,nsamples_module_withoutgene,allsamples_nogene)
    
    mtx_fisher<-matrix(vector_for_fisher,ncol=2,byrow=2)
    
    if(mtx_fisher[1,1]!=0){
      
      pvalue<-fisher.test(mtx_fisher,alternative="greater")$p.value
      
    } else {
      
      pvalue<-1
      
    }
    
    all_pvalue<-c(all_pvalue,pvalue)
    
  }
  
  names(all_pvalue)<-names(table_nsamples_for_state)
  
  all_pvalue<-p.adjust(all_pvalue)
  
  return(all_pvalue)
  
}

setwd("/data/pseudospace/")

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
cnv_markers_mes_vs_epi<- sapply(strsplit(grep(sub(grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement=""),pattern="focal",value=T),split="_"),"[[",1)

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
cnv_markers_pemt_vs_epi<-sapply(strsplit(grep(sub(grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement=""),pattern="focal",value=T),split="_"),"[[",1)

load("HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_mix_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
cnv_markers_mes_vs_mix<-sapply(strsplit(grep(sub(grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement=""),pattern="focal",value=T),split="_"),"[[",1)


cosmic_table<-read.csv(file="/data/pseudospace/Census_allMon Aug 17 13_43_26 2020.csv",stringsAsFactors=F)
cosmic_genes<-as.character(cosmic_table[,1])

#
# Read the data for the analysis
#

load("/data/pseudospace/A549_mapped_seurat_correct_all_timecourse.RData")
LUAD_scores_EMT<-df_scores_EMT

load("/data/pseudospace/MCF7_mapped_seurat_correct_all_timecourse.RData")
BRCA_scores_EMT<-df_scores_EMT

#
# Read input files for the analysis
#

setwd("/data/pseudospace/find_drivers_events")

LUAD_clusters <- read_excel("LUAD_groups_for_Lucie.xlsx")
BRCA_clusters <- read_excel("data_for_lucie_BRCA_PRAD_OV.xlsx","BRCA")

LUAD_input<-merge(LUAD_clusters,LUAD_scores_EMT[,c(1:2)],by.x="Samples",by.y="Samples")
BRCA_input<-merge(BRCA_clusters,BRCA_scores_EMT[,c(1:2)],by.x="Samples",by.y="Samples")

#
# Merge with hypoxia scores
#
load("/home/guidantoniomt/pseudospace/pathway_characterization/broad_genomic_features/PancancerHypoxia.RData")
hypoxia_results_df2[,1]<-unlist(lapply(strsplit(gsub(hypoxia_results_df2[,1],pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[2:5],collapse="-")}))

#
# Create a list containing the experiments 
#

list_input<-vector(mode="list",2)
list_input[[1]]<-LUAD_input
list_input[[2]]<-BRCA_input
names(list_input)<-c("LUAD","BRCA")

list_order<-vector(mode="list",2)
list_order[[1]]<-c("cluster3","cluster4","cluster5","cluster1","cluster2")
list_order[[2]]<-c("cluster4","cluster2","cluster1","cluster5","cluster3")

list_pvalue<-c(0.1,0.1)
thr_minimum<-c(0.1,0.1)

#
# Iteration for cancer type
#

for(scexp in 1:length(list_input)){
  
  current_results<-list_input[[scexp]]
  current_results[,1]<-unlist(lapply(strsplit(gsub(current_results[,1],pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[2:5],collapse="-")}))
  
  current_cancer<-names(list_input)[scexp]
  
  vector_order_to_use<-list_order[[scexp]]
  
  current_pvalue<-list_pvalue[scexp]
  
  thr_min<-thr_minimum[scexp]
  
  dir_cancer<-paste("/TCGA/harmonized/",current_cancer,sep="")
  
  setwd(dir_cancer)
  
  cnv_file<-fread(file=grep(dir(),pattern="_gene_level_copy.harmonized.txt",value=T),data.table=F)
  cnv_file[,2]<-sapply(strsplit(as.character(cnv_file[,2]),split="\\."),"[[",1)
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=cnv_file[,2],mart= mart)

  cnv_file2<-merge(G_list,cnv_file,by.x="ensembl_gene_id",by.y="Gene Symbol")[,-c(3:5)]
  cnv_file3<-cnv_file2[which(cnv_file2$hgnc_symbol %in% cosmic_genes),]
  colnames(cnv_file3)<-unlist(lapply(strsplit(colnames(cnv_file3),split="-"),FUN=function(X){paste(X[1:4],collapse="-")}))
  colnames(cnv_file3)[1:2]<-c("ENSG","gene_symbol")
  
  cnv_file4<-t(cnv_file3[,-c(1:2)])
  colnames(cnv_file4)<-cnv_file3[,2]
  
  cnv_file5<-data.frame(Samples=colnames(cnv_file3)[-c(1:2)],cnv_file4)
  
  input_cnv_to_use <- merge(current_results,cnv_file5,by.x="Samples",by.y="Samples",all.x=F,all.y=F)
  input_cnv_to_use$clusters<-as.factor(input_cnv_to_use$clusters)
  
  input_cnv_to_use2<-melt(input_cnv_to_use[,-c(3,4)])
  
  # 
  #  Get the data with copy number alterations (+1,/-1)
  # 
  
  input_cnv_to_use3<-input_cnv_to_use2[input_cnv_to_use2$value!=0,]
  list_dep_amp<-split(input_cnv_to_use3,input_cnv_to_use3$value)
  names(list_dep_amp)<-c("depletion","amplification")
  
  RES_AMP_DEP<-vector(mode="list",2)           # store the genes to consider
  RES_AMP_DEP_MATRICES<-vector(mode="list",2) #  store the number of patients with CNV alterations
  RES_AMP_DEP_SIG<-vector(mode="list",2) #  store the number of patients with CNV alterations
  
  for(events in 1:length(list_dep_amp)){
    
  current_df_event<-list_dep_amp[[events]]
  
  clusters_for_genes_df<-table(current_df_event$clusters,current_df_event$variable)

  # remove the common CNV alterations across the clusters
  clusters_for_genes_df_temp<-clusters_for_genes_df
  clusters_for_genes_df_temp[clusters_for_genes_df_temp>1]<-1

  idx_genes_in_common<-which(apply(clusters_for_genes_df_temp,2,sum) %in% nrow(clusters_for_genes_df_temp))

  clusters_for_genes_df<-clusters_for_genes_df[,-idx_genes_in_common]

  genes_min_samples<-names(which(apply(clusters_for_genes_df,2,sum)>=10))

  clusters_for_genes_df2<-clusters_for_genes_df[,which(colnames(clusters_for_genes_df)%in%genes_min_samples)]

  table_nsamples_for_state<-table(current_results[,2])
  
  fisher_res_table<-apply(clusters_for_genes_df2,2,FisherByModule,table_nsamples_for_state)
  
  #remove the genes that are not significant in any cluster
  check_significance<-function(x){ifelse(length(x[x<0.05])==0,"discard","save")}
  check_genes_module<-apply(fisher_res_table,2,check_significance)
  sig_genes2<-names(check_genes_module[check_genes_module%in%"save"])
  
  fisher_res_table2<-fisher_res_table[,sig_genes2]
  
  freq_for_genes<-t(apply(clusters_for_genes_df2[,colnames(clusters_for_genes_df2)%in%sig_genes2],2,FUN=function(X){X/sum(X)}))
  
  selected_freq_cnv_for_clusters<-apply(freq_for_genes,2,FUN=function(X){names(X[X>thr_min])})

  names(selected_freq_cnv_for_clusters)<-paste("cluster",colnames(selected_freq_cnv_for_clusters),sep="")
  
  RES_AMP_DEP[[events]]<-selected_freq_cnv_for_clusters
  RES_AMP_DEP_MATRICES[[events]]<-freq_for_genes
  RES_AMP_DEP_SIG[[events]]<-fisher_res_table2
    
  }
  
  names(RES_AMP_DEP)<-c("depletion","amplification")
  
  
  # 
  #  Get the significant genes
  # 
  
  # genes_for_lme<-unique(input_for_heatmap2$genes)
  
  SIGNIFICANT_CNV_GENES<-vector(mode="list",length(RES_AMP_DEP))
  
  # 
  #   Start iterations on the genes that are amplified or depleted
  # 
  
  for(cnv_events in 1:length(RES_AMP_DEP)){
  
  genes_for_lme<-unique(unlist(RES_AMP_DEP[[cnv_events]]))
  genes_for_lme<-genes_for_lme[!is.na(genes_for_lme)]
  
  current_cnv_df<-merge(current_results,input_cnv_to_use2,by.x=c("Samples","clusters"),by.y=c("Samples","clusters"))
  colnames(current_cnv_df)[5]<-"Hugo_Symbol"
  
  if(names(RES_AMP_DEP)[cnv_events]=="depletion"){
    
  current_cnv_df2<-current_cnv_df[current_cnv_df$value<= -1,]
    
  }else{
    
  current_cnv_df2<-current_cnv_df[current_cnv_df$value>=1,]
  
  }
  
 
  # use directly the genes found during the pre-processing
  genes_significant<-genes_for_lme

  SIGNIFICANT_CNV_GENES[[cnv_events]]<-genes_significant  

  emt_clusters<-current_cnv_df2[current_cnv_df2$Hugo_Symbol %in% genes_significant,]
  emt_clusters2<-aggregate(score_emt~clusters+Hugo_Symbol, emt_clusters, median)
  emt_clusters2$clusters<-paste("cluster",emt_clusters2$clusters,sep="")
  
  matrices_freq<-RES_AMP_DEP_MATRICES[[cnv_events]]
  colnames(matrices_freq)<-paste("cluster",colnames(matrices_freq),sep="")
  
  matrices_freq2<-melt(matrices_freq)
  colnames(matrices_freq2)<-c("genes","clusters","freq_perc")
  
  input_circular_heatmap<-merge(matrices_freq2,emt_clusters2,by.x=c("clusters","genes"),by.y=c("clusters","Hugo_Symbol"))
  
  fisher_res<-melt(RES_AMP_DEP_SIG[[cnv_events]])
  colnames(fisher_res)<-c("clusters","genes","sig")
  
  fisher_res$clusters<-paste("cluster",fisher_res$clusters,sep="")
  
  input_circular_heatmap2<-merge(input_circular_heatmap,fisher_res,by.x=c("clusters","genes"),by.y=c("clusters","genes"))
  
  # Create circular plot for copy number only
  input_circular_heatmap2$clusters<-factor(input_circular_heatmap2$clusters,levels=vector_order_to_use)
  
  input_circular_heatmap2$siglog<- -log(input_circular_heatmap2$sig,10)
  
  # Remove the modules not statistically significant
  # 0 is the -log(1,10) and -log(0.05,10)
  if(length(which(input_circular_heatmap2$siglog==0))!=0){
    
  input_circular_heatmap3<-input_circular_heatmap2[-which(input_circular_heatmap2$siglog==0),]
  
  } else {
  
  input_circular_heatmap3<-input_circular_heatmap2
    
  }
  
  if(length(which(input_circular_heatmap3$siglog < -log(0.05,10)))!=0){
    
  input_circular_heatmap3<-input_circular_heatmap3[-which(input_circular_heatmap3$siglog < -log(0.05,10)),]
  
  }else{
    
  input_circular_heatmap3<-input_circular_heatmap3
    
  }
  
  setwd("/data/pseudospace/")
  
  status<-names(RES_AMP_DEP)[cnv_events]
  
  output_pdf<-paste("CNV_circular_heatmap_emt.",current_cancer,".",status,"July.pdf",sep="") 
  
  pdf(output_pdf)
  
  if(names(RES_AMP_DEP)[cnv_events]=="amplification"){
    
  p<-ggplot(input_circular_heatmap3, aes(y = factor(clusters),
                       x = factor(genes))) +        ## global aes
    geom_tile(aes(fill = freq_perc))+
    scale_fill_gradientn(colours = (brewer.pal(11,"Reds"))) +         ## to get the rect filled
    geom_point(aes(colour = score_emt, 
                   size = siglog))  +    ## geom_point for circle illusion
    scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")))+
    theme_bw()+
    coord_polar(theta="x")
  
  print(p)
  
  }
  
  if(names(RES_AMP_DEP)[cnv_events]=="depletion"){
    
   p<-ggplot(input_circular_heatmap3, aes(y = factor(clusters),
                                           x = factor(genes))) +        ## global aes
      geom_tile(aes(fill = freq_perc))+
      scale_fill_gradientn(colours = (brewer.pal(11,"Greens"))) +         ## to get the rect filled
      geom_point(aes(colour = score_emt, 
                     size = siglog))  +    ## geom_point for circle illusion
      scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")))+
      theme_bw()+
      coord_polar(theta="x")
    
    print(p)
    
    
  }
  
  dev.off()
  
  }
  
  
  names(SIGNIFICANT_CNV_GENES)<-c("depletion","amplification")
  
  library(VennDetail)
  
  ven <- venndetail(list(depletion=SIGNIFICANT_CNV_GENES[["depletion"]],
                         amplification=SIGNIFICANT_CNV_GENES[["amplification"]],
                         mes_vs_epi = cnv_markers_mes_vs_epi,
                         pemt_vs_epi = cnv_markers_pemt_vs_epi,
                         mex_vs_mix = cnv_markers_mes_vs_mix))
  
  setwd("/data/pseudospace/")
  
  output_pdf_venn<-paste("CNV_venn_for_clusters",current_cancer,".July.pdf",sep="")
  output_txt_venn<-paste("CNV_venn_for_clusters",current_cancer,".July.txt",sep="")
  
  pdf(output_pdf_venn)
  
  p<-plot(ven, type = "upset")
  
  print(p)
  
  results2<-result(ven, wide = TRUE)
  
  write.table(results2,file=output_txt_venn,sep="\t",row.names=F,quote=F)
  
  dev.off()
  
}
