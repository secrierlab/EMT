library("readxl")
library("data.table")
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

setwd("/home/guidantoniomt/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

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

setwd("/home/guidantoniomt/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
mut_genes_mes_vs_epi<-sapply(strsplit(grep(markers_mes_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
mut_genes_pemt_vs_epi<-sapply(strsplit(grep(markers_pemt_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)

load("HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_mix_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
mut_genes_mes_vs_mix<-sapply(strsplit(grep(markers_mes_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)

#
# Read the data of cosmic
#

cosmic_table<-read.csv(file="/home/guidantoniomt/pseudospace/ml_for_ppt/Census_allMon Aug 17 13_43_26 2020.csv",stringsAsFactors=F)
cosmic_genes<-as.character(cosmic_table[,1])

#
# Read the data for the analysis
#

load("/home/guidantoniomt/pseudospace/res_multiple_pseudospace/A549_mapped_seurat_correct_all_timecourse.RData")
LUAD_scores_EMT<-df_scores_EMT

load("/home/guidantoniomt/pseudospace/res_multiple_pseudospace/MCF7_mapped_seurat_correct_all_timecourse.RData")
BRCA_scores_EMT<-df_scores_EMT

#
# Read input files for the analysis
#

setwd("/home/guidantoniomt/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")

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

list_pvalue<-c(0.1,0.05)
thr_minimum<-c(0.005,0.002)

for(scexp in 1:length(list_input)){
  
  current_results<-list_input[[scexp]]
  current_results[,1]<-unlist(lapply(strsplit(gsub(current_results[,1],pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[2:5],collapse="-")}))
  
  current_cancer<-names(list_input)[scexp]
  
  vector_order_to_use<-list_order[[scexp]]
  
  current_pvalue<-list_pvalue[scexp]
  
  thr_min<-thr_minimum[scexp]
  
  dir_cancer<-paste("/home/guidantoniomt/datasets/TCGA/harmonized/",current_cancer,"/","GDCdata",sep="")
  
  setwd(dir_cancer)

  mutect_file<-fread(file=grep(dir(),pattern="mutect",value=T))
  
  mutect_file_snvs <- data.frame(mutect_file[,c('Tumor_Sample_Barcode','Hugo_Symbol',
                                  'Chromosome','Start_Position','End_Position',
                                  'Variant_Classification','Variant_Type',
                                  'Reference_Allele','Tumor_Seq_Allele1',
                                  'Tumor_Seq_Allele2')])
  
  # Only interested in point mutations (i.e. VariantType == 'SNP' or Alleles %in% c('A','T','C','G'))
  mutect_file_snvs2 <- mutect_file_snvs[mutect_file_snvs$Reference_Allele %in% c('A','C','T','G') & mutect_file_snvs$Tumor_Seq_Allele2 %in% c('A','C','T','G'),]
  
  mutect_file_snvs3<-mutect_file_snvs2[mutect_file_snvs2$Hugo_Symbol %in% cosmic_genes,]
  
  mutect_file_snvs3$Tumor_Sample_Barcode<-unlist(lapply(strsplit(gsub(mutect_file_snvs3$Tumor_Sample_Barcode,pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[1:4],collapse="-")}))
  
  mutect_file_snvs4 <- merge(current_results,mutect_file_snvs3,by.x="Samples",by.y="Tumor_Sample_Barcode",all.x=F,all.y=F)
  
  number_patients_each_cluster<-as.numeric(table(mutect_file_snvs4$clusters))
  
  freq_mut_for_cluster<-table(mutect_file_snvs4$clusters,mutect_file_snvs4$Hugo_Symbol)
  percent_mut_for_cluster<-apply(freq_mut_for_cluster,1,FUN=function(X){X/sum(X)})
  percent_mut_for_genes<-apply(freq_mut_for_cluster,2,FUN=function(X){X/sum(X)})
  
  resFreq<-vector(mode="list",ncol(percent_mut_for_cluster))
  
  for(pcmut in 1:ncol(percent_mut_for_cluster)){
  
  sorted_genes<-percent_mut_for_cluster[order(percent_mut_for_cluster[,pcmut],decreasing=T),pcmut]
  
  selected_genes<-names(sorted_genes[sorted_genes>thr_min])
  
  resFreq[[pcmut]]<-selected_genes
  
  }
  
  names(resFreq)<-paste("cluster",colnames(percent_mut_for_cluster),sep="")
  
  # 
  # library(VennDetail)
  # ven <- venndetail(resFreq)
  # 
  # 
  #  Raw data for values heatmap
  # 
  rownames(freq_mut_for_cluster)<-paste("cluster",rownames(freq_mut_for_cluster),sep="")
  temp_matrix<-t(freq_mut_for_cluster[,colnames(freq_mut_for_cluster)%in%as.character(unlist(resFreq))])
  values_for_heatmap<-melt(temp_matrix)
  colnames(values_for_heatmap)<-c("genes","cluster","freq")
  
  sum_for_barplot<-data.frame(rowSums(temp_matrix))  
  sum_for_barplot2<-data.frame(genes=rownames(sum_for_barplot),sum=sum_for_barplot)
  colnames(sum_for_barplot2)[2]<-"sum"
  sum_for_barplot2<-sum_for_barplot2[order(sum_for_barplot2[,2],decreasing=T),]
  
  # 
  # data for colors heatmap
  # 
  rownames(percent_mut_for_genes)<-paste("cluster",rownames(percent_mut_for_genes),sep="")
  temp_matrix<-percent_mut_for_genes[,colnames(percent_mut_for_genes)%in%as.character(unlist(resFreq))]
  colors_for_heatmap<-melt(temp_matrix)
  colnames(colors_for_heatmap)<-c("cluster","genes","freq_perc")
  
  input_for_heatmap2<-merge(values_for_heatmap,colors_for_heatmap,by=c("genes","cluster"))
  

  library(forcats)

  #
  # Nested models
  # 
  genes_for_lme<-unique(input_for_heatmap2$genes)
  
  AOVGenes<- function(mutect_file_snvs4,genes_for_lme){
    
      all_pvalue<-NULL
      df_genes_emt_scores<-NULL
      res_AOV<-vector(mode="list",2)
      
      for(glme in genes_for_lme){
        
        emt_gene_with_mut<-unique(data.frame(status=1,mutect_file_snvs4[mutect_file_snvs4$Hugo_Symbol%in%glme,c("Samples","clusters","score_emt")]))
        emt_gene_without_mut<-unique(data.frame(status=0,mutect_file_snvs4[-which(mutect_file_snvs4$Hugo_Symbol%in%glme),c("Samples","clusters","score_emt")]))
        # full_mat<-rbind(emt_gene_with_mut,emt_gene_without_mut)
        # full_mat$status<-as.factor(full_mat$status)
        # full_mat$clusters<-as.factor(full_mat$clusters)
        # 
        # p<-ggboxplot(full_mat, "status", "score_emt", color = "clusters")+labs(title=glme)
        # print(p)
        # 
        #       fullmodel<-lmer(score_emt ~ status + (1|clusters2), data=full_mat)
        #       nullmodel<-lmer(score_emt ~ (1|clusters2), data=full_mat)
        #       pvalue_anova<-anova(fullmodel,nullmodel)["Pr(>Chisq)"][2,]      
        #       
        # pvalue_anova<-wilcox.test(emt_gene_with_mut$clusters,emt_gene_with_mut$score_emt)$p.value
        pvalue_anova <-  as.numeric(unlist(summary(aov(score_emt ~ clusters, data = emt_gene_with_mut)))["Pr(>F)1"])
        
        all_pvalue<-c(all_pvalue,pvalue_anova)
      
        df_genes_emt_scores<-rbind(df_genes_emt_scores,data.frame(gene=glme,emt_gene_with_mut))
      }
      
      names(all_pvalue)<-genes_for_lme
      
      res_AOV[[1]]<-all_pvalue
      res_AOV[[2]]<-df_genes_emt_scores
      
      return(res_AOV)
    }
  
  resAOV<-AOVGenes(mutect_file_snvs4,genes_for_lme)
  if(current_cancer=="LUAD"){
  genes_aov<-p.adjust(resAOV[[1]],"BH")
  }else{
  genes_aov<-resAOV[[1]]
  }
  genes_aov2<-genes_aov[genes_aov<=current_pvalue]
  # genes_aov2<-genes_aov[genes_aov<=0.30]
  
  input_for_heatmap2$genes<-factor(input_for_heatmap2$genes, levels=rev(sum_for_barplot2[,1]))
  
  input_for_heatmap3<-input_for_heatmap2[input_for_heatmap2$genes%in%names(genes_aov2),]
  
  library(VennDetail)
  
  ven <- venndetail(list(genes_Clusters=names(genes_aov2),
                         mes_vs_epi = mut_genes_mes_vs_epi,
                         pemt_vs_epi = mut_genes_pemt_vs_epi,
                         mex_vs_mix = mut_genes_mes_vs_mix))
  
  output_pdf_venn<-paste("venn_for_clusters",current_cancer,".pdf",sep="")
  output_txt_venn<-paste("venn_for_clusters",current_cancer,".txt",sep="")
  
  setwd("/home/guidantoniomt/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")
  
  pdf(output_pdf_venn)
  
  p<-plot(ven, type = "upset")
  print(p)
  
  results2<-result(ven, wide = TRUE)
  
  write.table(results2,file=output_txt_venn,sep="\t",row.names=F,quote=F)
  
  dev.off()
  
  #
  # Nested models
  # 
  
  # g2<-ggplot(input_for_heatmap3, aes(genes, cluster)) +
  #   geom_tile(aes(fill = freq_perc)) + 
  #   geom_text(aes(label = round(freq, 1))) +
  #   scale_fill_gradient(low = "white", high = "steelblue")+theme_minimal()+coord_polar(theta="x")
  
  emt_clusters<-resAOV[[2]][,c(1,3,4,5)]
  emt_clusters2<-aggregate(score_emt~clusters+gene, emt_clusters, median)
  # emt_clusters2$gene<-factor(input_for_heatmap2$genes, levels=rev(sum_for_barplot2[,1]))
  #   
  # emt_clusters2$gene<-intersect(emt_clusters2$gene,factor(input_for_heatmap2$genes, levels=rev(sum_for_barplot2[,1])))
  #   
  # p3<-ggplot(input_for_heatmap3, aes(x = genes, y = cluster,size = freq_perc,color=freq_perc)) +
  #   geom_point() + scale_colour_viridis(option = "plasma",direction = -1)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                                                                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+coord_polar(theta="x")
  
  emt_clusters3<-emt_clusters2
  emt_clusters3$clusters<-paste("cluster",emt_clusters2$clusters,sep="")
  
  test4<-merge(input_for_heatmap3,emt_clusters3,by.x=c("genes","cluster"),by.y=c("gene","clusters"))
  test4$cluster<-factor(test4$cluster,levels=vector_order_to_use)
  
  hypoxia_df<-merge(emt_clusters,hypoxia_results_df2[,c(1:2)],by.x="Samples",by.y="samples")
  
  df_for_boxplot_hypoxia<-unique(hypoxia_df[,c(1,3,5)])
  df_for_boxplot_hypoxia$clusters<-paste("cluster",df_for_boxplot_hypoxia$clusters,sep="")
  df_for_boxplot_hypoxia$clusters<-factor(df_for_boxplot_hypoxia$cluster,levels=vector_order_to_use)
  
  
  hypoxia_clusters2<-aggregate(Buffa_hypoxia_genes~clusters+gene, hypoxia_df, median)
  hypoxia_clusters2$clusters<-paste("cluster",hypoxia_clusters2$clusters,sep="")
  hypoxia_clusters2$clusters<-factor(hypoxia_clusters2$cluster,levels=vector_order_to_use)
  
  output_pdf_hypoxia<-paste("hypoxia_for_clusters",current_cancer,".pdf",sep="")
  print(output_pdf_hypoxia)
  print(getwd())
  
  setwd("/home/guidantoniomt/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")
  
  pdf(output_pdf_hypoxia)
  
  colors2<-brewer.pal(n = length(unique(df_for_boxplot_hypoxia$clusters)), name = "Reds")
  
  pboxplot<-ggviolin(df_for_boxplot_hypoxia, "clusters", "Buffa_hypoxia_genes", fill = "clusters",
           palette = colors2,
           add = "boxplot", add.params = list(fill = "white"))+theme(axis.text.x=element_text(angle = -90, hjust = 0))+  stat_compare_means()
           
  print(pboxplot)
  
  dev.off()
  
    
  test5<-merge(test4,hypoxia_clusters2,by.x=c("genes","cluster"),by.y=c("gene","clusters"))
  test5$cluster<-factor(test5$cluster,levels=vector_order_to_use)

  output_pdf<-paste("circular_heatmap_mutation_emt.",current_cancer,".pdf",sep="")
  
  pdf(output_pdf)
  
  p<-ggplot(test5, aes(y = factor(cluster),
                    x = factor(genes))) +        ## global aes
    geom_tile(aes(fill = freq_perc))+
    scale_fill_gradientn(colours = (brewer.pal(11,"Greys"))) +         ## to get the rect filled
    geom_point(aes(colour = score_emt, 
                   size = score_emt))  +    ## geom_point for circle illusion
    scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")))+
    theme_bw()+
    coord_polar(theta="x")
  
  print(p)
  
  dev.off()
}