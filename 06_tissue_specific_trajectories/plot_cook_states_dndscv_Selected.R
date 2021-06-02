library(VennDetail)
library(data.table)
library("readxl")
library(forcats)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

list_tum<-c("LUAD","BRCA")

for(tum in list_tum){
  
  setwd("/home/guidantoniomt/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")
  
  # Create a vector with the list of files to use, in order by state (1:5)
  files_current_tum<-paste(paste(tum,"dndscv_siggenes.0.1",1:5,sep="."),".txt",sep="")
  
  if(tum=="LUAD"){
    
  clusters_df<-data.frame(read_excel("LUAD_groups_for_Lucie.xlsx"))
  
  load("/home/guidantoniomt/pseudospace/res_multiple_pseudospace/A549_mapped_seurat_correct_all_timecourse.RData")
  
  df_scores_EMT2<-df_scores_EMT[,c(1:2)]
  
  current_results<-merge(clusters_df,df_scores_EMT2,by.x="Samples",by.y="Samples")
  
  vector_order_to_use<-c("cluster3","cluster4","cluster5","cluster2","cluster1")

  genes_degs<-c("ARID1A","ARID2","CDKN2A","CMTR2","EGFR","KEAP1","KRAS","MGA","NF1","RB1","RBM10","REG1B","REG3A","SETD2","SMARCA4","STK11","TP53","ZFP36L1","ZIC1")

  order_degs<-c("TP53","STK11","KRAS","EGFR","SETD2","RBM10","CDKN2A","RB1","MGA","ZIC1","NF1","REG3A","ARID2","ZFP36L1","ARID1A","CMTR2","SMARCA4","REG1B","KEAP1")

  }else{
    
  clusters_df<-data.frame(read_excel("data_for_lucie_BRCA_PRAD_OV.xlsx","BRCA"))
  
  load("/home/guidantoniomt/pseudospace/res_multiple_pseudospace/MCF7_mapped_seurat_correct_all_timecourse.RData")
  df_scores_EMT2<-df_scores_EMT[,c(1:2)]
  
  current_results<-merge(clusters_df,df_scores_EMT2,by.x="Samples",by.y="Samples")
  
  vector_order_to_use<-c("cluster4","cluster2","cluster1","cluster5","cluster3")
  
  genes_degs<-c("GATA3","MAP3K1","PIK3CA","TP53","CDH1","PTEN","RUNX1","TBX3","AKT1","SHISA4","NF1","CBFB","KMT2C","MAP2K4","PIK3R1","ARID1A","FOXA1","FBXW7","RB1","NCOR1","GPS2")
  order_degs<-c("GATA3","MAP3K1","PIK3CA","TP53","CDH1","PTEN","RUNX1","TBX3","KMT2C","MAP2K4","PIK3R1","CBFB","NCOR1","ARID1A","RB1","AKT1","SHISA4","NF1","MAPK2KA","FOXA1","FBXW7","GPS2")


  }
  
  res_global<-data.frame()
  
  for(i in 1:length(files_current_tum)){
  
  # The order of the state correspond with the order of the files
  current_tab<-read.delim(file=files_current_tum[i])
  tab_res_dndscv<-data.frame(state=i,current_tab)
  
  res_global<-rbind(res_global,tab_res_dndscv)
  
  }
  
  res_global<-res_global[grep(res_global$gene_name,pattern="^OR",invert=T),]
    
  results_dndscv<-lapply(split(res_global,f=res_global$state),FUN=function(X){X[,3]})
  
  output_venn<-paste(paste(tum,"dndscv_Cook_states","UpSet",sep="."),".txt",sep="")
  
  ven <- venndetail(results_dndscv)
  
  res_ven<-result(ven, wide = TRUE)
  
  #write.table(res_ven,file=output_venn,sep="\t",quote=F,row.names=F)
  
  gene_list<-intersect(unlist(results_dndscv),genes_degs)
    
  #
  # Prepare the data for the circular charts
  #
  
  dir_cancer<-paste("/home/guidantoniomt/datasets/TCGA/harmonized/",tum,"/","GDCdata",sep="")
  
  setwd(dir_cancer)
  
  mutect_file<-fread(file=grep(dir(),pattern="mutect",value=T))
  
  mutect_file_snvs <- data.frame(mutect_file[,c('Tumor_Sample_Barcode','Hugo_Symbol',
                                                'Chromosome','Start_Position','End_Position',
                                                'Variant_Classification','Variant_Type',
                                                'Reference_Allele','Tumor_Seq_Allele1',
                                                'Tumor_Seq_Allele2')])
  
  # Only interested in point mutations (i.e. VariantType == 'SNP' or Alleles %in% c('A','T','C','G'))
  mutect_file_snvs2 <- mutect_file_snvs[mutect_file_snvs$Reference_Allele %in% c('A','C','T','G') & mutect_file_snvs$Tumor_Seq_Allele2 %in% c('A','C','T','G'),]
  
  mutect_file_snvs3<-mutect_file_snvs2[mutect_file_snvs2$Hugo_Symbol %in% gene_list,]
  
  mutect_file_snvs3$Tumor_Sample_Barcode<-unlist(lapply(strsplit(gsub(mutect_file_snvs3$Tumor_Sample_Barcode,pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[1:4],collapse="-")}))
  
  current_results[,1]<-unlist(lapply(strsplit(gsub(current_results[,1],pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[2:5],collapse="-")}))
  
  mutect_file_snvs4 <- merge(current_results,mutect_file_snvs3,by.x="Samples",by.y="Tumor_Sample_Barcode",all.x=F,all.y=F)
  
  number_patients_each_cluster<-as.numeric(table(mutect_file_snvs4$clusters))
  

  freq_mut_for_cluster<-table(mutect_file_snvs4$clusters,mutect_file_snvs4$Hugo_Symbol)
  percent_mut_for_cluster<-apply(freq_mut_for_cluster,1,FUN=function(X){X/sum(X)})
  percent_mut_for_genes<-apply(freq_mut_for_cluster,2,FUN=function(X){X/sum(X)})
  
  #setwd("/home/guidantoniomt/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")
  #write.table(freq_mut_for_cluster,file=paste(tum,"raw_nsamples_modules_genes.txt",sep="_"),row.names=T,sep="\t")

  #check number of patients with mutations in one gene
  #apply(freq_mut_for_cluster,2,sum)
  #length(unique(mutect_file_snvs3[mutect_file_snvs3$Hugo_Symbol=="TP53",1]))
  
  resFreq<-vector(mode="list",ncol(percent_mut_for_cluster))
  
  for(pcmut in 1:ncol(percent_mut_for_cluster)){
    
    sorted_genes<-percent_mut_for_cluster[order(percent_mut_for_cluster[,pcmut],decreasing=T),pcmut]
    
    selected_genes<-names(sorted_genes)
    
    resFreq[[pcmut]]<-selected_genes
    
  }
  
  names(resFreq)<-paste("cluster",colnames(percent_mut_for_cluster),sep="")
  
  #
  #  Raw data for values heatmap
  # 
  
  rownames(freq_mut_for_cluster)<-paste("cluster",rownames(freq_mut_for_cluster),sep="")
  temp_matrix<-t(freq_mut_for_cluster[,colnames(freq_mut_for_cluster)%in%as.character(unlist(resFreq))])
  values_for_heatmap<-melt(temp_matrix)
  colnames(values_for_heatmap)<-c("genes","cluster","freq")
  
  # 
  # data for colors heatmap
  # 
  
  rownames(percent_mut_for_genes)<-paste("cluster",rownames(percent_mut_for_genes),sep="")
  temp_matrix<-percent_mut_for_genes[,colnames(percent_mut_for_genes)%in%as.character(unlist(resFreq))]
  colors_for_heatmap<-melt(temp_matrix)
  colnames(colors_for_heatmap)<-c("cluster","genes","freq_perc")
  
  input_for_heatmap3<-merge(values_for_heatmap,colors_for_heatmap,by=c("genes","cluster"))

  emt_clusters<-mutect_file_snvs4[mutect_file_snvs4[,5]%in%gene_list,c(1,2,4,5)]
  emt_clusters2<-aggregate(score_emt~clusters+Hugo_Symbol, emt_clusters, median)
  emt_clusters3<-emt_clusters2
  
  emt_clusters3$clusters<-paste("cluster",emt_clusters2$clusters,sep="")
  
  test4<-merge(input_for_heatmap3,emt_clusters3,by.x=c("genes","cluster"),by.y=c("Hugo_Symbol","clusters"))
  test4$cluster<-factor(test4$cluster,levels=vector_order_to_use)
  
  sign_genes_for_clusters<-res_global[,c(1,3,21)]
  sign_genes_for_clusters[,1]<-paste("cluster",sign_genes_for_clusters[,1],sep="")
  
  #all.x = T, i keep the clusters with the % but they can be also not significant.
  merge5<-merge(test4,sign_genes_for_clusters,by.x=c("genes","cluster"),by.y=c("gene_name","state"),all.x=T)
  
  merge5$qglobal_cv<- -log(merge5$qglobal_cv+1,10)
   
  merge5$genes<-factor(merge5$genes,levels=order_degs)

  output_pdf<-paste(paste("circular_heatmap_mutation_emt_dndscv_Cook_states",tum,sep="."),".SELECTED.pdf",sep="")
 
  # links_genes<-melt(res_ven[,c(1:6)])[,c(1,3)]
  # links_genes$links_to_use<-paste(links_genes[,1],links_genes[,2],sep="_")
  # links_genes[grep(links_genes$links_to_use,pattern="_0"),3]<-""
  # 
  # test5<-merge(test4,links_genes,by.x="genes",by.y="Detail")
  setwd("/home/guidantoniomt/pseudospace/explore_cook_sort_pseudotimes/find_drivers_events")
  
  pdf(output_pdf)
  
  p<-ggplot(merge5, aes(y = factor(cluster),
                       x = factor(genes)))+        ## global aes
    geom_tile(aes(fill = freq_perc))+
    scale_fill_gradientn(colours = (brewer.pal(11,"Greys"))) +         ## to get the rect filled
    geom_point(aes(colour = score_emt, 
                   size = qglobal_cv))+    ## geom_point for circle illusion
    scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")))+
    theme_bw()+
    coord_polar(theta="x")
  
  print(p)
  
  dev.off()
  
  #Integrate with the gene-expression data
  load("/home/guidantoniomt/pseudospace/input_pseudospace/TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")

  columns_to_use<-grep(colnames(TCGA_GEXP_ALL),pattern=tum,value=T)
  
  TCGA_GEXP_ALL_RID<-TCGA_GEXP_ALL[which(TCGA_GEXP_ALL[,1]%in%gene_list),]

  matlog<-log(TCGA_GEXP_ALL_RID[,colnames(TCGA_GEXP_ALL_RID)%in%columns_to_use]+1,2)
  colnames(matlog)<-unlist(lapply(strsplit(gsub(colnames(matlog),pattern="\\.",replacement="-"),split="-"),FUN=function(X){paste(X[2:5],collapse="-")}))

  #zscore<-function(x){(x-mean(x))/sd(x)}

  #matlog_zscore<-t(apply(matlog,1,zscore))

  #exp_matrix_tum<-data.frame(gene_symbol=TCGA_GEXP_ALL_RID[,1],matlog_zscore)
  
  exp_matrix_tum<-data.frame(gene_symbol=TCGA_GEXP_ALL_RID[,1],matlog)

  exp_matrix_tum_melt<-melt(exp_matrix_tum)
  
  exp_matrix_tum_melt$variable<-gsub(exp_matrix_tum_melt$variable,pattern="\\.",replacement="-")

  #merge with the clusters of gene-expression
  exp_matrix_tum_melt2<-merge(exp_matrix_tum_melt,current_results,by.x="variable",by.y="Samples")

  exp_matrix_tum_melt2$clusters<-factor(exp_matrix_tum_melt2$clusters,level=gsub(vector_order_to_use,pattern="cluster",replacement=""))

  list_genes<-unique(exp_matrix_tum_melt2$gene_symbol)
 
  #pdf(paste(tum,"expression_markers_clusters.pdf",sep="_"),width=10)
  
  #mat_rid<-exp_matrix_tum_melt2[,c(2,3)]

  #genes_mean<-unlist(lapply(split(mat_rid,mat_rid[,1]),FUN=function(X){mean(X[,2])}))
 
  #genes_to_save<- names(genes_mean[genes_mean>mean(genes_mean)])
   
  #exp_matrix_tum_melt2$status_deg<-rep("down",nrow(exp_matrix_tum_melt2))
   
  #exp_matrix_tum_melt2[which(exp_matrix_tum_melt2$gene_symbol %in% genes_to_save),"status_deg"]<-"up"

  # bp1<-ggboxplot(exp_matrix_tum_melt2,
  #            x="gene_symbol",
  #           y="value",
  #            add.params = list(fill = "white"),
  #            font.label="bold",
  #	      fill="status_deg")+stat_compare_means()+labs(title=tum)+ geom_hline(yintercept=mean(genes_mean), linetype="dashed", color = "red")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ylab("Expression values log(X+1,2)"
  #	      )
  #print(bp1)

  #plot_gene = FALSE

  #if(plot_gene!=FALSE){

  #for(i in list_genes){
 
  #print(i)

  #submat<-exp_matrix_tum_melt2[exp_matrix_tum_melt2$gene_symbol%in%i,]
  
  #if(isTRUE(max(submat$value)<6)){

  #bp2<-ggviolin(submat,
  #            x="clusters",
  #            y="value",
  #            add="boxplot",
  #            add.params = list(fill = "white"),
  #            font.label="bold",
  #            fill="clusters")+stat_compare_means()+labs(title=i)

  #print(bp2)
  
  #}else{

  #bp2<-ggviolin(submat,
  #           x="clusters",
  #           y="value",
  #	     add="boxplot",
  #	     add.params = list(fill = "white"),
  #	     font.label="bold",
  #	     fill="clusters",
  #	     ylim=c(min(submat$value),1.5))+stat_compare_means()+labs(title=i)

  #print(bp2)

  #}
  
  #}

  #}

  #dev.off()

}
