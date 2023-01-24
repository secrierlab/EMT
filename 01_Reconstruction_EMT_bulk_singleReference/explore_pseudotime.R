folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

markers_genes<-read.delim(file="EMT_and_pEMT_markers.txt",header=T)

list_results_clusters_markers<-vector(mode="list",2)

names(list_results_clusters_markers)<-colnames(markers_genes)[c(3:4)]

for(col_comparison in colnames(markers_genes)[c(3:4)]){
  
  print(col_comparison)
  
  markers_emt_met_pemt<-markers_genes[which(markers_genes[,col_comparison] %in% 1),]
  
  genes_metastatic_mock<-combat_mock_filter[combat_mock_filter[,1]%in%markers_emt_met_pemt[,1],1]
  genes_metastatic_tgfb<-combat_tgfb2_filter[combat_tgfb2_filter[,1]%in%markers_emt_met_pemt[,1],1]
  
  #
  # Make the clustering between the gene-expression and the pseudotime with also genes partial EMT
  #
  
  # Here i only check which are the metastatic and epithelial markers that were used for the analysis. However, I will plot all markers available.
  #
  
  common_genes<-intersect(genes_metastatic_mock,genes_metastatic_tgfb)
  
  tcga_met_table2<-t(tcga[which(tcga[,1]%in%common_genes),-1])
  colnames(tcga_met_table2)<-tcga[which(tcga[,1]%in%common_genes),1]
  patientsID<-gsub(make.names(sapply(strsplit(as.character(rownames(tcga_met_table2)),split="\\."),"[[",4),unique=T),pattern="\\.",replacement="_")
  
  #
  # Scale each gene between 0 to 1
  #
  
  normmax<-function(x){(x)/max(x)}
  tcga_met_table2<-apply(tcga_met_table2,2,normmax)
  
  markGenes<-data.frame(patients=patientsID,tcga_met_table2)
  
  #
  # Scale the pseudotime data between 0 to 1
  #
  
  epithelial_cancer<-c("ACC","BLCA","BRCA","CESC","CHOL",
                       "COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC",
                       "LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM",
                       "STAD","THYM","THCA","UCS","UCEC","UVM")
  
  input_cp_tum<-input_cp[which(input_cp$tumors%in%epithelial_cancer),]
  
  print(dim(input_cp_tum))
  
  input_cp_to_scale=input_cp_tum
  
  print(dim(input_cp_to_scale))
  
  input_cp_to_scale$mock<-input_cp_to_scale$mock/100
  input_cp_to_scale$tgfb<-input_cp_to_scale$tgfb/100
  
  #
  # Create the input data to do the clustering
  #
  
  input_cp_scale<-merge(input_cp_to_scale,markGenes,by="patients")[,c("patients","mock","tgfb",common_genes)]
  rownames(input_cp_scale)<-input_cp_scale$patients
  
  library(mclust)
  library(clustvarsel)
  
  print("Run clustering!!!")

  automatic=FALSE
  
  if(automatic!=TRUE){
  
  best_cluster=3
  
  results_mclust<-Mclust(input_cp_scale[,-1],G=best_cluster)
  
  }else{
    
  results_mclust<-Mclust(input_cp_scale[,-1])
  
  }

  samples_clusters<-rownames(input_cp_scale)
  
  dfclusters<-data.frame(patients=samples_clusters,
                         clusters=as.factor(results_mclust$classification))
  
  input_cp_scale<-merge(input_cp_scale,dfclusters,by="patients")
  
  
  require(ComplexHeatmap)
  
  
  df_clustering<-data.frame(patients=input_cp_scale$patients,
                            clusters=input_cp_scale$clusters
  )
  
  n_genes<-intersect(genes_metastatic_mock,genes_metastatic_tgfb)
  
  tcga_met_table<-t(tcga[which(tcga[,1]%in%common_genes),-1])
  colnames(tcga_met_table)<-tcga[which(tcga[,1]%in%common_genes),1]
  
  tcga_metastatic_genes<-data.frame(patients=as.character(colnames(tcga)[-1]),tcga_met_table)
  tcga_metastatic_genes[,1]<-gsub(make.names(sapply(strsplit(as.character(tcga_metastatic_genes[,1]),split="\\."),"[[",4),unique=T),pattern="\\.",replacement="_")
  
  z_score<-function(x){(x-mean(x))/sd(x)}
  tcga_metastatic_genes[,-1]<-apply(tcga_metastatic_genes[,-1],2,z_score)
  
  dfclustering2<-merge(df_clustering,input_cp_tum[,c("patients","tumors","tumor_stage")],by="patients")
  
  input_cp_scale2<-merge(dfclustering2[,which(colnames(dfclustering2)%in% c("patients","clusters","tumors","tumor_stage"))],tcga_metastatic_genes,by="patients")
  
  annotation_colors=list(tumors=c("ACC"="#996b2f","BLCA"="#724fb3","BRCA"="#ba449c",
                                  "CESC"="#adba3d","CHOL"="#8f2e43","COAD"="#191970",
                                  "ESCA"="#5060ac","HNSC"="#65b06f","KICH"="#B266FF",
                                  "KIRC"="#969841","KIRP"="#bf6d29","LIHC"="#d1972c",
                                  "LUAD"="#36dee6","LUSC"="#cc4b3e","OV"="#d97eb8",
                                  "PAAD"="#c43d7e","PRAD"="#d077d8","READ"="#FFB266","SKCM"="#452c7c",
                                  "STAD"="#43c29e","THYM"="#de6991","THCA"="#757ee8",
                                  "UCS"="#88321a","UCEC"="#75b550","UVM"="#8f2931"),
                         tumor_stage=c("low_stage"="green3","high_stage"="orange2","metastatic"="red2","not_reported"="gainsboro"),
                         clusters=c("1"="blue2","2"="green2","3"="orange2","4"="red2"))
  
  
  row_colors=list(type=c("Mesenchymal_marker"="firebrick4","Epithelial_marker"="dodgerblue3","pEMT"="orange2"))
  
  ann_samples<-input_cp_scale2[,c(2:4)]
  
  ann_samples$tumor_stage<-gsub(ann_samples$tumor_stage,pattern=" ",replacement="_")
  
  ann_genes<-markers_emt_met_pemt[markers_emt_met_pemt[,1]%in%colnames(tcga_metastatic_genes),]
  colnames(ann_genes)<-c("genes","type")
  
  col_ann<-HeatmapAnnotation(df=ann_samples,col=annotation_colors)
  
  row_colors=list(type=c("Mesenchymal_marker"="firebrick4","Epithelial_marker"="dodgerblue3","pEMT"="orange2"))
  
  max_col<-ncol(input_cp_scale2)
  
  order_samples_colors_genes<-match(colnames(input_cp_scale2[,c(5:max_col)]),ann_genes[,1])
  
  ann_genes2<-data.frame(ann_genes[order_samples_colors_genes,2])
  colnames(ann_genes2)<-"type"
  
  row_ann<-rowAnnotation(df=ann_genes2,col=row_colors)
  
  output<-paste(paste("heatmaps_clusters_from_pseudospace","Mclust","UNION_CLUSTERS","k",best_cluster,col_comparison,sep="."),paste("epithelialcancer.pdf",sep=""),sep=".")
  
  order_columns=factor(ann_samples$clusters,levels=c(1:best_cluster))
  
  pdf(output,width=10)
  ha<-Heatmap(t(input_cp_scale2[,c(5:max_col)]),
              column_split=order_columns,
              left_annotation = row_ann,
              top_annotation = col_ann,
              border=T,
              cluster_column_slices=F)
  print(ha)
  dev.off()
  
  output<-paste(paste("heatmaps_clusters_from_pseudospace","Mclust","UNION_CLUSTERS","k",best_cluster,col_comparison,sep="."),paste("epithelialcancer.KMEANS.pdf",sep=""),sep=".")
  
  order_columns=factor(ann_samples$clusters,levels=c(1:best_cluster))
  
  pdf(output,width=10)
  ha<-Heatmap(t(input_cp_scale2[,c(5:max_col)]),
              column_split=order_columns,
              left_annotation = row_ann,
              top_annotation = col_ann,
              border=T,
              cluster_column_slices=F,
              row_km=3)
  print(ha)
  dev.off()
  
  #
  # Step1: Plot the results of clustering on the pseudospace
  #
  
  pdf(paste('KNN_projection_TCGA_to_MCF10A_mock_vs_tgfb_mclust_scaling',"_k_",best_cluster,col_comparison,'.epithelialcancer.pdf',sep=""))
  
  knn_df_tcga_mock2<-knn_df_tcga_mock[which(knn_df_tcga_mock$tumors%in%epithelial_cancer),]
  knn_df_tcga_tgfb2<-knn_df_tcga_tgfb[which(knn_df_tcga_tgfb$tumors%in%epithelial_cancer),]
  
  mock_scale<-knn_df_tcga_mock2
  tgfb_scale<-knn_df_tcga_tgfb2
  
  mock_scale$pseudospace<-knn_df_tcga_mock2$pseudospace/max(knn_df_tcga_mock2$pseudospace)
  tgfb_scale$pseudospace<-knn_df_tcga_tgfb2$pseudospace/max(knn_df_tcga_tgfb2$pseudospace)
  
  mock<-ggplot(mock_scale, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")
  
  tgfb<-ggplot(tgfb_scale, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C")+ coord_flip()
  
  if(best_cluster==3){

  sp <- ggplot(input_cp_scale, aes(x=mock, y=tgfb,col=clusters)) +geom_point(alpha = 5/10)+scale_color_manual(breaks = c("3", "2", "1"), values=c("red2", "orange2", "blue2"))

  }

  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
    )
  
  grid.arrange(mock,blankPlot,sp, tgfb,
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  
  dev.off()
  
  #
  # Step2: Plot the results of clustering on the pseudospace facetwrap
  #
  
  pdf(paste('KNN_projection_TCGA_to_MCF10A_mock_vs_tgfb_mclust_scaling_facetwrap',"_k_",best_cluster,col_comparison,'.epithelialcancer.pdf',sep=""))
  
  knn_df_tcga_mock2<-knn_df_tcga_mock[which(knn_df_tcga_mock$tumors%in%epithelial_cancer),]
  knn_df_tcga_tgfb2<-knn_df_tcga_tgfb[which(knn_df_tcga_tgfb$tumors%in%epithelial_cancer),]
  
  mock_scale<-knn_df_tcga_mock2
  tgfb_scale<-knn_df_tcga_tgfb2
  
  mock_scale$pseudospace<-knn_df_tcga_mock2$pseudospace/max(knn_df_tcga_mock2$pseudospace)
  tgfb_scale$pseudospace<-knn_df_tcga_tgfb2$pseudospace/max(knn_df_tcga_tgfb2$pseudospace)
  
  mock<-ggplot(mock_scale, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#0075F2")
  
  tgfb<-ggplot(tgfb_scale, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#70163C")+ coord_flip()
  
  sp <- ggplot(input_cp_scale, aes(x=mock, y=tgfb,col=clusters)) +geom_point(alpha = 5/10)+facet_wrap(~clusters)
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
    )
  
  grid.arrange(mock,blankPlot,sp, tgfb,
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  
  dev.off()
  
  #
  # Step3: Plot the the expression of markers genes in the pseudospace
  #
  
  n_genes<-intersect(genes_metastatic_mock,genes_metastatic_tgfb)
  
  tcga_met_table<-t(tcga[which(tcga[,1]%in%common_genes),-1])
  colnames(tcga_met_table)<-tcga[which(tcga[,1]%in%common_genes),1]
  
  tcga_metastatic_genes<-data.frame(patients=as.character(colnames(tcga)[-1]),tcga_met_table)
  tcga_metastatic_genes[,1]<-gsub(make.names(sapply(strsplit(as.character(tcga_metastatic_genes[,1]),split="\\."),"[[",4),unique=T),pattern="\\.",replacement="_")
  
  z_score<-function(x){(x-mean(x))/sd(x)}
  tcga_metastatic_genes[,-1]<-apply(tcga_metastatic_genes[,-1],2,z_score)
  
  input_cp_scale2<-merge(input_cp_scale[,c("patients","mock","tgfb","clusters")],tcga_metastatic_genes,by="patients")
  
  facet=TRUE
  
  for(metgenes in common_genes){
    
    output_file<-paste(paste("KNN_projection_TCGA_to_MCF10A_mock_vs_tgfb_mclust_scaling_facet","_k_",best_cluster,"_",col_comparison,"_",metgenes,sep=""),".epithelialcancer.pdf",sep="")
    
    pdf(output_file)
    
    input_for_plot<-input_cp_scale2[,c("mock","tgfb",metgenes,"clusters")]
    colnames(input_for_plot)[3]<-"gene"
    
    input_for_plot$clusters<-as.factor(input_for_plot$clusters)
    
    if(facet==TRUE){
      
      sp <- ggplot(input_for_plot,
                   aes(x=mock, y=tgfb,col=gene,shape=clusters))+geom_point(alpha = 5/10)+scale_colour_gradient2(low = "blue", high = "red")+facet_wrap(~clusters)
    }else{
      
      sp <- ggplot(input_for_plot, aes(x=mock, y=tgfb,col=gene,shape=clusters))+geom_point(alpha = 5/10)+scale_colour_gradient2(low = "blue", high = "red")
      
    }
    
    blankPlot <- ggplot()+geom_blank(aes(1,1))+
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
      )+labs(title=metgenes)
    
    grid.arrange(mock,blankPlot,sp, tgfb,
                 ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
    
    dev.off()
    
  }
  
  require(reshape2)
  require(ggpubr)
  
  input_for_boxplot<-melt(input_cp_scale2[,c(4:ncol(input_cp_scale2))])
  
  dfcomp<-expand.grid(c(1:best_cluster),c(1:best_cluster))
  dfcomp2<-dfcomp[which(dfcomp[,1]!=dfcomp[,2]),]
  dfcomp3<-dfcomp2[which(duplicated(rowSums(dfcomp2))%in%FALSE),]
  
  my_comparisons <-  lapply(split(dfcomp3,seq(nrow(dfcomp3))),as.character) 
  
  pdf(paste('KNN_projection_TCGA_to_MCF10A_mock_vs_tgfb_mclust_boxplot',"_k_",best_cluster,col_comparison,'.epithelialcancer.pdf',sep=""))

  if(best_cluster>3){
  palette_to_use = c("green2","blue2", "orange2", "red2")
  }else{
  palette_to_use = c("orange2", "blue2", "red2")
  }

  p<-ggviolin(input_for_boxplot, x = "clusters", y = "value", fill = "clusters",
              palette = palette_to_use,                                                          
              add = "boxplot", add.params = list(fill = "white"))+
              stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+                             
	      stat_summary(fun.y=median, geom="line", aes(group=1))
 
  p<-facet(p, facet.by = "variable")
  print(p)
  dev.off()
  
  #
  # Step4: Create heatmaps for each marker genes and clusters
  #
  
  
  annotation_colors=list(tumors=c("ACC"="#996b2f","BLCA"="#724fb3","BRCA"="#ba449c",
                                  "CESC"="#adba3d","CHOL"="#8f2e43","COAD"="#191970",
                                  "ESCA"="#5060ac","HNSC"="#65b06f","KICH"="#B266FF",
                                  "KIRC"="#969841","KIRP"="#bf6d29","LIHC"="#d1972c",
                                  "LUAD"="#36dee6","LUSC"="#cc4b3e","OV"="#d97eb8",
                                  "PAAD"="#c43d7e","PRAD"="#d077d8","READ"="#FFB266","SKCM"="#452c7c",
                                  "STAD"="#43c29e","THYM"="#de6991","THCA"="#757ee8",
                                  "UCS"="#88321a","UCEC"="#75b550","UVM"="#8f2931"),
                         tumor_stage=c("low_stage"="green3","high_stage"="orange2","metastatic"="red2","not_reported"="gainsboro"),
                         status_genes=c("Mesenchymal_marker"="firebrick4",
                                        "Epithelial_marker"="dodgerblue3",
                                        "pEMT"="orange2")
  )
  
  
  ann_genes<-markers_emt_met_pemt[markers_emt_met_pemt[,1]%in%colnames(tcga_metastatic_genes),]
  colnames(ann_genes)<-c("genes","type")
  
  input_cp_scale3<-merge(input_cp_scale2,input_cp_tum[,c("patients","tumors","tumor_stage")],by="patients")
  input_cp_scale3$tumors<-as.character(input_cp_scale3$tumors)
  
  for(myclu in unique(input_cp_scale3$clusters)){
    
    current_patients<-input_cp_scale3[input_cp_scale3$clusters%in%myclu,1]
    
    input_met_heatmap<-t(tcga_metastatic_genes[which(tcga_metastatic_genes[,1]%in%current_patients),-1])
    
    input_met_heatmap_sort<-input_met_heatmap[match(ann_genes[,1],rownames(input_met_heatmap)),]
    colnames(input_met_heatmap_sort)<-gsub(make.names(sapply(strsplit(colnames(input_met_heatmap_sort),split="\\."),"[[",4),unique=T),pattern="\\.",replacement="_")
    
    ann_samples<-input_cp_scale3[,c("tumors","tumor_stage")]
    rownames(ann_samples)<-input_cp_scale3$patients
    ann_samples<-ann_samples[which(rownames(ann_samples)%in%current_patients),]
    ann_samples$tumor_stage<-gsub(ann_samples$tumor_stage,pattern="not reported",replacement="not_reported")
    
    ann_genes2<-data.frame(ann_genes[,2])
    colnames(ann_genes2)<-"status_genes"
    
    rownames(ann_genes2)<-ann_genes[,1]
    
    require(pheatmap)
    
    output<-paste(paste("heatmaps_clusters_from_pseudospace_mclust","k",best_cluster,col_comparison,sep="."),paste(myclu,".epithelialcancer.pdf",sep=""),sep=".")
    
    pdf(output,width=10,height=10)
    
    pheatmap(input_met_heatmap_sort,
             annotation_row = ann_genes2,
             annotation_col = ann_samples,
             show_colnames=FALSE,
             scale="none",
             annotation_colors=annotation_colors)
    
    dev.off()
  }
  
  ann_genes<-markers_emt_met_pemt[markers_emt_met_pemt[,1]%in%colnames(tcga_metastatic_genes),]
  colnames(ann_genes)<-c("genes","type")
  
  canonical_markers<-c("VIM","CDH2","FOXC2","SNAI1","SNAI2","TWIST1","FN1","ITGB6","MMP2","MMP3","MMP9","SOX10","GCS","CDH1","DSP","OCLN","ZEB1","ZEB2",                                                                            
                       "TWIST2","CRB3")
  
  ann_genes<-ann_genes[ann_genes[,1]%in%canonical_markers,]
  colnames(ann_genes)<-c("genes","type")
  
  #
  # Step4: Create heatmaps for each marker genes and canonical markers genes
  #
  
  
  for(myclu in unique(input_cp_scale2$clusters)){
    
    current_patients<-input_cp_scale2[input_cp_scale2$clusters%in%myclu,1]
    
    input_met_heatmap<-t(tcga_metastatic_genes[which(tcga_metastatic_genes[,1]%in%current_patients),-1])
    
    input_met_heatmap<-input_met_heatmap[which(rownames(input_met_heatmap)%in%canonical_markers),]
    
    input_met_heatmap_sort<-input_met_heatmap[match(ann_genes[,1],rownames(input_met_heatmap)),]
    colnames(input_met_heatmap_sort)<-gsub(make.names(sapply(strsplit(colnames(input_met_heatmap_sort),split="\\."),"[[",4),unique=T),pattern="\\.",replacement="_")
    
    ann_genes2<-data.frame(ann_genes[,2])
    colnames(ann_genes2)<-"status_genes"
    
    ann_samples<-input_cp_scale3[,c("tumors","tumor_stage")]
    rownames(ann_samples)<-input_cp_scale3$patients
    ann_samples<-ann_samples[which(rownames(ann_samples)%in%current_patients),]
    ann_samples$tumor_stage<-gsub(ann_samples$tumor_stage,pattern="not reported",replacement="not_reported")
    
    rownames(ann_genes2)<-ann_genes[,1]
    
    require(pheatmap)
    
    output<-paste(paste("heatmaps_clusters_from_pseudospace_mclust.canonical_emt","k",best_cluster,col_comparison,sep="."),paste(myclu,".epithelialcancer.pdf",sep=""),sep=".")
    
    pdf(output,width=10,height=10)
    
    pheatmap(input_met_heatmap_sort,
             annotation_row = ann_genes2,
             annotation_col = ann_samples,
             show_colnames=FALSE,
             scale="none",
             annotation_colors=annotation_colors)
    
    dev.off()
  }
  
  #
  #
  #
  input_cp_save<-merge(input_cp_scale2[,c("patients","clusters")],input_cp,by="patients")
  
  
  list_results_clusters_markers[[col_comparison]]<-input_cp_save
  
}


