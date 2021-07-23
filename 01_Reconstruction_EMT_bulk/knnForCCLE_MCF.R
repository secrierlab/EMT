knnForCCLE_MCF<-function(FM,knn_k = 10, num_dim = 25,ann_pseudotime,list_ann){
  
  # message("Performing dimensionality reduction with PCA")
  # 
  # irlba_res <- prcomp(t(FM))
  # 
  irlba_pca_res <- FM$x[,c(1:num_dim)]
  # rownames(irlba_pca_res)<-colnames(FM)
  
  mcf_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='\\CCLE_',invert=T),]
  ccle_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='\\CCLE_'),]
  
  all_id<-NULL
  all_pseudotime<-NULL

  for(ccle_exp in 1:length(rownames(ccle_pca))){
    
    print(ccle_exp)
    
    ccle_id2<-rownames(ccle_pca)[ccle_exp]
    all_id<-c(all_id,ccle_id2)
    
    ccle.mtx<-ccle_pca[row.names(ccle_pca)==ccle_id2,]
    
    res<-get.knnx(mcf_pca,t(data.frame(ccle.mtx)), k=knn_k, algorithm="kd_tree")
    
    id_mcf<-gsub(rownames(mcf_pca)[res$nn.index],pattern="\\.",replacement="-")

    mean_pseudotime<-mean(ann_pseudotime[ann_pseudotime$cell %in% id_mcf,"Pseudotime"])
    
    all_pseudotime<-c(all_pseudotime,mean_pseudotime)
    
    }
  
    resFinal<-data.frame(all_id,all_pseudotime)
  
  return(resFinal)
}
