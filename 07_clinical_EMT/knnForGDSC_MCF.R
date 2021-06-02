knnForGDSC_MCF<-function(FM,knn_k = 10, num_dim = 25,ann_pseudotime){
  
  # message("Performing dimensionality reduction with PCA")
  # 
  # irlba_res <- prcomp(t(FM))
  # 
  irlba_pca_res <- FM$x[,c(1:num_dim)]
  # rownames(irlba_pca_res)<-colnames(FM)
  
  mcf_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='\\Mock_TGFB_agregated_samples'),]
  xmatrix_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='\\Mock_TGFB_agregated_samples',invert=T),]
  
  all_id<-NULL
  all_pseudotime<-NULL

  for(xmatrix_exp in 1:length(rownames(xmatrix_pca))){
    
    print(xmatrix_exp)
    
    xmatrix_id2<-rownames(xmatrix_pca)[xmatrix_exp]
    all_id<-c(all_id,xmatrix_id2)
    
    xmatrix.mtx<-xmatrix_pca[row.names(xmatrix_pca)==xmatrix_id2,]
    
    res<-get.knnx(mcf_pca,t(data.frame(xmatrix.mtx)), k=knn_k, algorithm="kd_tree")
    
    id_mcf<-gsub(rownames(mcf_pca)[res$nn.index],pattern="\\.",replacement="-")

    mean_pseudotime<-mean(ann_pseudotime[ann_pseudotime$cell %in% id_mcf,"Pseudotime"])
    
    all_pseudotime<-c(all_pseudotime,mean_pseudotime)
    
    }
  
    resFinal<-data.frame(all_id,all_pseudotime)
  
  return(resFinal)
}
