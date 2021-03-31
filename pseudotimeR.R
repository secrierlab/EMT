#
# kNN method 
#

knnProject<-function(FM,knn_k = 10, num_dim = 25,...){
  
  message("Performing dimensionality reduction with PCA")
  
  irlba_res <- prcomp(t(FM[,-1]))
  
  irlba_pca_res <- irlba_res$x[,c(1:num_dim)]
  rownames(irlba_pca_res)<-colnames(FM[,-1])
  
  bulk_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='DO'),]
  scrna_pca<-irlba_pca_res[grep(rownames(irlba_pca_res),pattern='DO',invert=T),]
  
  knn.list <- list()
  
  for(patient in 1:length(rownames(bulk_pca))){
    print(patient)
    
    patient2<-rownames(bulk_pca)[patient]
    
    bulk.mtx<-bulk_pca[row.names(bulk_pca)==patient2,]
    
    res<-get.knnx(scrna_pca, t(data.frame(bulk.mtx)), k=knn_k, algorithm="kd_tree")
    
    knn.list[[patient2]]<-res
    row.names(knn.list[[patient2]]$nn.index)<-patient2
    colnames(knn.list[[patient2]]$nn.index)<-rownames(scrna_pca)[res$nn.index]
    row.names(knn.list[[patient2]]$nn.dist) <- patient2
    
  }
  return(knn.list)
  
}

#
# Combat 
#

prepareCombat<-function(dat,batch){        
  
  dat <- as.matrix(dat)
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(dat[, batch == batch_level], 1,
                         function(x) {
                           var(x) == 0
                         })))
    }else {
      
      return(which(rep(1, 3) == 2))
      
    }       
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n",
                length(zero.rows)))
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
}


knnWithHMM<-function(matrix_bulk,matrix_scrnaseq,pseudotime_scrnaseq,file_markers="/home/guidantoniomt/pseudospace/HMM/EMT_and_pEMT_markers.txt",hmm_states=3,...){

  matrix_for_combat<-merge(matrix_bulk,matrix_scrnaseq[,-1],by.x='hgnc_symbol')

  batch<-factor(c(rep(1,ncol(matrix_for_combat[,-c(1:2)])),rep(2,ncol(matrix_scrnaseq[,-c(1:2)]))))
  
  print("Step 1: Combat")
  
  matrix_for_combat2<-as.data.frame(setDT(matrix_for_combat[,-2])[, lapply(.SD, mean), by = hgnc_symbol])
  rownames(matrix_for_combat2)<-matrix_for_combat2[,1]
  combat_temp<-prepareCombat(dat=matrix_for_combat2[,-1],batch)
  combat_temp2<- ComBat(dat=combat_temp, batch=batch)
  combat_temp3<-data.frame(genes=rownames(combat_temp2),combat_temp2)
  combat_temp2[,1]<-as.character(combat_temp2[,1])
  combat_temp2[combat_temp2<0]<-0
  
  print("Step 2: Remove low variable genes")
  
  vardata<-apply(combat_temp2[,-1],1,var)
  combat_temp2$vardata<-vardata
  combat_temp2<-combat_temp2[order(combat_temp2$vardata,decreasing=T),]
  combat_filter<- combat_temp2[1:round((nrow(combat_temp2)*30)/100),-ncol(combat_temp2)]
  
  print("Step 3:  Principal component analysis")
  
  irlba_res <- prcomp(t(combat_mock_filter[,-1]))
  
  print("Step 4:  Run kNN")
  
  knn_results<-knnProject(combat_mock_filter,num_dim=25)
  
  knn_pseudospace_final <- NULL
  
  for(patient in names(knn_mock_pseudospace)){
    
    scRNAseq<-gsub(colnames(knn_mock_pseudospace[[patient]]$nn.index),pattern='\\.',replacement='-')
    knn_pseudospace<-mean(pseudotime_scrnaseq[pseudotime_scrnaseq[,1]%in%scRNAseq,'Pseudotime'])
    knn_pseudospace_final<-c(knn_pseudospace_final,knn_pseudospace)
    
  }
  
  knn_pseudospace_final_withID<-data.frame(patients=names(knn_mock_pseudospace),pseudospace=knn_pseudospace_final)
  

  pseudospace_input_sort<-knn_pseudospace_final_withID[order(knn_pseudospace_final_withID$pseudospace,decreasing=F),]
  
  genes_to_segment_df<-read.table(file=file_markers,header=T)
  markers_genes<-genes_to_segment_df[,1]
  
  matrix_bulk<-as.data.frame(setDT(matrix_bulk[,-1])[, lapply(.SD, mean), by = hgnc_symbol])
  rownames(matrix_bulk)<-matrix_bulk[,1]
  
  matrix_for_hmm<-matrix_bulk[matrix_bulk[,1]%in%markers_genes,colnames(matrix_for_hmm)[-1]]
  matrix_for_hmm<-matrix_for_hmm[,match(pseudospace_input_sort$patients,colnames(matrix_bulk))]
  
  print("Step 5: LASSO")
  
  input_glmnet<- t(matrix_for_hmm)
  
  set.seed(12345)
  
  cvfit <-cv.glmnet(x = input_glmnet, y = pseudospace_input_sort$pseudospace, alpha = 0)
  
  dfcoef<-as.matrix(coef(cvfit,s= "lambda.min"))
  genes_to_use<- names(dfcoef[dfcoef[,1]!=0,])[-1]
  
  print("Step 6:  Run HMM")
  
  list_HMM<- lapply(as.list(paste(colnames(input_glmnet[,genes_to_use]),"~1")),as.formula)
  
  list_families<-vector(mode="list",length(genes_to_use))
  
  for(lf in 1:length(genes_to_use)){
    
    list_families[[lf]]<-gaussian()
  }
  
  nstates = hmm_states
  
  input_HMM<-data.frame(time=pseudospace_input_sort$pseudospace,input_glmnet[,genes_to_use])
  
  HMM<-depmix(list_HMM,input_HMM[,-1],nstates=nstates,family=list_families)
  
  HMMfit<-fit(HMM, verbose = FALSE) #fit our model to the data set
  
  idx_trans<-which(names(getpars(HMMfit))=="")
  
  transition_matrix<-matrix(getpars(HMMfit)[idx_trans],byrow=T,nrow=nstates)
  rownames(transition_matrix)<-paste("from",1:nstates,sep="")
  colnames(transition_matrix)<-paste("to",1:nstates,sep="")

  gp <- graph.adjacency(transition_matrix, mode = "directed", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  
  print("Step 6.1:  Plot output HMM")
  
  pdf(paste("HMM_graph",nstates,"pdf",sep="."),width=12)
  plot(gp,edge.label=round(E(gp)$weight, 3))
  dev.off()
  
  pdf(paste("HMM_graph_minimum_spanning_tree",nstates,"pdf",sep="."),width=12)
  plot(dp_mst,edge.label=round(E(dp_mst)$weight, 3))
  dev.off()
  
  print("Step 6.2:  Compute Posterior Probability")
  
  HMMpost<-posterior(HMMfit)
  
  output_HMM<-data.frame(samples=rownames(input_HMM),HMM_states=as.character(HMMpost[,1]),pseudospace=pseudospace_input_sort$pseudospace)
  
  return(output_HMM)
  
}
