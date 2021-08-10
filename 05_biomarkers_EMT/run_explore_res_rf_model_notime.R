library(data.table)
library(gridExtra)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

#Parameters
output_dir<-output_dir
tissue_correction = TRUE
string_output = "cosmic_arms_focal"
#

setwd(input_dir)

input_ml<-fread("input_for_ml_hmm_states_3_mock_as_gd.txt",data.table=F)
input_ml[is.na(input_ml)]<-0

idx_arm_gd_aneuploidy<-grep(colnames(input_ml),pattern=paste(c("arm","aneuploidy","genome_doubling"),collapse="|"),value=F)

input_ml_aga<-input_ml[,idx_arm_gd_aneuploidy]

#
# Select only the columns with dndscv and focal events
#
idx_columns<-grep(colnames(input_ml),pattern=paste(c("dndscv","focal"),collapse="|"),value=T)

#
# Load COSMIC genes for the analysis
#

cosmic_table<-read.csv(file="Census_allMon Aug 17 13_43_26 2020.csv",stringsAsFactors=F)
cosmic_genes<-paste(cosmic_table[,1],collapse="|")

#cosmic_genes<-paste(cosmic_table[,1],collapse="|")

#idx_columns2<-grep(idx_columns,pattern=cosmic_genes,value=T)

#
#update 21/10/2020
#

feat_cosmic_cnv<-paste(cosmic_table[,1],"_focal",sep="")
feat_cosmic_dndscv<-paste(cosmic_table[,1],"_dndscv",sep="")

idx_columns2<-c(feat_cosmic_cnv,feat_cosmic_dndscv)

#
# Create the new input for the algorithm
#

input_ml<-cbind(input_ml[,c(1:6)],input_ml[,which(colnames(input_ml)%in%idx_columns2)],input_ml_aga)

#
# Compare rf with random forest
#

setwd(output_dir)

list_rdata_rf<-c("RF_HMM_nstates3_mock.mes.vs.epi.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData","RF_HMM_nstates3_mock.mes.vs.mix.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData","RF_HMM_nstates3_mock.mix.vs.epi.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData")

list_class1<-c("mes","mes","mix")
list_class2<-c("epi","mix","epi")

list_rf_res<-vector(mode="list",2)

for(i in 1:length(list_rdata_rf)){

	print(i)

	class1<-list_class1[i]
	class2<-list_class2[i]

  input_ml2<-input_ml[input_ml$biological_states %in% c(class1,class2),c(5,6:ncol(input_ml))]

	#
	# Get top  genes from rf
	#
	setwd(input_dir)
	
	load(list_rdata_rf[i])

	rf_features_top<-features_selected_rf[,1]

	rf_columns<-which(colnames(input_ml2) %in% c("biological_states","tumors",as.character(rf_features_top)))

  input_ml2$biological_states<-as.factor(gsub(gsub(input_ml2[,1],pattern=class1,replacement=1),pattern=class2,replacement=0))
	
	idx_for_training<-sample(round(nrow(input_ml2)*80/100))

	training<-input_ml2[idx_for_training,rf_columns]
	test<-input_ml2[-idx_for_training,rf_columns]

	tissue_common<-intersect(training[,2],test[,2])

  training<-training[training[,2]%in%tissue_common,]
  test<-test[test[,2]%in%tissue_common,]

	library(ROCR)
	library(caret)

	#
	# Model 
	#

	glmnet <- glmnetUtils::cv.glmnet(biological_states~.,training,family="binomial",nfolds=5,alpha=1,type.measure = "auc")
	
	glmnet_prediction<-predict(glmnet,newdata=test,type="response")
	glmnet_prediction2<-predict(glmnet,newdata=test,type="class")

	per<-data.frame(confusionMatrix(as.factor(glmnet_prediction2),as.factor(test$biological_states),positive="1")$byClass)
	per<-cbind(rownames(per),per)
	colnames(per)<-c("perf","values")

	pred <- prediction(glmnet_prediction, test$biological_states)
	auc<-round(unlist(performance(pred,"auc")@y.values),3)
	perf <- performance(pred,"tpr","fpr")
	
	setwd(output_dir)
	
	pdf(paste("ROCR_selected_RF_cvglmnet_",class1,"_vs_",class2,".pdf",sep=""))
	plot(perf,main="CV.GLMNET") 
  lines(c(0,1),c(0,1),col = "gray", lty = 4)
	text(x=0.7,y=0.1,paste("AUC-:",auc))
	grid.table(per,rows=NULL)
	dev.off()

	ctrl <- trainControl(method="repeatedcv",repeats=10,number=5, summaryFunction=twoClassSummary,p=0.80,classProbs=T, savePredictions = T)

	train_mtx<-training
  test_mtx<-test

	train_mtx[,1]<-ifelse(train_mtx[,1]==1,"yes","no")
  levels(train_mtx[,1])<-c("yes","no")

  test_mtx[,1]<-ifelse(test_mtx[,1]==1,"yes","no")
  levels(test_mtx[,1])<-c("yes","no")

  # Random Forest
  rfmodel <- train(biological_states ~ .,data=train_mtx,method="ranger",trControl=ctrl)
  pred <- predict(rfmodel, newdata=test_mtx[,-1], type="prob")
  rf_pred_df<-data.frame(pred, test_mtx$biological_states,"ranger")
  colnames(rf_pred_df)<-c("no","yes","obs","Group")
  
  # GBM
  gbmmodel <- train(biological_states ~ .,data=train_mtx,method="gbm",trControl=ctrl)
  pred <- predict(gbmmodel, newdata=test_mtx[,-1], type="prob")
  gbm_pred_df<-data.frame(pred, test_mtx$biological_states,"gbm")
  colnames(gbm_pred_df)<-c("no","yes","obs","Group")
  
  
  # glm
  glmmodel <- train(biological_states ~ .,data=train_mtx,method="glm",trControl=ctrl)
  pred <- predict(glmmodel, newdata=test_mtx[,-1], type="prob")
  glm_pred_df<-data.frame(pred, test_mtx$biological_states,"glm")
  colnames(glm_pred_df)<-c("no","yes","obs","Group")
  
  
  # nb
  nbmodel <- train(biological_states ~ .,data=train_mtx,method="nb",trControl=ctrl)
  pred <- predict(nbmodel, newdata=test_mtx[,-1], type="prob")
  nb_pred_df<-data.frame(pred, test_mtx$biological_states,"nb")
  colnames(nb_pred_df)<-c("no","yes","obs","Group")
  
  
  ## run MLeval
  library(MLeval)
  
  input_evalm<-rbind(rf_pred_df,gbm_pred_df,glm_pred_df,nb_pred_df)
  
  pdf(paste("EVALM_5cv_selected_RF_",class1,"_vs_",class2,".pdf",sep=""))
  res <- evalm(input_evalm)
  dev.off()

}
