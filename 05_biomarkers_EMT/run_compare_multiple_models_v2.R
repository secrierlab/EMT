library(data.table)
library(gridExtra)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

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
setwd(input_dir)
cosmic_table<-read.csv(file="Census_allMon Aug 17 13_43_26 2020.csv",stringsAsFactors=F)
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
# Compare lasso with random forest
#

setwd(output_dir)

list_rdata_lasso<-c("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData","HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData",
"HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

list_class1<-c("mes","mes","mix")
list_class2<-c("epi","mix","epi")

list_lasso_res<-vector(mode="list",2)

for(i in 1:length(list_rdata_lasso)){

	print(i)

	class1<-list_class1[i]
	class2<-list_class2[i]

  input_ml2<-input_ml[input_ml$biological_states %in% c(class1,class2),c(5,6:ncol(input_ml))]

	#
	# Get top 50 genes from lasso
	#

	setwd("/data/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")
	
	load(list_rdata_lasso[i])

	df_coef2<-do.call(rbind,res_lasso)
	
	df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

  features<-data.frame(table(df_coef2[,1]))

  lasso_features_top50<-features[order(features[,2],decreasing=T),1][1:50]
	
	#
	# Get top 80% genes from lasso
	#

	# get the most frequent (80%) genes from lasso
  df_coef2<-do.call(rbind,res_lasso)

  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

  features<-data.frame(table(df_coef2[,1]))

  features_freq_df<-features[order(features[,2],decreasing=T),]

  lasso_features_80<-features_freq_df[features_freq_df[,2]>=800,1]

	

	lasso_50_columns<-which(colnames(input_ml2) %in% c("biological_states","tumors",as.character(lasso_features_top50)))
  lasso_80_columns<-which(colnames(input_ml2) %in% c("biological_states","tumors",as.character(lasso_features_80)))

        				
  input_ml2$biological_states<-as.factor(gsub(gsub(input_ml2[,1],pattern=class1,replacement=1),pattern=class2,replacement=0))
	
	idx_for_training<-sample(round(nrow(input_ml2)*80/100))

	training_50<-input_ml2[idx_for_training,lasso_50_columns]
	test_50<-input_ml2[-idx_for_training,lasso_50_columns]

	tissue_common_50<-intersect(training_50[,2],test_50[,2])

  training_50<-training_50[training_50[,2]%in%tissue_common_50,]
  test_50<-test_50[test_50[,2]%in%tissue_common_50,]

  training_80<-input_ml2[idx_for_training,lasso_80_columns]
  test_80<-input_ml2[-idx_for_training,lasso_80_columns]

  tissue_common_80<-intersect(training_80[,2],test_80[,2])

  training_80<-training_80[training_80[,2]%in%tissue_common_80,]
  test_80<-test_80[test_80[,2]%in%tissue_common_80,]

	library(ROCR)
	library(caret)

	#
	# Model 50
	#

	glmnet50 <- glmnetUtils::cv.glmnet(biological_states~.,training_50,family="binomial",nfolds=5,alpha=1,type.measure = "auc")
	
	glmnet50_prediction<-predict(glmnet50,newdata=test_50,type="response")
	glmnet50_prediction2<-predict(glmnet50,newdata=test_50,type="class")

	per<-data.frame(confusionMatrix(as.factor(glmnet50_prediction2),as.factor(test_50$biological_states),positive="1")$byClass)
	per50<-cbind(rownames(per),per)
	colnames(per50)<-c("perf","values")

	pred <- prediction(glmnet50_prediction, test_50$biological_states)
	auc50<-round(unlist(performance(pred,"auc")@y.values),3)

	perf50 <- performance(pred,"tpr","fpr")
	

	#
	# Model 80 
	#

	glmnet80 <- glmnetUtils::cv.glmnet(biological_states~.,training_80,family="binomial",nfolds=5,alpha=1,type.measure = "auc")
		
  glmnet80_prediction<-predict(glmnet80,newdata=test_80,type="response")
  glmnet80_prediction2<-predict(glmnet80,newdata=test_80,type="class")

  per<-data.frame(confusionMatrix(as.factor(glmnet80_prediction2),as.factor(test_80$biological_states),positive="1")$byClass)
  per80<-cbind(rownames(per),per)
  colnames(per80)<-c("perf","values")


  pred <- prediction(glmnet80_prediction, test_80$biological_states)
  auc80<-round(unlist(performance(pred,"auc")@y.values),3)

  perf80 <- performance(pred,"tpr","fpr")

	perall<-cbind(per50[,-1],per80[,-1])
	rownames(perall)<-per50[,1]
	colnames(perall)<-c("50","80")
  
	setwd(output_dir)
  
	pdf(paste("ROCR_lasso_top50_80perc_",class1,"_vs_",class2,".pdf",sep=""))
	plot(perf50,main="CV.GLMNET") 
	lines(perf80@x.values[[1]], perf80@y.values[[1]], col = 2)
        lines(c(0,1),c(0,1),col = "gray", lty = 4)
	text(x=0.7,y=0.1,paste("AUC-50:",auc50))
	text(x=0.7,y=0,paste("AUC-80:",auc80),col=2)
	grid.table(perall)
	dev.off()

	ctrl <- trainControl(method="repeatedcv",repeats=10,number=5, summaryFunction=twoClassSummary,p=0.80,classProbs=T, savePredictions = T)
	
	temp_input50<-input_ml2[,lasso_50_columns]
	levels(temp_input50[,1]) <- c("no", "yes")

	list_training<-c("training_50","training_80")
	list_testing<-c("test_50","test_80")
	list_output<-c("EVALM_5cv_lasso_top50","EVALM_5cv_lasso_perc80")	

	for(ml in 1:length(list_training)){

	train_mtx<-get(list_training[ml])
	test_mtx<-get(list_testing[ml])
	string_output<-list_output[ml]

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

	pdf(paste(string_output,class1,"_vs_",class2,".pdf",sep=""))
	res <- evalm(input_evalm)	
	dev.off()
	
	}
	
	setwd(input_dir)
	
}
