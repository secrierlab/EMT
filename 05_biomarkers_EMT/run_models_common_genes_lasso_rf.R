library(data.table)
library(gridExtra)

setwd("/data/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")
load("common_genomic_features_lasso_random_forest.RData")

#Parameters
output_dir<-"/data/pseudospace/ml_for_ppt/cosmic_focal_broad_variants"
tissue_correction = FALSE
string_output = "cosmic_arms_focal"
#

setwd("/data/pseudospace/ml_for_ppt")

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

setwd(output_dir)

list_class1<-c("mes","mes","mix")
list_class2<-c("epi","mix","epi")


for(i in 1:length(list_class1)){

	print(i)

	class1<-list_class1[i]
	class2<-list_class2[i]

	input_ml2<-input_ml[input_ml$biological_states %in% c(class1,class2),c(5,6:ncol(input_ml))]


#
# Upload the list of genes in common between lasso and random forest
#

	list_genes_to_use<-c("biological_states",list_common_genes_between_experiments[[paste(class1,"_vs_",class2,sep="")]])

	columns_to_consider<-which(colnames(input_ml2) %in% c("biological_states",as.character(list_genes_to_use)))

	input_ml2$biological_states<-as.factor(gsub(gsub(input_ml2[,1],pattern=class1,replacement=1),pattern=class2,replacement=0))

	idx_for_training<-sample(round(nrow(input_ml2)*70/100))

        training<-input_ml2[idx_for_training,columns_to_consider]
        test<-input_ml2[-idx_for_training,columns_to_consider]

        tissue_common<-intersect(training[,2],test[,2])

        training<-training[training[,2]%in%tissue_common,]
        test<-test[test[,2]%in%tissue_common,]

	library(glmnetUtils)
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


        pdf(paste("ROCR_common_genes_LASSO_with_RF_cvglmnet_",class1,"_vs_",class2,".pdf",sep=""))
        plot(perf,main="CV.GLMNET")
        lines(c(0,1),c(0,1),col = "gray", lty = 4)
        text(x=0.7,y=0.1,paste("AUC-:",auc))
        grid.table(per,rows=NULL)
        dev.off()

        ctrl <- trainControl(method="cv",number=5, sampling = "down",summaryFunction=twoClassSummary,p=0.70,classProbs=T, savePredictions = T)

	train_mtx<-training
	test_mtx<-test	
	
	train_mtx[,1]<-ifelse(train_mtx[,1]==1,"yes","no")
	test_mtx[,1]<-ifelse(test_mtx[,1]==1,"yes","no")

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
        nbmodel <- train(biological_states ~ .,data=train_mtx,method="glmnet",trControl=ctrl)
        pred <- predict(nbmodel, newdata=test_mtx[,-1], type="prob")
        nb_pred_df<-data.frame(pred, test_mtx$biological_states,"nb")
        colnames(nb_pred_df)<-c("no","yes","obs","Group")


        ## run MLeval
        library(MLeval)

        input_evalm<-rbind(rf_pred_df,gbm_pred_df,glm_pred_df,nb_pred_df)

        ## run MLeval
        library(MLeval)

        pdf(paste("EVALM_5cv_common_genes_LASSO_with_RF_",class1,"_vs_",class2,".pdf",sep=""))        
	res <- evalm(input_evalm)
	dev.off()

}
