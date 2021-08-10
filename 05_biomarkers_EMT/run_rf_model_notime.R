library(data.table)
library(caret)
library(ranger)
library(fastshap)

source("autoplot2.R")

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

#Parameters
output_dir<-output_dir
tissue_correction = FALSE
string_output = "cosmic_arms_focal_train_0.80"
top_genes = 800
#

setwd(input_dir)

input_ml<-fread("input_for_ml_hmm_states_3_mock_as_gd.txt",data.table=F)
input_ml[is.na(input_ml)]<-0

#
# Select a data.frame with only arm, aneuploidy and genome doubling
#

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

group1<-c("mes","mix","mes")
group2<-c("epi","epi","mix")

for(i in 1:length(group1)){

	class1<-group1[i]
	class2<-group2[i]
	
	print(paste(class1,class2))

	input_ml2<-input_ml[input_ml$biological_states %in% c(class1,class2),c(5,6:ncol(input_ml))]
  input_ml2$biological_states<-as.factor(gsub(gsub(input_ml2[,1],pattern=class1,replacement=1),pattern=class2,replacement=0))

	#
	# Run Ranger feature selection + model 
	#

	trainIndex <- createDataPartition(input_ml2$biological_states, p = .8, 
                                  list = FALSE, 
                                  times = 1000)

	#
	# Define objects for random forest results
	#

	res_rf<-vector(mode="list",ncol(trainIndex))
	overall_performance_rf<-data.frame()
	overall_bygroup_rf<-data.frame()

	for(ti in 1:ncol(trainIndex)){
	
	print(paste("Iteration number:",ti))

	vector_partition<-trainIndex[,ti]
	
	training<-input_ml2[vector_partition,]
	test<-input_ml2[-vector_partition,]

	#take only the samples in which there is the same tissue between training and test - prevent errors during the prediction
	tissue_common<-intersect(training[,2],test[,2])

	training<-training[training[,2]%in%tissue_common,]
	test<-test[test[,2]%in%tissue_common,]

	#
	# RF: select most relevant variables
	#
	
	model_rf <- ranger(biological_states ~ ., data = training,importance="impurity_corrected",classification=T)
	rf_pvalues<-importance_pvalues(model_rf, method = "janitza")
	variables_ok<-rf_pvalues[rf_pvalues[,2]<=0.00001,]
  res_rf[[ti]]<-data.frame(variable=rownames(variables_ok),variables_ok)

	#
	# RF: run model with best variables
	#

	training_sv<-training[,which(colnames(training)%in%c("biological_states",rownames(variables_ok)))]

	model_rf_sv<-ranger(biological_states ~ ., data = training_sv,importance='impurity',num.threads=60,classification=T)
  predictions_rf<-predict(model_rf_sv,data= test[,-1])    
	confMatrix_rf<-table(predictions_rf$prediction,test[,1])[2:1,2:1]

	#source https://towardsdatascience.com/understanding-confusion-matrix-a9ad42dcfd62
	#source https://machinelearningmastery.com/confusion-matrix-machine-learning/

	#expected <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0))
	#predicted <- factor(c(1, 0, 0, 1, 0, 0, 1, 1, 1, 0))
	#results <- confusionMatrix(data=predicted, reference=expected,positive="1")
	#confpart <- table(predicted,expected)[2:1,2:1]

  overall_per_rf<-confusionMatrix(confMatrix_rf,positive="1")$overal
  byclass_per_rf<-confusionMatrix(confMatrix_rf,positive="1")$byClass

  overall_performance_rf<-rbind(overall_performance_rf,data.frame(ti,overall_per_rf))
  overall_bygroup_rf<-rbind(overall_bygroup_rf,data.frame(ti,byclass_per_rf))

	}

	setwd(output_dir)

	overall_bygroup_res<-do.call(cbind,split(overall_bygroup_rf,overall_bygroup_rf$ti))
	overall_bygroup_res2<-overall_bygroup_res[,grep(colnames(overall_bygroup_res),pattern="ti",invert=T)]

	pdf(paste("RF_HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"byGroup.pdf",sep="."))
	boxplot(t(overall_bygroup_res2),las=2)
	dev.off()

  overall_performance_res<-do.call(cbind,split(overall_performance_rf,overall_performance_rf$ti))
  overall_perfomance_res2<-overall_performance_res[,grep(colnames(overall_performance_res),pattern="ti",invert=T)]

  pdf(paste("RF_HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"topgenes",top_genes,"performance.pdf",sep="."))
  boxplot(t(overall_perfomance_res2),las=2)
  dev.off()

	df_coef2<-do.call(rbind,res_rf)
	df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

	features<-data.frame(table(df_coef2[,1]))

	features_selected_rf<-features[features[,2]>=top_genes,]
	colnames(features_selected_rf)<-c("variables","freq")

	pdf(paste("RF_HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"topgenes",top_genes,"top_markers_feature_selected.pdf",sep="."))
	p<-ggplot(data=features_selected_rf, aes(x=variables, y=freq))+geom_bar(stat="identity")+theme(axis.text=element_text(size=5))
	# Horizontal bar plot
	p2<-p + coord_flip()
	print(p2)
	dev.off()

	#
	# Shap analysis on the predicted genes
	#
	
  train_index <- sample(1:nrow(input_ml2), 0.8 * nrow(input_ml2))
  
  train_data<-input_ml2[train_index,which(colnames(input_ml2)%in%c("biological_states",as.character(features_selected_rf$variables)))]
  test_data<-input_ml2[-train_index,which(colnames(input_ml2)%in%c("biological_states",as.character(features_selected_rf$variables)))]
        
  model_rf <- ranger(biological_states ~ . , data = train_data,probability = TRUE)
             
  pfun <- function(object, newdata) { 
  predict(object, data = newdata)$predictions[,2]
  }

  exp_rf<- fastshap::explain(model_rf, X = test_data[,-1], pred_wrapper = pfun, nsim = 10,.parallel = TRUE)

	pdf(paste("RF_HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"topgenes",top_genes,"top_markers_feature_selected.SHAP.pdf",sep="."))
	p3<-autoplot2(exp_rf,type="contribution")
	print(p3)
  dev.off()
	
	save(list=c("res_rf","overall_performance_rf","overall_bygroup_rf","features_selected_rf"),
  file=paste("RF_HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"topgenes",top_genes,string_output,"RData",sep="."))



}


