library(data.table)
library(caret)
library(biglasso)
library(glmnetUtils)
library(doParallel)


folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

#Parameters
output_dir<-output_dir
tissue_correction = TRUE
string_output = "cosmic_arms_focal"


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
	
	trainIndex <- createDataPartition(input_ml2$biological_states, p = .8, 
                                  list = FALSE, 
                                  times = 1000)

	res_lasso<-vector(mode="list",ncol(trainIndex))

	overall_performance<-data.frame()
	overall_bygroup<-data.frame()

	for(ti in 1:ncol(trainIndex)){
	
	print(paste("Iteration number:",ti))

	vector_partition<-trainIndex[,ti]
	
	training<-input_ml2[vector_partition,]
	test<-input_ml2[-vector_partition,]

	#take only the samples in which there is the same tissue between training and test - prevent errors during the prediction
	tissue_common<-intersect(training[,2],test[,2])

	training<-training[training[,2]%in%tissue_common,]
	test<-test[test[,2]%in%tissue_common,]

	if(tissue_correction == TRUE){
    	
	cvfit <- glmnetUtils::cv.glmnet(biological_states~.,training,family="binomial",nfolds=5,alpha=1,parallel=T)

  variables_ok<- data.frame(coef(cvfit)[which(coef(cvfit) != 0),])

  res_lasso[[ti]]<-data.frame(coef=rownames(variables_ok),values=variables_ok)
  predictions<-as.numeric(predict(cvfit,newdata=test[,-1],type="class"))

  overall_per<-confusionMatrix(as.factor(predictions),as.factor(test$biological_states),positive="1")$overal
  byclass_per<-confusionMatrix(as.factor(predictions),as.factor(test$biological_states),positive="1")$byClass

  overall_performance<-rbind(overall_performance,data.frame(ti,overall_per))
  overall_bygroup<-rbind(overall_bygroup,data.frame(ti,byclass_per))
	
	}else{

	biological_states<-training[,"biological_states"]

  X.bm <- as.big.matrix(as.matrix(training[,-(colnames(training)%in%c("biological_states","tumors"))]))

	cvfit <-cv.biglasso(X.bm, biological_states, nfolds = 5, ncores = 10,penalty="lasso",family="binomial",alpha=1)
	
  variables_ok<- data.frame(coef(cvfit)[which(coef(cvfit) != 0),])
  
  res_lasso[[ti]]<-data.frame(coef=rownames(variables_ok),values=variables_ok)
  predictions<-as.numeric(predict(cvfit,X=as.big.matrix(test[,-1]),type="class"))

  overall_per<-confusionMatrix(as.factor(predictions),as.factor(test$biological_states),positive="1")$overal
  byclass_per<-confusionMatrix(as.factor(predictions),as.factor(test$biological_states),positive="1")$byClass

  overall_performance<-rbind(overall_performance,data.frame(ti,overall_per))
  overall_bygroup<-rbind(overall_bygroup,data.frame(ti,byclass_per))

	}

	}

	setwd(output_dir)

	overall_bygroup_res<-do.call(cbind,split(overall_bygroup,overall_bygroup$ti))
	overall_bygroup_res2<-overall_bygroup_res[,grep(colnames(overall_bygroup_res),pattern="ti",invert=T)]

	pdf(paste("HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"byGroup.pdf",sep="."))
	boxplot(t(overall_bygroup_res2),las=2)
	dev.off()

  overall_performance_res<-do.call(cbind,split(overall_performance,overall_performance$ti))
  overall_perfomance_res2<-overall_performance_res[,grep(colnames(overall_performance_res),pattern="ti",invert=T)]

  pdf(paste("HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"performance.pdf",sep="."))
  boxplot(t(overall_perfomance_res2),las=2)
  dev.off()

	df_coef2<-do.call(rbind,res_lasso)
	df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

	features<-data.frame(table(df_coef2[,1]))

	features_to_select<-features[order(features[,2],decreasing=T),1][1:50]

	input_for_boxplot<-df_coef2[which(df_coef2[,1]%in%features_to_select),]
	colnames(input_for_boxplot)[2]<-"values"

	pdf(paste("HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"top_markers_50.pdf",sep="."))
	p<-ggplot(input_for_boxplot,aes(x=reorder(coef,values,.desc=TRUE, FUN = median), y=values))+geom_boxplot()+coord_flip()+geom_vline(xintercept = 0, linetype="dotted",color = "red", size=1)+theme(text = element_text(size=8))+geom_vline(xintercept=0, linetype="dashed", color = "red")
	print(p)
	dev.off()
	
	pdf(paste("HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),"top_markers_50_v2.pdf",sep="."))
        p<-ggplot(input_for_boxplot,aes(x=reorder(coef,values,.desc=TRUE, FUN = median), y=values))+geom_boxplot()+coord_flip(ylim=c(-5,5))+geom_vline(xintercept = 0, linetype="dotted",color = "red", size=1)+theme(text = element_text(size=8))+geom_vline(xintercept=0, linetype="dashed", color = "red")
	print(p)
	dev.off()

	save(list=c("res_lasso","overall_performance","overall_bygroup"),file=paste("HMM_nstates3_mock",class1,"vs",class2,"tissue",tissue_correction,ncol(trainIndex),string_output,"RData",sep="."))


}


