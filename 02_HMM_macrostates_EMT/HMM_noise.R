#
# Add noise HMM
#
library(depmixS4)
library(glmnet)
library(caret)
library(reshape2)
library(mltest)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

load("TCGA_matrix_gene_expression_signals_ALLGENES_29_01_2020.RData")

setwd(input_dir)
pseudospace_input<-read.delim(file="proj_pseudospace_mcf_tgfb_mock_Mclust3_with_pemt.txt")

pseudospace_input_sort<-pseudospace_input[order(pseudospace_input$mock,decreasing=F),]

setwd(input_dir)
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
markers_genes<-markers_genes_read[,1]

TCGA_GEXP_ALL_rid<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,]

TCGA_GEXP_ALL_sort<-log(t(TCGA_GEXP_ALL_rid[,match(pseudospace_input_sort$patients2,colnames(TCGA_GEXP_ALL_rid))])+1,2)
colnames(TCGA_GEXP_ALL_sort)<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,1]

#
# Original HMM
#

input_file<-c("HMM_results_nstates_3.txt")
knn_df_tcga<-read.delim(file=input_file)
knn_df_tcga<-knn_df_tcga[,which(colnames(knn_df_tcga)%in%c("samples","biological_states"))]

#
# Start the procedure
#

nstates <- 3
noise_level<-c("default",10,60,125,250,500,1000,2000,3000,4000,5000,5500)

all_res_noise<-data.frame()
all_accuracy_noise<-data.frame()
all_delta_noise<-data.frame()
all_tissues_noise<-data.frame()

for(nl in noise_level){
  
all_res<-data.frame()
all_accuracy<-data.frame()
all_delta<-data.frame()
all_tissue<-data.frame()

for(i in 1:100){
    
    if(nl=="default"){
    
    TCGA_GEXP_ALL_noise<-apply(TCGA_GEXP_ALL_sort,2,jitter)

    }else{
      
    nl<-as.numeric(nl)

    TCGA_GEXP_ALL_noise<-apply(TCGA_GEXP_ALL_sort,2,jitter,nl)  
    
    }
  
    cvfit <-cv.glmnet(x = TCGA_GEXP_ALL_noise, y = pseudospace_input_sort$mock, alpha = 1)
    
    dfcoef<-as.matrix(coef(cvfit,s= "lambda.min"))
    genes_to_use<- names(dfcoef[dfcoef[,1]!=0,])[-1]
    
    list_HMM<- lapply(as.list(paste(colnames(TCGA_GEXP_ALL_noise[,genes_to_use]),"~1")),as.formula)
    
    list_families<-vector(mode="list",length(genes_to_use))

    for(lf in 1:length(genes_to_use)){
      
      list_families[[lf]]<-gaussian()
    }


    input_HMM<-data.frame(time=pseudospace_input_sort$mock,TCGA_GEXP_ALL_noise[,genes_to_use])
    
    HMM<-depmix(list_HMM,input_HMM[,-1],nstates=nstates,family=list_families)
    
    HMMfit<-fit(HMM, verbose = FALSE) #fit our model to the data set

    #use also summary summary(HMMfit, which = "transition") to get the transitions

    idx_trans<-which(names(getpars(HMMfit))=="")
    
    #https://stackoverflow.com/questions/50327876/obtain-transition-matrix-in-r-depmixs4-package
    
    transition_matrix<-matrix(getpars(HMMfit)[idx_trans],byrow=T,nrow=nstates)
    rownames(transition_matrix)<-paste("from",1:nstates,sep="")
    colnames(transition_matrix)<-paste("to",1:nstates,sep="")
    
    HMMpost<-posterior(HMMfit)
    HMMpost_withsamples<-data.frame(patientsID=rownames(TCGA_GEXP_ALL_noise),states=HMMpost$state)
 
    # Compare with the original dataset
    
    HMM_res_comp<-merge(knn_df_tcga,HMMpost_withsamples,by.x="samples",by.y="patientsID")
    
    res_table<-table(HMM_res_comp$biological_states,HMM_res_comp$states)
    
    samples_for_states<-data.frame(noise=as.character(nl),t(apply(res_table,1,max)))
    all_res<-rbind(all_res,samples_for_states)
    
    res_table_for_stat<-res_table
    
    # take the indices in which there is the maximum number of samples in each column
    # this is to determine if a state e.g. 1/2/3 is more likelihood to be epi/mix/mes
    idx_to_assign_states<-apply(res_table_for_stat,2,which.max)
    
    if(length(unique(idx_to_assign_states))==3){
    
    # rename the column names of the contingency table
    # use the rownames of the contingency table sorted by the indices get before
    colnames(res_table_for_stat)<-rownames(res_table_for_stat)[idx_to_assign_states]

    stats_confmatrix<-confusionMatrix(res_table_for_stat[c("epi","mix","mes"),c("epi","mix","mes")])
    accuracy<-data.frame(noise=as.character(nl),accuracy=stats_confmatrix$overall[1])
    
    print(res_table_for_stat)
    
    original_states<-table(HMM_res_comp$biological_states)
    df_prediction_states<-data.frame(res_table_for_stat)
    
    epi_count<-df_prediction_states[df_prediction_states[,1]=="epi"&df_prediction_states[,2]=="epi",3]
    mix_count<-df_prediction_states[df_prediction_states[,1]=="mix"&df_prediction_states[,2]=="mix",3]
    mes_count<-df_prediction_states[df_prediction_states[,1]=="mes"&df_prediction_states[,2]=="mes",3]
    
    predicted_states<-c(epi_count,mix_count,mes_count)
    names(predicted_states)<-c("epi","mix","mes")
    
    # delta_states<-data.frame(iteration=i,noise=as.character(nl),data.frame(original_states[c("epi","mix","mes")]/predicted_states[c("epi","mix","mes")]))
    delta_states<-data.frame(iteration=i,noise=as.character(nl),data.frame(predicted_states[c("epi","mix","mes")]/original_states[c("epi","mix","mes")]))
    
    colnames(delta_states)<-c("iteration","noise","states","delta")
    
    #tissues levels statistics
    names(idx_to_assign_states)<-rownames(res_table_for_stat)[idx_to_assign_states]
    
    HMM_res_comp$states<-gsub(HMM_res_comp$states,pattern=idx_to_assign_states[1],replacement=names(idx_to_assign_states)[1])
    HMM_res_comp$states<-gsub(HMM_res_comp$states,pattern=idx_to_assign_states[2],replacement=names(idx_to_assign_states)[2])
    HMM_res_comp$states<-gsub(HMM_res_comp$states,pattern=idx_to_assign_states[3],replacement=names(idx_to_assign_states)[3])
    HMM_res_comp$tissue<-sapply(strsplit(as.character(HMM_res_comp[,1]),split="\\."),"[[",1)
      
    df_tissue<-data.frame(iteration=i,noise=as.character(nl),table(HMM_res_comp$states,HMM_res_comp$tissue))
    colnames(df_tissue)<-c("iteration","noise","states","tissue","freq_state_tissue")
    
    }else{
      
    accuracy<-data.frame(noise=as.character(nl),accuracy=NA) 
    delta_states<-data.frame(iteration=i,noise=as.character(nl),states=NA,delta=NA)
    df_tissue<-data.frame(iteration=i,noise=as.character(nl),states=NA,tissue=NA,freq_state_tissue=NA)
    
    }
    
    all_accuracy<-rbind(all_accuracy,accuracy)
    all_delta<-rbind(all_delta,delta_states)
    all_tissue<-rbind(all_tissue,df_tissue)
    
}

    all_res_noise<-rbind(all_res_noise,all_res)
    all_accuracy_noise<-rbind(all_accuracy_noise,all_accuracy)
    all_delta_noise<-rbind(all_delta_noise,all_delta)
    all_tissues_noise<-rbind(all_tissues_noise,all_tissue)
}

setwd(output_dir)

save(list=c("all_res_noise","all_accuracy_noise","all_delta_noise","all_tissues_noise"),file="HMM_noise.RData")

all_res_noise_melt<-data.frame(iterations=rep(1:100,12),melt(all_res_noise))

pdf("HMM_noise_number_samples_for_states.pdf",width=15,height=8)
all_res_noise_melt$variable<-factor(all_res_noise_melt$variable,levels=c("epi","mix","mes"))
p<-ggplot(all_res_noise_melt,aes(x=noise,y=value,colour=noise)) +
  geom_boxplot()+facet_wrap(~variable)+ labs(y="Number of samples in a state")
print(p)
p2<-ggplot(all_res_noise_melt,aes(x=noise,y=value,colour=noise)) +
  geom_boxplot()+geom_jitter(color="black", size=0.4, alpha=0.9)+facet_wrap(~variable)+ labs(y="Number of samples in a state")
print(p2)
dev.off()

all_accuracy_boxplot<-data.frame(iterations=rep(1:100,12),all_accuracy_noise)

pdf("HMM_noise_accuracy.pdf")
p2<-ggplot(all_accuracy_noise,aes(x=noise,y=accuracy,colour=noise)) +
  geom_boxplot()+ labs(y="Accuracy")
print(p2)
dev.off()

pdf("HMM_noise_delta.pdf",width=15,height=8)
all_delta_noise2<-all_delta_noise[!is.na(all_delta_noise$states),]
p2<-ggplot(all_delta_noise2,aes(x=noise,y=delta,colour=noise)) +
  geom_boxplot()+facet_wrap(~states)+labs(y="NumberSamplesSimulation/NumberSamplesOriginal")
print(p2)
dev.off()


pdf("HMM_noise_number_samples_for_states.TISSUE.pdf",width=15,height=8)

tissues_agg<-all_tissues_noise[,c(2,3,4,5)]
mean.df <- aggregate(freq_state_tissue ~ tissue+noise+states, tissues_agg, median)

mean.df$states<-factor(mean.df$states,levels=c("epi","mix","mes"))
p<-ggplot(mean.df,aes(x=noise,y=freq_state_tissue,colour=noise)) +
  geom_boxplot()+facet_wrap(~states)+ labs(y="Median Number of samples in a state by tissue")
print(p)

tissues_agg<-all_tissues_noise[,c(2,3,4,5)]
# mean.df <- aggregate(freq_state_tissue ~ tissue+noise+states, tissues_agg, sum)
mean.df <- aggregate(freq_state_tissue ~ tissue+noise+states, tissues_agg, mean)

mean.df$states<-factor(mean.df$states,levels=c("epi","mix","mes"))
p2<-ggplot(mean.df,aes(x=noise,y=freq_state_tissue,colour=noise)) +
  geom_boxplot()+facet_wrap(~states)+ labs(y="Mean Number of samples in a state by tissue")
print(p2)


dev.off()




