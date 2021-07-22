#
# Integrate the information between the markers, EMT, Fibroblasts, hypoxia 
#
library(ggpubr)
library(data.table)
library(gridExtra)
library(viridis)
library(plyr)

# The list of fibroblasts are those one 
load("/home/data/pseudospace/pathway_characterization/broad_genomic_features/aneuploidy_stem_hypo_TME_Fibroblasts_TCGA.RData")
input_gf_tme_fb2<-input_gf_tme_fb
input_gf_tme_fb2$states<-gsub(input_gf_tme_fb2$states,pattern="1",replacement="pEMT")
input_gf_tme_fb2$states<-gsub(input_gf_tme_fb2$states,pattern="2",replacement="epi")
input_gf_tme_fb2$states<-gsub(input_gf_tme_fb2$states,pattern="3",replacement="mes")

input_ml<-fread("/home/data/pseudospace/ml_for_ppt/input_for_ml_hmm_states_3_mock_as_gd.txt",data.table=F)
colnames(input_ml)<-gsub(colnames(input_ml),pattern="^X",replacement="")
input_ml$patients<-unlist(lapply(strsplit(as.character(input_ml$patients),split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
input_ml2<-input_ml

input_ml<-input_ml[input_ml$patients%in%input_gf_tme_fb2$samplesID,]

# The list of fibroblasts are those one that have been obtained from the glm analysis
list_fibroblasts_to_consider_mes_vs_epi<-c("Fibroblasts","COL11A1_FS","Esophageal_CAF_markers_carcinogenesis","BLADDER_CAF_markers")

setwd("/home/data/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

getFeaturesSS<-function(res_lasso,ss=800){
  
  df_coef2<-do.call(rbind,res_lasso)
  
  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]
  
  features<-data.frame(table(df_coef2[,1]))
  
  features_freq_df<-features[order(features[,2],decreasing=T),]
  
  features_to_select<-features_freq_df[features_freq_df[,2]>=ss,1]
  
  return(features_to_select)
}


getCoefSS<-function(res_lasso,features_to_select){
  
  df_coef2<-do.call(rbind,res_lasso)
  df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]
  
  features<-data.frame(table(df_coef2[,1]))
  
  input_for_boxplot<-df_coef2[which(df_coef2[,1]%in%features_to_select),]
  colnames(input_for_boxplot)[2]<-"values"
  
  return(input_for_boxplot)
  
}

setwd("/home/data/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
mut_genes_mes_vs_epi<-sapply(strsplit(grep(markers_mes_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)
cnv_markers_mes_vs_epi<-gsub(grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
cnv_markers_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T)
coef_mes_vs_epi<-getCoefSS(res_lasso,c(mut_genes_mes_vs_epi_raw,cnv_markers_mes_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_mes_vs_epi2<-merge(coef_mes_vs_epi,mu,by="coef")
coef_mes_vs_epi2[,1]<-gsub(coef_mes_vs_epi2[,1],pattern="^X",replacement="")


load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
mut_genes_pemt_vs_epi<-sapply(strsplit(grep(markers_pemt_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)
cnv_markers_pemt_vs_epi<-gsub(grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
cnv_markers_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T)
coef_pemt_vs_epi<-getCoefSS(res_lasso,c(mut_genes_pemt_vs_epi_raw,cnv_markers_pemt_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_pemt_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_pemt_vs_epi2<-merge(coef_pemt_vs_epi,mu,by="coef")
coef_coef_pemt_vs_epi2[,1]<-gsub(coef_coef_pemt_vs_epi2[,1],pattern="^X",replacement="")
#
raw_markers_specifics_mes_vs_epi<-setdiff(c(mut_genes_mes_vs_epi_raw,cnv_markers_mes_vs_epi),c(mut_genes_pemt_vs_epi_raw,cnv_markers_pemt_vs_epi))
raw_markers_specifics_pemt_vs_epi<-setdiff(y=c(mut_genes_mes_vs_epi_raw,cnv_markers_mes_vs_epi),x=c(mut_genes_pemt_vs_epi_raw,cnv_markers_pemt_vs_epi))
raw_markers_common<-intersect(c(mut_genes_mes_vs_epi_raw,cnv_markers_mes_vs_epi),c(mut_genes_pemt_vs_epi_raw,cnv_markers_pemt_vs_epi))

common_markers<-intersect(c(mut_genes_mes_vs_epi_raw,cnv_markers_mes_vs_epi),c(mut_genes_pemt_vs_epi_raw,cnv_markers_pemt_vs_epi))

multiOmicsCharts<-function(X,genomic_matrix,input_ml,cnv_thr=0.7,suffix_output,heigth=10){
  
                group1<-X[[4]][1]
                group2<-X[[4]][2]
  
                genomic_matrix_subset<-genomic_matrix[genomic_matrix$states%in%c(group1,group2),colnames(genomic_matrix)%in%c("samplesID",X[[1]])]
                
                #Discretize the mutations data 1 (yes) 0 (no)
                mutations<-paste(X[[2]],"_dndscv",sep="")
                cnv<-X[[3]]
                  
                input_ml_mut<-input_ml[input_ml$biological_states %in% c(group1,group2),colnames(input_ml)%in% c("patients","biological_states",mutations)]
                input_ml_cnv<-input_ml[input_ml$biological_states %in% c(group1,group2),colnames(input_ml)%in% c("patients","biological_states",cnv)]
                
                input_ml_cnv[,-c(1:2)]<-apply(input_ml_cnv[,-c(1:2)],2,FUN=function(X){ifelse(abs(X)>=cnv_thr,1,0)})
                
                
                #
                # Prepare mutation data 
                #
                
                feat_results_mut<-data.frame()
                
                for(markers in 3:ncol(input_ml_mut)){
                
                marker_name<-colnames(input_ml_mut)[markers]
                  
                all_samples<-input_ml_mut[,c(1,markers)]
                
                samples_with_mutations<-all_samples[all_samples[,2]%in%1,1]
                
                samples_without_mutations<-all_samples[all_samples[,2]%in%0,1]
                
                feat_results_temp<-data.frame()
                
                for(feature in 2:ncol(genomic_matrix_subset)){
                  
                print(feature)
                feature_name<-colnames(genomic_matrix_subset)[feature]
                
                samples_with_current_marker<-genomic_matrix_subset[genomic_matrix_subset$samplesID %in%samples_with_mutations,feature]
                samples_without_current_marker<-genomic_matrix_subset[genomic_matrix_subset$samplesID %in%samples_without_mutations,feature]
                

                if(length(samples_with_current_marker)>=1&length(samples_without_current_marker)>=1){
                  pval<-wilcox.test(samples_with_current_marker,samples_without_current_marker)$p.value
                }
                
                # else{
                # #   pval<-0.01
                # }
                
                if(length(samples_with_current_marker)==0&length(samples_without_current_marker)>=1){
                  
                  pval<-NA
                
                } 
                
                if(length(samples_with_current_marker)>=1&length(samples_without_current_marker)>=1){
                  or<-median(samples_with_current_marker)-median(samples_without_current_marker)
                }
                
                if(length(samples_without_current_marker)>1&length(samples_with_current_marker)==0){
                  or<-median(samples_without_current_marker)
                }
                
                if(length(samples_without_current_marker)==0 & length(samples_with_current_marker)>1){
                  or<-median(samples_with_current_marker)
                }
                
                feat_samples_with_mutations<-data.frame(marker_name,feature_name,pval=pval,or=or)
                colnames(feat_samples_with_mutations)<-c("genes","feature_value","pval","or")
                feat_results_temp<-rbind(feat_results_temp,feat_samples_with_mutations)

                }
                
                feat_results_mut<-rbind(feat_results_mut,feat_results_temp)
                
                }
                
  
                
                #
                # Prepare CNV data 
                #
                
                
                feat_results_cnv<-data.frame()
                
                for(markers in 3:ncol(input_ml_cnv)){
                  
                  marker_name<-colnames(input_ml_cnv)[markers]
                  
                  all_samples<-input_ml_cnv[,c(1,markers)]
                  samples_with_cnv<-all_samples[all_samples[,2]%in%1,1]
                  samples_without_cnv<-all_samples[all_samples[,2]%in%0,1]
                  
                  feat_results_temp<-data.frame()
                  
                  for(feature in 2:ncol(genomic_matrix_subset)){
                    
                    print(feature)
                    feature_name<-colnames(genomic_matrix_subset)[feature]
                    
                    samples_with_current_marker<-genomic_matrix_subset[genomic_matrix_subset$samplesID %in%samples_with_cnv,feature]
                    samples_without_current_marker<-genomic_matrix_subset[genomic_matrix_subset$samplesID %in%samples_without_cnv,feature]
                    # if(length(samples_with_current_marker)>1&length(samples_with_current_marker)>1){
                    # pval<-wilcox.test(samples_with_current_marker,samples_without_current_marker)$p.value
                    # }else{
                    #   pval<-0.01
                    # }
                    
                    if(length(samples_with_current_marker)>=1&length(samples_without_current_marker)>=1){
                      pval<-wilcox.test(samples_with_current_marker,samples_without_current_marker)$p.value
                    }
                    
                    # else{
                    # #   pval<-0.01
                    # }
                    
                    if(length(samples_with_current_marker)==0&length(samples_without_current_marker)>=1){
                      
                      pval<-NA
                      
                    } 
                    
                    if(length(samples_with_current_marker)>1&length(samples_without_current_marker)>1){
                    or<-median(samples_with_current_marker)-median(samples_without_current_marker)
                    }
                    
                    if(length(samples_without_current_marker)>1&length(samples_with_current_marker)==0){
                      or<-median(samples_without_current_marker)
                    }
                    
                    if(length(samples_without_current_marker)==0&length(samples_with_current_marker)>1){
                      or<-median(samples_with_current_marker)
                    }
                    
                    
                    feat_samples_with_cnv<-data.frame(marker_name,feature_name,pval=pval,or=or)
                    colnames(feat_samples_with_cnv)<-c("genes","feature_value","pval","or")
                    feat_results_temp<-rbind(feat_results_temp,feat_samples_with_cnv)
                    
                  }
                  
                  feat_results_cnv<-rbind(feat_results_cnv,feat_results_temp)
                  
                }
                
                res_tot<-rbind(feat_results_mut,feat_results_cnv)
                res_tot$status_sig<-as.factor(ifelse(res_tot$pval<=0.05,"sig","no_sig"))
                res_tot$pval<- -log(res_tot$pval,10)
                #
                # Balloon chart
                #

                # dotchart
                
                list_dotcharts<-vector(mode="list",length(unique(res_tot$feature_value)))
                list_features<-unique(res_tot$feature_value)
                
                for(lfi in list_features){
                  res_tot[res_tot$feature_value%in%lfi,"or"]<-res_tot[res_tot$feature_value%in%lfi,"or"]/max(res_tot[res_tot$feature_value%in%lfi,"or"])
                }
                # res_tot$or<-res_tot$or/max(res_tot$or)
                  
                # res_tot[res_tot$or>0.3,"or"]<- 0.3
                # res_tot[res_tot$or< -0.3,"or"]<- -0.3
                
                # res_tot$markers<-factor(res_tot$genes,levels=rev(order_genes))
                
                # hec_sp <- ggplot(res_tot, aes(x = feature_value, y = genes)) +
                #   geom_point(aes(size = pval, shape =status_sig,colour = or)) + scale_colour_viridis(option = "plasma",direction = -1)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))          
                
                return(res_tot)
                
}

#
# Define a function to define is a marker is present or not and create a barplot
#

percentage_markers_for_barplot<-function(input_ml2,classes,markers_to_use,cnv_thr=0.7){
  
  subset_ml<-input_ml2[which(input_ml2$biological_states%in%classes),which(colnames(input_ml2)%in%c("biological_states",markers_to_use))]
  
  nsamplesgroup1<-table(subset_ml$biological_states)[classes[1]]
  nsamplesgroup2<-table(subset_ml$biological_states)[classes[2]]
  
  #melt everything
  melt_subset_ml<-melt(subset_ml)
  melt_subset_ml$value2<-rep(0,nrow(melt_subset_ml))
  #is amplified
  melt_subset_ml$value2[melt_subset_ml$value>=cnv_thr]<-1 #is amplified
  #is depleted
  melt_subset_ml$value2[melt_subset_ml$value<= -cnv_thr]<-2 #is depleted, i used 2 to indicate depletion, to avoid problems with columns names after
  #just to be sure, change the values also for the mutations
  melt_subset_ml[grep(melt_subset_ml$variable,pattern="dndscv"),"value2"]<-  melt_subset_ml[grep(melt_subset_ml$variable,pattern="dndscv"),"value"]
  
  #split the data.frame for variable (marker, e.g focal events)
  list_subset_ml<-split(melt_subset_ml,melt_subset_ml$variable)
  
  # define a function to compute in each data.frame of the list the frequencies of the samples in the groups and with amplification/depletion or not
  table_custom<-function(X){
    
                freq_for_marker<-table(X[,c(1,4)])
                freq_row<-data.frame(freq_for_marker)[,3] # the frequencies
                groups_found<-data.frame(freq_for_marker)[,1] # the groups
                status_found<-data.frame(freq_for_marker)[,2] # the status (0 no, 1 amplified or mutated, 2 depleted)
                
                names(freq_row)<-paste(groups_found,status_found,sep="_")
                
                return(data.frame(t(freq_row)))
    }
  
  resFreq<-do.call(rbind.fill,lapply(list_subset_ml,table_custom))
  df_freq2<-cbind(names(list_subset_ml),resFreq)
  df_freq2[is.na(df_freq2)]<-0
  
  # df_freq3<-round(df_freq2[,-1]/apply(df_freq2[,-1],1,sum),3)
  # df_freq4<-cbind(markers=df_freq2[,1],df_freq3)
  colnames(df_freq2)[1]<-"mks"
  df_freq5<-df_freq2[,grep(colnames(df_freq2),pattern="_0",invert="T")]
  
  #compute % by group
  freq_group1<-df_freq5[,grep(colnames(df_freq5),pattern=classes[1])]/nsamplesgroup1
  freq_group2<-df_freq5[,grep(colnames(df_freq5),pattern=classes[2])]/nsamplesgroup2
  
  df_freq6<-data.frame(df_freq2[,1],cbind(freq_group1,freq_group2))
  colnames(df_freq6)[1]<-"markers"
  
  df_freq_melt<-melt(df_freq6)
  markers_final<-df_freq_melt[,1]
  groups_final<-sapply(strsplit(as.character(df_freq_melt[,2]),split="_"),"[[",1)
  status_final<-sapply(strsplit(as.character(df_freq_melt[,2]),split="_"),"[[",2)
  freq_final<-df_freq_melt[,3]
  
  df_final<-data.frame(markers=markers_final,
                       groups_final=groups_final,
                       status_final=status_final,
                       freq_final=freq_final)
  return(df_final)
  
}

#
# Define a function for chart
#

integration_multichart<-function(coef_data_matrix,list_res_comparison,genomic_matrix=input_gf_tme_fb2,input_ml,input_ml2,classes=classes,cnv_thr=0.7,manual_order=F,variables_sorted=NULL,height=10,width=16,output_file="cond1_vs_cond2_Coef_HistFreq_GenomicMarkers.pdf"){
  
  # coef_data_matrix$new_labels<-coef_data_matrix$coef
  # coef_data_matrix$new_labels<-gsub(coef_data_matrix$new_labels,pattern="_dndscv",replacement="")
  # coef_data_matrix$new_labels<-gsub(coef_data_matrix$new_labels,pattern="_focal",replacement="")
  # 
  #
  # Step1: boxplot of the coefficients
  #
  
  coef_data_matrix_subset<-coef_data_matrix[which(coef_data_matrix$coef%in%list_res_comparison[[6]]),]
  coef_uniq_mu<-unique(coef_data_matrix_subset[,c(1,3)])
  
  if(manual_order==FALSE){
    
  variables_sorted<-coef_uniq_mu[order(coef_uniq_mu$mu,decreasing=T),1]
  
  } else{
    
  variables_sorted<-variables_sorted
    
  }
  
  coef_data_matrix_subset$coef<-factor(coef_data_matrix_subset$coef,levels=rev(variables_sorted))
    
  p<-ggplot(coef_data_matrix_subset,aes(x=coef, y=values,fill=color))+
    geom_boxplot()+scale_fill_manual(values=c(positive="firebrick1",negative="steelblue"))+scale_y_continuous(limits = c(min(coef_data_matrix_subset$values), max(coef_data_matrix_subset$values)), breaks=seq(min(coef_data_matrix_subset$values),max(coef_data_matrix_subset$values),by=1), labels=round(seq(min(coef_data_matrix_subset$values),max(coef_data_matrix_subset$values),by=1)))+
    geom_hline(yintercept = 0, linetype="dotted",color = "red", size=1)+
    theme(text = element_text(size=8))+coord_flip()+theme_bw()+ theme(legend.position="top")+labs(title="Boxplots")
  

  #
  # Step2: Prepare the data for the barplot
  #
  
  coef_uniq_mu<-unique(coef_data_matrix_subset[,c(1,3)])
  
  if(manual_order==FALSE){
    
  variables_sorted<-coef_uniq_mu[order(coef_uniq_mu$mu,decreasing=T),1]
  
  }else{
  
  variables_sorted<-variables_sorted
    
  }
  
  freq_barplot_cond1_vs_cond2<-percentage_markers_for_barplot(input_ml2,classes,markers_to_use=intersect(list_res_comparison[[6]],coef_data_matrix_subset$coef),cnv_thr=cnv_thr)
  freq_barplot_cond1_vs_cond2<-freq_barplot_cond1_vs_cond2[freq_barplot_cond1_vs_cond2$markers%in%intersect(list_res_comparison[[6]],coef_data_matrix_subset$coef),]
  
  freq_barplot_cond1_vs_cond2$markers<-factor(freq_barplot_cond1_vs_cond2$markers,levels=rev(variables_sorted))
  
  #
  # Step3: Create a column to highlight if a feature is mutation, copy number or others
  #
  
  coef_data_matrix_subset$type_element<-rep(0,nrow(coef_data_matrix_subset))
  coef_data_matrix_subset$type_element[grep(coef_data_matrix_subset$coef,pattern="focal")]<-1
  coef_data_matrix_subset$type_element[grep(coef_data_matrix_subset$coef,pattern="dndscv")]<-2
  coef_data_matrix_subset$type_element[grep(coef_data_matrix_subset$coef,pattern="arm")]<-3
  coef_data_matrix_subset$type_element[grep(coef_data_matrix_subset$coef,pattern="aneuploidy")]<-4
  
  coef_data_matrix_subset$type_element2<-as.factor(coef_data_matrix_subset$type_element)
  
  coef_data_matrix_subset_forann<-coef_data_matrix_subset[,colnames(coef_data_matrix_subset)%in%c("coef","type_element2")]
  
  coef_data_matrix_subset_forann$coef<-factor(coef_data_matrix_subset_forann$coef,levels=rev(variables_sorted))
  
  if(length(table(coef_data_matrix_subset$type_element))==3){
    
  annotation_elements_cond1_vs_cond2 <- ggplot(coef_data_matrix_subset_forann, aes(x = 1, y = coef))+
                                        geom_point(aes(shape =type_element2,color=type_element2),size=3)+scale_shape_manual(values=c(15,17,18))  + 
                                        scale_color_manual(values=c("#fdc500", "#577590", "#f94144"))+
                                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
                                        theme(legend.position="top",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(title="annotation")
  
  annotation_elements_cond1_vs_cond2_debug <- ggplot(coef_data_matrix_subset_forann, aes(x = 1, y = coef))+
                                              geom_point(aes(shape =type_element2,color=type_element2),size=3)+scale_shape_manual(values=c(15,17,18)) + 
                                              scale_color_manual(values=c("#fdc500", "#577590", "#f94144"))+
                                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")+
                                              theme_bw()+labs(title="annotation")
  
  }
  
  if(length(table(coef_data_matrix_subset$type_element))==4){
    
    annotation_elements_cond1_vs_cond2 <- ggplot(coef_data_matrix_subset_forann, aes(x = 1, y = coef))+
      geom_point(aes(shape =type_element2,color=type_element2),size=3)+scale_shape_manual(values=c(15,17,18,16))  + 
      scale_color_manual(values=c("#fdc500", "#577590", "#f94144","#bd9eff"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
      theme(legend.position="top",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(title="annotation")
    
    annotation_elements_cond1_vs_cond2_debug <- ggplot(coef_data_matrix_subset_forann, aes(x = 1, y = coef))+
      geom_point(aes(shape =type_element2,color=type_element2),size=3)+scale_shape_manual(values=c(15,17,18,16)) + 
      scale_color_manual(values=c("#fdc500", "#577590", "#f94144","#bd9eff"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="top")+
      theme_bw()+labs(title="annotation")
    
  }
  
  
  barplot_cond1_vs_cond2 <- ggplot(freq_barplot_cond1_vs_cond2, aes(x = markers, y = freq_final))+
                            geom_col(aes(fill = status_final), width = 0.7)+theme_bw()+
                            facet_wrap(~groups_final) + coord_flip()+
                            theme(legend.position="top",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  barplot_cond1_vs_cond2_debug <- ggplot(freq_barplot_cond1_vs_cond2, aes(x = markers, y = freq_final))+
                                  geom_col(aes(fill = status_final), width = 0.7)+theme_bw()+
                                  facet_wrap(~groups_final) + coord_flip()+
                                  theme(legend.position="top")
  
  fold_freq<-do.call(cbind,split(freq_barplot_cond1_vs_cond2,freq_barplot_cond1_vs_cond2$groups_final))
  
  column_class1_fold<-paste(classes[2],".freq_final",sep="")
  column_class2_fold<-paste(classes[1],".freq_final",sep="")
  
  column_class1_status_final<-paste(classes[2],".status_final",sep="")
  
  fold_freq2<-data.frame(markers=fold_freq[,1],
                         status_final=fold_freq[,column_class1_status_final],
                         freq_fold= fold_freq[,column_class1_fold] - fold_freq[,column_class2_fold])
  
  fold_freq2$direction<-ifelse(fold_freq2$freq_fold>0,"positive_fold","negative_fold")
  
  freqgp<-fold_freq2 %>%
    ggplot(aes(x = markers, y = freq_fold, fill = status_final))+
    geom_col(stat = "identity")+
    coord_flip()+theme_bw()+theme(legend.position="top",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+ggtitle("FoldFreq")
  
  freqgp_debug<-fold_freq2 %>%
    ggplot(aes(x = markers, y = freq_fold, fill = status_final))+
    geom_col(stat = "identity")+
    coord_flip()+theme_bw()
  
  X=list_res_comparison
  
  #
  # Step4: Create a balloonchart
  #
  
  ballonchart_cond1_vs_cond2<-multiOmicsCharts(X,genomic_matrix,input_ml,cnv_thr=0.7,suffix_output="markers_gs_december",heigth=15)
  
  ballonchart_cond1_vs_cond2<-ballonchart_cond1_vs_cond2[ballonchart_cond1_vs_cond2$genes%in%intersect(list_res_comparison[[6]],coef_data_matrix_subset$coef),]
  
  ballonchart_cond1_vs_cond2$genes<-factor(ballonchart_cond1_vs_cond2$genes,levels=rev(variables_sorted))
  
  
  ballonchart_cond1_vs_cond2_to_plot <- ggplot(ballonchart_cond1_vs_cond2, aes(x = feature_value, y = genes)) +
                                        geom_point(aes(size = pval, shape =status_sig,colour = or))+scale_shape_manual(values=c(26, 18))+ scale_colour_viridis(option = "plasma",direction = -1)+
                                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
                                        theme(legend.position="top",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
                                      
  ballonchart_cond1_vs_cond2_to_plot_debug <- ggplot(ballonchart_cond1_vs_cond2, aes(x = feature_value, y = genes)) +
                                              geom_point(aes(size = pval, shape =status_sig,colour = or))+scale_shape_manual(values=c(26, 18))+ scale_colour_viridis(option = "plasma",direction = -1)+
                                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
                                              theme(legend.position="top")
  
  pdf(output_file,height=height,width=width)
  p4<-ggarrange(p,annotation_elements_cond1_vs_cond2,barplot_cond1_vs_cond2,freqgp,ballonchart_cond1_vs_cond2_to_plot, widths= c(0.4,0.05,0.20,0.10,0.1),nrow=1,ncol=5)
  print(p4)
  p5<-ggarrange(p,annotation_elements_cond1_vs_cond2_debug,barplot_cond1_vs_cond2_debug,freqgp_debug,ballonchart_cond1_vs_cond2_to_plot_debug, widths= c(1,0.1,0.1,0.1,1),nrow=1,ncol=5)
  print(p5)
  print(ballonchart_cond1_vs_cond2_to_plot_debug)
  dev.off()
  
}  
setwd("/home/data/pseudospace/pathway_characterization/broad_genomic_features")

# 
# Build plot MES-EPI 
#

input_mes_vs_epi<-vector(mode="list",6)

input_mes_vs_epi[[1]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_mes_vs_epi[[2]]<-setdiff(mut_genes_mes_vs_epi,mut_genes_pemt_vs_epi)
input_mes_vs_epi[[3]]<-setdiff(cnv_markers_mes_vs_epi,cnv_markers_pemt_vs_epi)
input_mes_vs_epi[[4]]<-c("epi","mes")
input_mes_vs_epi[[5]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_mes_vs_epi[[6]]<-raw_markers_specifics_mes_vs_epi
names(input_mes_vs_epi)<-c("columns_to_use","mutations","copy_number_alterations","comparison","boxplot_features","raw_markers_specifics")

integration_multichart(coef_data_matrix=coef_mes_vs_epi2,
                       list_res_comparison=input_mes_vs_epi,
                       genomic_matrix=input_gf_tme_fb2,
                       input_ml=input_ml,
                       input_ml2=input_ml2,
                       classes=c("epi","mes"),
                       cnv_thr=0.7,
                       output_file="MES_vs_EPI_Coef_HistFreq_GenomicMarkers_rawdndscv.pdf")


input_mes_vs_epi<-vector(mode="list",6)

input_mes_vs_epi[[1]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_mes_vs_epi[[2]]<-setdiff(mut_genes_mes_vs_epi,mut_genes_pemt_vs_epi)
input_mes_vs_epi[[3]]<-setdiff(cnv_markers_mes_vs_epi,cnv_markers_pemt_vs_epi)
input_mes_vs_epi[[4]]<-c("epi","mes")
input_mes_vs_epi[[5]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_mes_vs_epi[[6]]<-raw_markers_specifics_mes_vs_epi
names(input_mes_vs_epi)<-c("columns_to_use","mutations","copy_number_alterations","comparison","boxplot_features","raw_markers_specifics")

integration_multichart(coef_data_matrix=coef_mes_vs_epi2[which(coef_mes_vs_epi2$mu>=1 | coef_mes_vs_epi2$mu<= -1),],
                       list_res_comparison=input_mes_vs_epi,
                       genomic_matrix=input_gf_tme_fb2,
                       input_ml=input_ml,
                       input_ml2=input_ml2,
                       classes=c("epi","mes"),
                       cnv_thr=0.7,
                       manual_order = F,
                       height= 4,
                       width = 10, 
                       output_file="MES_vs_EPI_Coef_HistFreq_GenomicMarkers_v2_rawdndscv.pdf")

# 
# Build plot  pEMTvsEPI 
#

input_epi_vs_pemt<-vector(mode="list",6)

input_epi_vs_pemt[[1]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_epi_vs_pemt[[2]]<-setdiff(mut_genes_pemt_vs_epi,mut_genes_mes_vs_epi)
input_epi_vs_pemt[[3]]<-setdiff(cnv_markers_pemt_vs_epi,cnv_markers_mes_vs_epi)
input_epi_vs_pemt[[4]]<-c("epi","mix")
input_epi_vs_pemt[[5]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_epi_vs_pemt[[6]]<-raw_markers_specifics_pemt_vs_epi



integration_multichart(coef_data_matrix=coef_coef_pemt_vs_epi2[which(coef_coef_pemt_vs_epi2$mu>=1 | coef_coef_pemt_vs_epi2$mu<= -1),],
                       list_res_comparison=input_epi_vs_pemt,
                       genomic_matrix=input_gf_tme_fb2,
                       input_ml=input_ml,
                       input_ml2=input_ml2,
                       classes=c("epi","mix"),
                       cnv_thr=0.7,
                       output_file="pEMT_vs_EPI_Coef_HistFreq_GenomicMarkers_rawdndscv.pdf")
      
input_epi_vs_pemt<-vector(mode="list",6)

input_epi_vs_pemt[[1]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_epi_vs_pemt[[2]]<-setdiff(mut_genes_pemt_vs_epi,mut_genes_mes_vs_epi)
input_epi_vs_pemt[[3]]<-setdiff(cnv_markers_pemt_vs_epi,cnv_markers_mes_vs_epi)
input_epi_vs_pemt[[4]]<-c("epi","mix")
input_epi_vs_pemt[[5]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_epi_vs_pemt[[6]]<-raw_markers_specifics_pemt_vs_epi

integration_multichart(coef_data_matrix=coef_coef_pemt_vs_epi2[which(coef_coef_pemt_vs_epi2$mu>=1.5 | coef_coef_pemt_vs_epi2$mu<= -1.5),],
                       list_res_comparison=input_epi_vs_pemt,
                       genomic_matrix=input_gf_tme_fb2,
                       input_ml=input_ml,
                       input_ml2=input_ml2,
                       classes=c("epi","mix"),
                       cnv_thr=0.7,
                       height= 5,
                       width = 10, 
                       output_file="pEMT_vs_EPI_Coef_HistFreq_GenomicMarkers_V2_rawdndscv.pdf")
# 
# Common between MES vs EPI and pEMT vs EPI
#

input_common<-vector(mode="list",6)

input_common[[1]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_common[[2]]<-intersect(mut_genes_pemt_vs_epi,mut_genes_mes_vs_epi)
input_common[[3]]<-intersect(cnv_markers_pemt_vs_epi,cnv_markers_mes_vs_epi)
input_common[[4]]<-c("epi","mix")
input_common[[5]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_common[[6]]<-raw_markers_common


integration_multichart(coef_data_matrix=coef_coef_pemt_vs_epi2,
                       list_res_comparison=input_common,
                       genomic_matrix=input_gf_tme_fb2,
                       input_ml=input_ml,
                       input_ml2=input_ml2,
                       classes=c("epi","mix"),
                       cnv_thr=0.7,
                       output_file="Common_Coef_HistFreq_GenomicMarkers_V1_rawdndscv.pdf")


input_common<-vector(mode="list",6)

input_common[[1]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_common[[2]]<-intersect(mut_genes_pemt_vs_epi,mut_genes_mes_vs_epi)
input_common[[3]]<-intersect(cnv_markers_pemt_vs_epi,cnv_markers_mes_vs_epi)
input_common[[4]]<-c("epi","mes")
input_common[[5]]<-c("score_emt","Buffa_hypoxia_genes","mRNAsi","AneuploidyScore.AS.","CA20")
input_common[[6]]<-raw_markers_common


order_variables<-c("11p_arm_deletion","16q_arm_deletion","22q_arm_deletion","15q_arm_deletion","HRAS_dndscv","13q_arm_deletion","9q_arm_deletion",
                   "RB1_focal","FLCN_focal","RB1_dndscv","17p_arm_deletion","TP53_dndscv",
                   "PIK3CA_dndscv","KRAS_dndscv","10p_arm_amplification","GATA3_focal","CCNE1_focal",
                   "20q_arm_amplification","SDHC_focal","NOTCH2_focal","8q_arm_amplification","SPOP_dndscv",
                   "1q_arm_amplification","CDKN2A_focal","17q_arm_amplification","ASXL1_focal",
                   "CCND1_dndscv","NRAS_dndscv","JAK1_dndscv","ZFHX3_dndscv","RNF43_dndscv",
                   "FGFR2_dndscv","CHD4_dndscv","FGFR3_dndscv")
  
integration_multichart(coef_data_matrix=coef_mes_vs_epi2,
                       list_res_comparison=input_common,
                       genomic_matrix=input_gf_tme_fb2,
                       input_ml=input_ml,
                       input_ml2=input_ml2,
                       classes=c("epi","mes"),
                       cnv_thr=0.7,
                       manual_order = TRUE,
                       variables_sorted = order_variables,
                       output_file="Common_Coef_HistFreq_GenomicMarkers_V2_rawdndscv.pdf")



# freq_barplot_mes_vs_epi<-percentage_markers_for_barplot(input_ml2,c("mes","epi"),markers_to_use=input_mes_vs_epi[[6]],cnv_thr=0.7)
