#
# Read the depmap data
#

library(depmap)
library(data.table)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

load("CCLE_pseudotime_MCF.RData")

annotation_drug<-read.csv(file="primary-screen-replicate-collapsed-treatment-info.csv",stringsAsFactors = F)
annotation_cellline<-read.csv(file="primary-screen-cell-line-info.csv",stringsAsFactors = F)

tab_input<-fread(file="primary-screen-replicate-collapsed-logfold-change.csv")
colnames(tab_input)[1]<-"code_cellline"

annotation_cellline2<-merge(annotation_cellline[,c(2:4)],tab_input,by.x="depmap_id",by.y="code_cellline")

CCLE_with_drugs<-merge(CCLE_pseudotime_EMT_metpot,annotation_cellline2,by.x="CCLE_ID",by.y="ccle_name",all.x=F,all.y=F)
CCLE_with_drugs2<-CCLE_with_drugs[,c(9,12:ncol(CCLE_with_drugs))]

ccle_mutation<-fread(file="CCLE_mutations.csv")
  
#
# Now integrate the genomic features
#

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

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T)
mut_genes_mes_vs_epi<-sapply(strsplit(grep(markers_mes_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
mut_genes_pemt_vs_epi<-sapply(strsplit(grep(markers_pemt_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)

setwd(output_dir)

list_analysis<-vector(mode="list",2)
list_analysis[[1]]<-setdiff(mut_genes_mes_vs_epi,mut_genes_pemt_vs_epi)
list_analysis[[2]]<-setdiff(mut_genes_pemt_vs_epi,mut_genes_mes_vs_epi)

names(list_analysis)<-c("mes_vs_epi","pemt_vs_epi")

AOVlist<-function(X,status_mut){
  
  X2<-data.frame(cbind(X,status_mut))
  colnames(X2)<-c("drug","status_mut")
  
  res.aov2 <- aov(drug ~ status_mut, data = X2)
  pvalue_aov = unlist(summary(res.aov2))[["Pr(>F)1"]]
  
  nomut<-mean(X2[X2[,2]%in%0,1],na.rm=T)
  mut<-mean(X2[X2[,2]%in%1,1],na.rm=T)
  lfc<-mut-nomut
  
  res<-c(pvalue_aov,mut,nomut,lfc)
  names(res)<-c("pvalue","mut","nomut","lfc")
  
  return(res)
  
}

results_aov<-vector(mode="list",length(list_analysis))

for(i in 1:length(list_analysis)){
  
  print(i)
  
  current_mut_genes<-list_analysis[[i]]
  
  current_comparison<-names(list_analysis)[i]
  
  samples_with_mutations<-data.frame(status_mut=1,ccle_mutation[ccle_mutation$Hugo_Symbol%in%current_mut_genes,c(1,16)])
  
  samples_without_mutations<-data.frame(status_mut=0,ccle_mutation[-which(ccle_mutation$Hugo_Symbol%in%current_mut_genes),c(1,16)])
  
  samples_without_mutations$Hugo_Symbol<-""
  
  samples_mutation_status_all<-rbind(samples_with_mutations,samples_without_mutations)
  
  samples_mutation_status_all2<-merge(annotation_cellline,samples_mutation_status_all,by.x="depmap_id",by.y="DepMap_ID")
  
  CCLE_pseudotime_EMT_HMM$CCLE_ID<-as.character(CCLE_pseudotime_EMT_HMM$CCLE_ID)
  
  samples_mutation_status_all2_with_emt<-merge(samples_mutation_status_all2,CCLE_pseudotime_EMT_metpot,by.x="ccle_name",by.y="CCLE_ID",all.x=F,all.y=F)
  
  samples_mutation_status_all2_with_emt<-samples_mutation_status_all2_with_emt[,c(1,2,4,11,12,14,20)]
  
  samples_mutation_status_all2_with_emt<-unique(samples_mutation_status_all2_with_emt)
  
  # samples_mutation_status_all2_with_emt_and_drugs<-merge(samples_mutation_status_all2_with_emt,tab_input,by.x="depmap_id",by.y="code_cellline")
  
  #
  # See if the emt scores and mutations are able to predict the drug status of each 
  #
  gene_symbol_to_analyze<-unique(samples_mutation_status_all2_with_emt$Hugo_Symbol)
  gene_symbol_to_analyze<-gene_symbol_to_analyze[!gene_symbol_to_analyze==""]
  
  df_all_genes<-data.frame()
  
  for(g in gene_symbol_to_analyze){
    
  print(g)
    
  df_genes<-samples_mutation_status_all2_with_emt[samples_mutation_status_all2_with_emt$Hugo_Symbol%in%g,]
  df_nogenes<-samples_mutation_status_all2_with_emt[-which(samples_mutation_status_all2_with_emt$Hugo_Symbol%in%g),]
  
  df_tot_genes<-rbind(df_genes,df_nogenes)
  
  df_tot_genes[-which(df_tot_genes$Hugo_Symbol%in%g),"status_mut"]<-0
  
  df_genes_and_drugs<-merge(df_tot_genes,tab_input,by.x="depmap_id",by.y="code_cellline",all.x=F,all.y=F)
  
  only_drugs_df<-df_genes_and_drugs[,8:ncol(df_genes_and_drugs)]
  
  res_temp<-t(apply(only_drugs_df,2,AOVlist,df_genes_and_drugs$status_mut))
  
  res_temp2<-data.frame(genes=g,drugs=rownames(res_temp),res_temp,padjust=p.adjust(res_temp[,"pvalue"],"BH"))
  
  df_all_genes<-rbind(df_all_genes,res_temp2)
  
  }
  
  results_aov[[i]]<-df_all_genes

}

save(results_aov,file="anova_depmap_featuresML.RData")

#
# Volcano plot
#
names(results_aov)[1]<-c("mes_vs_epi")
names(results_aov)[2]<-c("pEMT_vs_epi")

pdf("volcano_plot_aov_0.01_depmap_mlfeatures.pdf")
for(vp in 1:length(results_aov)){
  
  current_comparison<-names(results_aov)[vp]
  
  temp_aov_vp<-results_aov[[vp]]
  
  temp_aov_vp2<-temp_aov_vp[temp_aov_vp$pvalue<=0.01,]
  
  temp_aov_vp3<-merge(annotation_drug,temp_aov_vp2,by.x="column_name",by.y="drugs")
  
  temp_aov_vp3$status<-rep("nosig",nrow(temp_aov_vp3))
  temp_aov_vp3$status[which(temp_aov_vp3$lfc>=0.5)]<-"up"
  temp_aov_vp3$status[which(temp_aov_vp3$lfc<=-0.5)]<-"dw"
  
  temp_aov_vp3$label<-paste(temp_aov_vp3$genes,"\n",temp_aov_vp3$name)
  
  p2<-ggplot(temp_aov_vp3, aes(x = lfc, y = -log10(pvalue))) +
    geom_point(aes(color = status)) +
    scale_color_manual(values = c("blue2", "grey","red2")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel( data = subset(temp_aov_vp3, abs(lfc)>0.5),
                     aes(label = label),
                     size = 5,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"),max.overlaps=20
    )+ geom_vline(xintercept = c(-0.5,0,0.5),linetype="dotted",size=0.5)+labs(title=current_comparison)
  
  print(p2)
  
}
dev.off()
