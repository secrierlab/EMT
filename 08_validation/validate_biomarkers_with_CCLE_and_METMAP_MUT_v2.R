library(data.table)
library(openxlsx)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(vcd)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)

RadarTheme<-theme(
panel.background = element_rect(fill = "white", colour = "white", linetype = "black"),
panel.grid.major = element_line(size = 0.25, linetype = 'black', colour = "white"), 
panel.grid.minor = element_line(size = 0, linetype = 'black',colour = "black"))


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
coef_mes_vs_epi<-getCoefSS(res_lasso,c(mut_genes_mes_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_mes_vs_epi2<-merge(coef_mes_vs_epi,mu,by="coef")
coef_mes_vs_epi2[,1]<-gsub(coef_mes_vs_epi2[,1],pattern="^X",replacement="")

mut_genes_mes_vs_epi_raw<-intersect(mut_genes_mes_vs_epi_raw,unique(coef_mes_vs_epi2[coef_mes_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
mut_genes_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T)
mut_genes_pemt_vs_epi<-sapply(strsplit(grep(markers_pemt_vs_epi,pattern="dndscv",value=T),split="_"),"[[",1)
coef_pemt_vs_epi<-getCoefSS(res_lasso,c(mut_genes_pemt_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_pemt_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_pemt_vs_epi2<-merge(coef_pemt_vs_epi,mu,by="coef")
coef_coef_pemt_vs_epi2[,1]<-gsub(coef_coef_pemt_vs_epi2[,1],pattern="^X",replacement="")

mut_genes_pemt_vs_epi_raw<-intersect(mut_genes_pemt_vs_epi_raw,unique(coef_coef_pemt_vs_epi2[coef_coef_pemt_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_mix<-getFeaturesSS(res_lasso,ss=800)
mut_genes_mes_vs_mix_raw<-grep(markers_mes_vs_mix,pattern="dndscv",value=T)
mut_genes_mes_vs_mix<-sapply(strsplit(grep(markers_mes_vs_mix,pattern="dndscv",value=T),split="_"),"[[",1)
coef_mes_vs_mix<-getCoefSS(res_lasso,c(mut_genes_mes_vs_mix_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_mix, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_mes_vs_mix2<-merge(coef_mes_vs_mix,mu,by="coef")
coef_coef_mes_vs_mix2[,1]<-gsub(coef_coef_mes_vs_mix2[,1],pattern="^X",replacement="")

mut_genes_mes_vs_mix_raw<-intersect(mut_genes_mes_vs_mix_raw,unique(coef_coef_mes_vs_mix2[coef_coef_mes_vs_mix2$color%in%"positive",1]))


raw_markers_specifics_mes_vs_epi<-setdiff(c(mut_genes_mes_vs_epi_raw),c(mut_genes_pemt_vs_epi_raw,mut_genes_mes_vs_mix_raw))
raw_markers_specifics_pemt_vs_epi<-setdiff(y=c(mut_genes_mes_vs_epi_raw,mut_genes_mes_vs_mix_raw),x=c(mut_genes_pemt_vs_epi_raw))
raw_markers_common<-intersect(c(mut_genes_mes_vs_epi_raw),c(mut_genes_pemt_vs_epi_raw))

raw_markers_specifics_mes_vs_mix<-setdiff(x=c(mut_genes_mes_vs_mix_raw),y=c(mut_genes_mes_vs_epi_raw,mut_genes_pemt_vs_epi_raw))

list_common_all<-list(c(mut_genes_mes_vs_epi_raw),c(mut_genes_mes_vs_mix_raw),c(mut_genes_pemt_vs_epi_raw))

common_between_three_groups<-Reduce(intersect,list_common_all)

#
# Upload the biomarkers
#

tab_sampleinfo<-fread("primary-screen-cell-line-info.csv",data.table=F)

#
# Upload MET500
#

brain_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 1)
colnames(brain_metmap500)[1]<-"cell_line"

lung_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 2)
colnames(lung_metmap500)[1]<-"cell_line"

liver_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 3)
colnames(liver_metmap500)[1]<-"cell_line"

bone_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 4)
colnames(bone_metmap500)[1]<-"cell_line"

kidney_metmap500<-read.xlsx("MetMap500_met_potential.xlsx",sheet = 5)
colnames(kidney_metmap500)[1]<-"cell_line"

met_potential_list<-vector(mode="list",length=5)
met_potential_list[[1]]<-brain_metmap500
met_potential_list[[2]]<-lung_metmap500
met_potential_list[[3]]<-liver_metmap500
met_potential_list[[4]]<-bone_metmap500
met_potential_list[[5]]<-kidney_metmap500
names(met_potential_list)<-c("brain","lung","liver","bone","kidney")


met_pot_all<-do.call(rbind,met_potential_list)
met_pot_all2<-data.frame(tissue=rownames(met_pot_all),met_pot_all)
met_pot_all2$tissue<-sapply(strsplit(as.character(met_pot_all2$tissue),split="\\."),"[[",1)

codes_to_use_metmap<-tab_sampleinfo[which(tab_sampleinfo[,3]%in%met_pot_all2$cell_line),2]


#
# Load cna data CCLE
#

tab_mut_ccle<-fread(file="CCLE_mutations.csv",header=T,fill = TRUE,data.table=F)
tab_mut_ccle_subselect<-tab_mut_ccle[tab_mut_ccle$DepMap_ID%in%codes_to_use_metmap,]
tab_mut_ccle_subselect<-tab_mut_ccle_subselect[which(tab_mut_ccle_subselect$Variant_annotation%in%"damaging"),]
colnames(tab_mut_ccle_subselect)[16]<-"Samples"

#
# See the effect of the CNA in CCLE and on the metastatic potential
#

list_features<-vector(mode="list",7)
list_features[[1]]<-raw_markers_specifics_mes_vs_epi
list_features[[2]]<-raw_markers_specifics_pemt_vs_epi
list_features[[3]]<-raw_markers_specifics_mes_vs_mix
list_features[[4]]<-raw_markers_common
list_features[[5]]<-common_between_three_groups
list_features[[6]]<-c("EGFR_dndscv","SETD2_dndscv","RBM10_dndscv","CDKN2A_dndscv","RB1_dndscv","MGA_dndscv","ZIC1_dndscv","NF1_dndscv","REG3A_dndscv","ARID2_dndscv","ZFP36L1_dndscv","ARID1A_dndscv")
list_features[[7]]<-c("CBFB_dndscv","NCOR1_dndscv","AKT1_dndscv","SHISA4_dndscv","NF1_dndscv","FOXA1_dndscv")

names(list_features)<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common","common_all","luad","breast")

setwd(output_dir)

stat.test_sig_all<-data.frame()
stat.test_sig_all_pancancer<-data.frame()
pancancerboxplot_tot<-data.frame()

for(i in 1:length(list_features)){
  
  print(i)
  
  name_comparison<-names(list_features)[i]
  features<-sapply(strsplit(grep(list_features[[i]],pattern="dndscv",value=T),split="\\_"),"[[",1)
  
  tab_mut_features<-tab_mut_ccle_subselect[which(tab_mut_ccle_subselect$Hugo_Symbol %in%features),]
  
  tab_mut_features2<-melt(table(tab_mut_features$Samples, tab_mut_features$Hugo_Symbol))
  tab_mut_features2[,1]<-as.character(tab_mut_features2[,1])
  colnames(tab_mut_features2)[1]<-"ID"
  
  tab_mut_with_ann<-merge(tab_sampleinfo[,c(2,3,4)],tab_mut_features2,by.x="depmap_id",by.y="ID")
  tab_mut_with_ann_and_metpot<-merge(met_pot_all2,tab_mut_with_ann,by.x="cell_line",by.y="ccle_name")
  
  tab_mut_with_ann_and_metpot$value2<-factor(ifelse(tab_mut_with_ann_and_metpot$value>=1,"withMUT","withoutMUT"),levels=c("withoutMUT","withMUT"))
  
  input_for_stat_pancancer<-tab_mut_with_ann_and_metpot 
  freq_genes_status_pancancer<- table(input_for_stat_pancancer$Var2,input_for_stat_pancancer$value2)
  genes_to_consider_pancancer<- rownames(freq_genes_status_pancancer[freq_genes_status_pancancer[,1]>1 & freq_genes_status_pancancer[,2]>1,])

  stat.test_pancancer <-input_for_stat_pancancer[which(input_for_stat_pancancer$Var2 %in%genes_to_consider_pancancer),]%>%
    group_by(Var2) %>%
    wilcox_test(mean ~ value2)
  
  stat.test_pancancer_df_temp<-data.frame(stat.test_pancancer)
  
  input_for_mean_pancancer<- input_for_stat_pancancer[which(input_for_stat_pancancer$Var2 %in%genes_to_consider_pancancer),]
  
  input_for_mean_pancancer$mean_unlog<-10^input_for_mean_pancancer$mean
  
  stat_means_pancancer<-aggregate(mean_unlog ~ Var2 + value2 , data = input_for_mean_pancancer, FUN = median, na.rm = TRUE)
  
  logFCmut_no_mut<-stat_means_pancancer[stat_means_pancancer$value2%in%"withMUT",3]/stat_means_pancancer[stat_means_pancancer$value2%in%"withoutMUT",3]
  
  logFCmut_no_mut_df<-data.frame(genes=stat_means_pancancer[stat_means_pancancer$value2%in%"withMUT","Var2"],logFCmut_no_mut=logFCmut_no_mut)
  
  res_stat_with_mean_pancancer <- merge(stat.test_pancancer_df_temp,logFCmut_no_mut_df,by.x="Var2",by.y="genes")
  
  stat.test_pancancer<-res_stat_with_mean_pancancer[res_stat_with_mean_pancancer$p<=0.05,]
  
  stat.test_sig2_pancancer<-res_stat_with_mean_pancancer[which(res_stat_with_mean_pancancer$p<=0.05),]
  print(paste("The number of rows in the pan-cancer analysis is:",nrow(stat.test_sig2_pancancer)))
  
  pancancerboxplot<-input_for_stat_pancancer[which(input_for_stat_pancancer$genes%in%stat.test_sig2_pancancer[,1]),]
  
  #
  # Compare the metastatic potential
  #
  res_stat_df<-data.frame()
 
  for(tss in unique(tab_mut_with_ann_and_metpot$primary_tissue)){
    
  print(tss)
    
  input_for_stat<-tab_mut_with_ann_and_metpot[which(tab_mut_with_ann_and_metpot$primary_tissue%in%tss),] 
  freq_genes_status<- table(input_for_stat$Var2,input_for_stat$value2)
  genes_to_consider<- rownames(freq_genes_status[freq_genes_status[,1]>1 & freq_genes_status[,2]>1,])

  if(length(genes_to_consider)!=0){
    
  stat.test <-input_for_stat[which(input_for_stat$Var2 %in%genes_to_consider),] %>%
    group_by(primary_tissue, Var2) %>%
    wilcox_test(mean ~ value2)
  
  stat.test<-stat.test[which(stat.test$n1 >= 10 & stat.test$n2 >= 10),]
  
  if(nrow(stat.test)>=1){
    
  stats_df_temp<-data.frame(stat.test)
  
  input_for_mean<- input_for_stat[which(input_for_stat$Var2 %in%genes_to_consider),]
  input_for_mean$mean_unlog<-10^input_for_mean$mean
    
  stat_means<-aggregate(mean_unlog ~ primary_tissue + Var2 + value2 , data = input_for_mean, FUN = median, na.rm = TRUE)
  
  logFCmut_no_mut<-stat_means[stat_means$value2%in%"withMUT",4]/stat_means[stat_means$value2%in%"withoutMUT",4]
  
  logFCmut_no_mut_df<-data.frame(genes=stat_means[stat_means$value2%in%"withMUT","Var2"],logFCmut_no_mut=logFCmut_no_mut)

  res_stat_with_mean <- merge(stats_df_temp,logFCmut_no_mut_df,by.x="Var2",by.y="genes")
  
  res_stat_df<-rbind(res_stat_df,data.frame(res_stat_with_mean))
  
  }
  
  }
  
  }
  
  stat.test_sig<-res_stat_df[which(res_stat_df$p<=0.05),]
  
  pdf(paste(name_comparison,".HEATMAP.metastatic_potential.MUT.pdf",sep=""))
  stat.test_sig_heat<-stat.test_sig
  stat.test_sig_heat$logFCmut_no_mut<-log(stat.test_sig_heat$logFCmut_no_mut+1,2)
  ph<-ggplot(stat.test_sig_heat, aes(primary_tissue, Var2, fill= logFCmut_no_mut)) +
    geom_tile(color = "black")+scale_fill_gradient2(low = "white", mid = "orange2", high = "blue2", midpoint = .25)+theme_bw()
  print(ph)

  dev.off()
  
  #
  # Use an heatmap to represent the met pot fold-change
  #
  
  res_stat_df$statistic<-as.numeric(res_stat_df$statistic)

  stat.test_sig<-stat.test_sig[order(stat.test_sig[,"p"],decreasing=F),]
    
  stat.test_sig$rank<-1:nrow(stat.test_sig)
  
  output_txt<-paste(name_comparison,".metastatic_potential.MUT.pval=0.05.txt",sep="")
  
  write.table(stat.test_sig,file=output_txt,sep="\t",row.names=F,quote=F)
  
  if(nrow(stat.test_sig)!=0){
    
  input_boxplot<-tab_mut_with_ann_and_metpot[which(tab_mut_with_ann_and_metpot$primary_tissue %in% stat.test_sig[,2] & tab_mut_with_ann_and_metpot$Var2 %in% stat.test_sig[,1]),]
  
  pdf(paste(name_comparison,".metastatic_potential.MUT.pdf",sep=""))
  
  stat.test_sig$fold_change_status<-ifelse(stat.test_sig$logFCmut_no_mut>1,"UP","DOWN")
  
  # ps<-ggplot(stat.test_sig, aes(x=rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
  #   geom_point(size=3) + 
  #   geom_text_repel(data=stat.test_sig, mapping=aes(x=rank, y=-log10(p), label=Var2), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")  
  # 
  # 
  # print(ps)
  # 
   ps2<-ggplot(stat.test_sig, aes(x=rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
        geom_point(size=3) + 
        geom_text_repel(data=stat.test_sig, mapping=aes(x=rank, y=-log10(p), label=Var2), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")+coord_polar()+theme_bw()
  
   print(ps2)
  
  stat.test_sig_temp<-stat.test_sig[which(stat.test_sig$fold_change_status%in%"MUT"),]
  stat.test_sig_temp<-stat.test_sig_temp[order(stat.test_sig_temp$p),]

  if(nrow(stat.test_sig_temp)>=1){
    
  stat.test_sig_temp$new_rank<-1:nrow(stat.test_sig_temp)
    
  # ps<-ggplot(stat.test_sig_temp, aes(x=new_rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
  #     geom_point(size=3) + 
  #     geom_text_repel(data=stat.test_sig_temp, mapping=aes(x=new_rank, y=-log10(p), label=Var2), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")  
  #   
  # print(ps)
    
  ps2<-ggplot(stat.test_sig_temp, aes(x=new_rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
    geom_point(size=3) +  scale_x_continuous(breaks=seq(1,max(stat.test_sig_temp$new_rank,by=1)))+
    geom_text_repel(data=stat.test_sig_temp, mapping=aes(x=new_rank, y=-log10(p), label=Var2), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")+coord_polar()+theme_bw()
  
  print(ps2)
  
  }
  
  print("inside")
  
  # stat.test_sig$tissue<-as.character(stat.test_sig$tissue)
  stat.test_sig$primary_tissue<-as.character(stat.test_sig$primary_tissue)
  
  for(tis in unique(stat.test_sig$primary_tissue)){
  
  sig_genes_tissue<-unique(stat.test_sig[stat.test_sig[,2]%in%tis,1])
  
  input_boxplot2<-input_boxplot[which(input_boxplot$Var2 %in% sig_genes_tissue),]

  sig_genes_tissue2<-names(which(rowSums(table(input_boxplot2$Var2,input_boxplot2$value2)) >= 2))
  
  if(length(unique(input_boxplot2$value2))==2){
    
  input_boxplot3<-input_boxplot2[which(input_boxplot2$primary_tissue %in%tis & input_boxplot2$Var2%in% sig_genes_tissue2),]
  colnames(input_boxplot3)[9]<-"genes"
  
  stat.test3 <- input_boxplot3 %>%
    group_by(genes) %>%
    wilcox_test(mean ~ value2)
  
  stat.test3 <- stat.test3 %>% add_xy_position(x = "genes")
  
  p<-ggplot(input_boxplot3, aes(x=genes, y=mean, color=value2)) +
    geom_boxplot()+
    geom_point(width=0.1,aes(color=value2), position = position_jitterdodge())+
    ggtitle(tis)+ylab("Metastatic potential")
  
  p2<-p+stat_pvalue_manual(stat.test3, xmin = "genes", xmax = NULL)
  print(p2)
  
  }
  
  }
  
  dev.off()
  
  
  }
  
  if(nrow(stat.test_sig)){
  
  stat.test_sig2<-data.frame(ID=name_comparison,stat.test_sig)
  stat.test_sig_all<-rbind(stat.test_sig_all,stat.test_sig2)
  
  }
  
  if(nrow(stat.test_pancancer)!=0){
    
  stat.test_sig2_pancancer<-data.frame(ID=name_comparison,res_stat_with_mean_pancancer)
  stat.test_sig_all_pancancer<-rbind(stat.test_sig_all_pancancer,stat.test_sig2_pancancer)
  
  }
  
  pancancerboxplot_tot<-rbind(pancancerboxplot_tot,pancancerboxplot)
  
}

pdf("HEATMAP.metastatic_potential.MUT.pdf")

stat.test_sig_all_heat<-stat.test_sig_all[-which(stat.test_sig_all$ID%in%c("luad","common","breast")),]

stat.test_sig_all_heat$logFCmut_no_mut<-log(stat.test_sig_all_heat$logFCmut_no_mut,2)
stat.test_sig_all_heat$p<- -log(stat.test_sig_all_heat$p,10)

ph<-ggplot(stat.test_sig_all_heat, aes(primary_tissue, Var2, fill= logFCmut_no_mut)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown2", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))

p_annotation<-ggplot(stat.test_sig_all_heat, aes(1, Var2))+geom_tile(aes(fill=ID), colour="white")+
  coord_equal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6))

print(ph+p_annotation)

dev.off()


pdf("PANCANCER_HEATMAP.metastatic_potential.MUT.pdf")

stat.test_sig_pancancer_all_heat<-stat.test_sig_all_pancancer[-which(stat.test_sig_all_pancancer$ID%in%c("luad","common","breast")),]

stat.test_sig_pancancer_all_heat$logFCmut_no_mut<-log(stat.test_sig_pancancer_all_heat$logFCmut_no_mut,2)
stat.test_sig_pancancer_all_heat$p<- -log(stat.test_sig_pancancer_all_heat$p,10)

ph<-ggplot(stat.test_sig_pancancer_all_heat, aes(ID, Var2, fill= logFCmut_no_mut)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))

p_annotation<-ggplot(stat.test_sig_pancancer_all_heat, aes(1, Var2))+geom_tile(aes(fill=ID), colour="white")+
  coord_equal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6))

print(ph+p_annotation)

#
# Plot boxplot
#

genes_for_boxplot<- unique(pancancerboxplot_tot$genes)

for(gfb in genes_for_boxplot){
  
  current_gene_pancancer<-pancancerboxplot_tot[which(pancancerboxplot_tot$genes%in%gfb),]  
  current_pvalue<-stat.test_sig_pancancer_all_heat[which(stat.test_sig_pancancer_all_heat$genes%in%gfb),"p"]
  
  p<-ggplot(current_gene_pancancer, aes(x=genes, y=mean, color=value2)) +
    geom_boxplot()+
    ggtitle(paste(gfb,"p-value:",current_pvalue))+ylab("Metastatic potential")
  
  p2<-ggplot(current_gene_pancancer, aes(x=genes, y=value, color=value2)) +
    geom_boxplot()+
    ggtitle(paste(gfb,"p-value:",current_pvalue))+ylab("Copy Number Levels")
  
  print(p+p2)
  
}

dev.off()

pdf("HEATMAP.metastatic_potential.MUT.BRCA_LUAD.pdf")

stat.test_sig_all_heat<-stat.test_sig_all[which(stat.test_sig_all$ID%in%c("luad","breast")),]

stat.test_sig_all_heat$logFCmut_no_mut<-log(stat.test_sig_all_heat$logFCmut_no_mut,2)
stat.test_sig_all_heat$p<- -log(stat.test_sig_all_heat$p,10)

p_annotation<-ggplot(stat.test_sig_all_heat, aes(1, Var2))+geom_tile(aes(fill=ID), colour="white")+
  coord_equal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6))

ph<-ggplot(stat.test_sig_all_heat, aes(primary_tissue, Var2, fill= logFCmut_no_mut)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))

print(ph+p_annotation)

dev.off()


pdf("PANCANCER_HEATMAP.metastatic_potential.MUT.BRCA_LUAD.pdf")

stat.test_sig_pancancer_all_heat<-stat.test_sig_all_pancancer[which(stat.test_sig_all_pancancer$ID%in%c("luad","breast")),]

stat.test_sig_pancancer_all_heat$logFCmut_no_mut<-log(stat.test_sig_pancancer_all_heat$logFCmut_no_mut,2)
stat.test_sig_pancancer_all_heat$p<- -log(stat.test_sig_pancancer_all_heat$p,10)

p_annotation<-ggplot(stat.test_sig_pancancer_all_heat, aes(1, Var2))+geom_tile(aes(fill=ID), colour="white")+
  coord_equal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6))

ph<-ggplot(stat.test_sig_pancancer_all_heat, aes(ID, Var2, fill= logFCmut_no_mut)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))

print(ph+p_annotation)

#
# Plot boxplot
#

genes_for_boxplot<- unique(pancancerboxplot_tot$genes)

for(gfb in genes_for_boxplot){
  
  current_gene_pancancer<-pancancerboxplot_tot[which(pancancerboxplot_tot$genes%in%gfb),]  
  current_pvalue<-stat.test_sig_pancancer_all_heat[which(stat.test_sig_pancancer_all_heat$genes%in%gfb),"p"]
  
  p<-ggplot(current_gene_pancancer, aes(x=genes, y=mean, color=value2)) +
    geom_boxplot()+
    ggtitle(paste(gfb,"p-value:",current_pvalue))+ylab("Metastatic potential")
  
  p2<-ggplot(current_gene_pancancer, aes(x=genes, y=value, color=value2)) +
    geom_boxplot()+
    ggtitle(paste(gfb,"p-value:",current_pvalue))+ylab("Copy Number Levels")
  
  print(p+p2)
  
}

dev.off()


