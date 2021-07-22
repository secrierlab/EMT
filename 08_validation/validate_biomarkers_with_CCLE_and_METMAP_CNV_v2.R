library(data.table)
library(openxlsx)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(vcd)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

RadarTheme<-theme(panel.background=element_blank(),
                  legend.position="bottom",legend.title=element_blank(),legend.direction="vertical",
                  axis.text.x = element_blank(),
                  axis.line.x=element_line(size=0.5),
                  panel.grid.major=element_line(size=0.3,linetype = 2,colour="grey"))

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

load("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_epi<-getFeaturesSS(res_lasso,ss=800)
cnv_markers_mes_vs_epi<-gsub(grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
 cnv_markers_mes_vs_epi_raw<-grep(markers_mes_vs_epi,pattern="dndscv",value=T,invert=T)
coef_mes_vs_epi<-getCoefSS(res_lasso,c(cnv_markers_mes_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_mes_vs_epi2<-merge(coef_mes_vs_epi,mu,by="coef")
coef_mes_vs_epi2[,1]<-gsub(coef_mes_vs_epi2[,1],pattern="^X",replacement="")

cnv_markers_mes_vs_epi<-intersect(cnv_markers_mes_vs_epi,unique(coef_mes_vs_epi2[coef_mes_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_pemt_vs_epi<-getFeaturesSS(res_lasso,ss=800)
cnv_markers_pemt_vs_epi<-gsub(grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
 cnv_markers_pemt_vs_epi_raw<-grep(markers_pemt_vs_epi,pattern="dndscv",value=T,invert=T)
coef_pemt_vs_epi<-getCoefSS(res_lasso,c(cnv_markers_pemt_vs_epi_raw))
mu<-aggregate(values ~ coef, data=coef_pemt_vs_epi, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_pemt_vs_epi2<-merge(coef_pemt_vs_epi,mu,by="coef")
coef_coef_pemt_vs_epi2[,1]<-gsub(coef_coef_pemt_vs_epi2[,1],pattern="^X",replacement="")

cnv_markers_pemt_vs_epi<-intersect(cnv_markers_pemt_vs_epi,unique(coef_coef_pemt_vs_epi2[coef_coef_pemt_vs_epi2$color%in%"positive",1]))

load("HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData")

markers_mes_vs_mix<-getFeaturesSS(res_lasso,ss=800)
cnv_markers_mes_vs_mix<-gsub(grep(markers_mes_vs_mix,pattern="dndscv",value=T,invert=T),pattern="^X",replacement="")
cnv_markers_mes_vs_mix_raw<-grep(markers_mes_vs_mix,pattern="dndscv",value=T,invert=T)
coef_mes_vs_mix<-getCoefSS(res_lasso,c(cnv_markers_mes_vs_mix_raw))
mu<-aggregate(values ~ coef, data=coef_mes_vs_mix, median)
mu$color<-as.factor(ifelse(mu$values>0,"positive","negative"))
colnames(mu)[2]<-"mu"
coef_coef_mes_vs_mix2<-merge(coef_mes_vs_mix,mu,by="coef")
coef_coef_mes_vs_mix2[,1]<-gsub(coef_coef_mes_vs_mix2[,1],pattern="^X",replacement="")

cnv_markers_mes_vs_mix<-intersect(cnv_markers_mes_vs_mix,unique(coef_coef_mes_vs_mix2[coef_coef_mes_vs_mix2$color%in%"positive",1]))



raw_markers_specifics_mes_vs_epi<-setdiff(c(cnv_markers_mes_vs_epi),c(cnv_markers_pemt_vs_epi)) #there are not significant raw_markers for mes vs epi during the downstream analysis, so i don't consider mix vs mes
raw_markers_specifics_pemt_vs_epi<-setdiff(y=c(cnv_markers_mes_vs_epi,cnv_markers_mes_vs_mix),x=c(cnv_markers_pemt_vs_epi))
raw_markers_common<-intersect(c(cnv_markers_mes_vs_epi),c(cnv_markers_pemt_vs_epi))

raw_markers_specifics_mes_vs_mix<-setdiff(x=c(cnv_markers_mes_vs_mix),y=c(cnv_markers_mes_vs_epi,cnv_markers_pemt_vs_epi))

list_common_all<-list(c(cnv_markers_mes_vs_epi),c(cnv_markers_mes_vs_mix),c(cnv_markers_pemt_vs_epi))

common_between_three_groups<-Reduce(intersect,list_common_all)

#
# Upload the biomarkers
#

tab_sampleinfo<-fread("/home/data/pseudospace/depmap/primary-screen-cell-line-info.csv",data.table=F)

#
# Upload MET500
#

brain_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 1)
colnames(brain_metmap500)[1]<-"cell_line"

lung_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 2)
colnames(lung_metmap500)[1]<-"cell_line"

liver_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 3)
colnames(liver_metmap500)[1]<-"cell_line"

bone_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 4)
colnames(bone_metmap500)[1]<-"cell_line"

kidney_metmap500<-read.xlsx("/home/data/pseudospace/CCLE/MetMap500_met_potential.xlsx",sheet = 5)
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

tab_cna_ccle<-fread(file="/home/data/pseudospace/CCLE/CCLE_gene_cn.csv",header=T,fill = TRUE,data.table=F)
tab_cna_ccle_subselect<-tab_cna_ccle[tab_cna_ccle[,1]%in%codes_to_use_metmap,]
colnames(tab_cna_ccle_subselect)<-gsub(sapply(strsplit(colnames(tab_cna_ccle_subselect),split="\\("),"[[",1),pattern=" ",replacement="")
colnames(tab_cna_ccle_subselect)[1]<-"samples"

#
# See the effect of the CNA in CCLE and on the metastatic potential
#

list_features<-vector(mode="list",4)
list_features[[1]]<-raw_markers_specifics_mes_vs_epi
list_features[[2]]<-raw_markers_specifics_pemt_vs_epi
list_features[[3]]<-raw_markers_specifics_mes_vs_mix
list_features[[4]]<-raw_markers_common
#list_features[[5]]<-c("NONO_focal","MPL_focal","KIAA1459_focal","ETNK1_focal","FIP1L1_focal","CNTNAP2_focal","CREB3L2_focal")
#list_features[[6]]<-c("PRDM2_focal")
list_features[[5]]<-c("RUNX1_focal","NCOR1_focal","EPHA7_focal","TSHR_focal","CUX1_focal")
list_features[[6]]<-c("MDS2_focal","PAX7_focal","PRDM2_focal","ID3_focal","ARID1A_focal","BAP1_focal","PBRM1_focal")


# list_features[[5]]<-common_between_three_groups
# list_features[[6]]<-c("EGFR_dndscv","SETD2_dndscv","RBM10_dndscv","CDKN2A_dndscv","RB1_dndscv","MGA_dndscv","ZIC1_dndscv","NF1_dndscv","REG3A_dndscv","ARID2_dndscv","ZFP36L1_dndscv","ARID1A_dndscv")
# list_features[[7]]<-c("CBFB_dndscv","NCOR1_dndscv","AKT1_dndscv","SHISA4_dndscv","NF1_dndscv","FOXA1_dndscv")

# names(list_features)<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common","common_all","luad","breast")
names(list_features)<-c("mes_vs_epi","hemt_vs_epi","mes_vs_mix","common","luad","breast")

setwd("/home/data/pseudospace/CCLE")

stat.test_sig_all<-data.frame()
stat.test_sig_all_pancancer<-data.frame()
pancancerboxplot_tot<-data.frame()

for(i in 1:length(list_features)){
  
  print(i)
  
  name_comparison<-names(list_features)[i]
  features<-sapply(strsplit(grep(grep(list_features[[i]],pattern="dndscv",invert=T,value=T),pattern="arm",invert=T,value=T),split="\\_"),"[[",1)
  
  if(length(features)>=1){
    
  tab_cna_features<-tab_cna_ccle_subselect[,which(colnames(tab_cna_ccle_subselect) %in% c(features))]
  
  tab_cna_features2<-melt(data.frame(ID=tab_cna_ccle_subselect[,1],tab_cna_features))
  
  if(length(features)==1){
    
  tab_cna_features2[,2]<-features
  
  }
  
  tab_cna_features2[,1]<-as.character(tab_cna_features2[,1])
  
  tab_cna_features2$value2<-rep(0,nrow(tab_cna_features2))
  tab_cna_features2[which(tab_cna_features2$value >= 1),"value2"]  <- 1
  tab_cna_features2[which(tab_cna_features2$value <= -1),"value2"] <- 1
  
  tab_cna_with_ann<-merge(tab_sampleinfo[,c(2:4)],tab_cna_features2,by.x="depmap_id",by.y="ID")
  tab_cna_with_ann_and_metpot<-merge(met_pot_all2,tab_cna_with_ann,by.x="cell_line",by.y="ccle_name")
  colnames(tab_cna_with_ann_and_metpot)[9]<-"genes"
  
  tab_cna_with_ann_and_metpot$value2<-factor(ifelse(tab_cna_with_ann_and_metpot$value>=1,"withCNA","withoutCNA"),levels=c("withoutCNA","withCNA"))
  
  input_for_stat_pancancer<-tab_cna_with_ann_and_metpot 
  freq_genes_status_pancancer<- table(input_for_stat_pancancer$genes,input_for_stat_pancancer$value2)
  
  if(nrow(freq_genes_status_pancancer)>1){
    
  genes_to_consider_pancancer<- rownames(freq_genes_status_pancancer[freq_genes_status_pancancer[,1]>1 & freq_genes_status_pancancer[,2]>1,])
  
  } else {
    
  genes_to_consider_pancancer<-  rownames(freq_genes_status_pancancer)
    
  }
  
  stat.test_pancancer <-input_for_stat_pancancer[which(input_for_stat_pancancer$genes %in%genes_to_consider_pancancer),] %>%
                      group_by(genes) %>%
                      wilcox_test(mean ~ value2)
  
  stat.test_pancancer_df_temp<-data.frame(stat.test_pancancer)
  
  input_for_mean_pancancer<- input_for_stat_pancancer[which(input_for_stat_pancancer$genes %in%genes_to_consider_pancancer),]
  
  input_for_mean_pancancer$mean_unlog<-10^input_for_mean_pancancer$mean
  
  stat_means_pancancer<-aggregate(mean_unlog ~ genes + value2 , data = input_for_mean_pancancer, FUN = median, na.rm = TRUE)
  
  logFCmut_no_mut<-stat_means_pancancer[stat_means_pancancer$value2%in%"withCNA",3]/stat_means_pancancer[stat_means_pancancer$value2%in%"withoutCNA",3]
  
  logFCmut_no_mut_df<-data.frame(genes=stat_means_pancancer[stat_means_pancancer$value2%in%"withCNA","genes"],logFCmut_no_mut=logFCmut_no_mut)
  
  res_stat_with_mean_pancancer <- merge(stat.test_pancancer_df_temp,logFCmut_no_mut_df,by.x="genes",by.y="genes")
  
  stat.test_sig2_pancancer<-res_stat_with_mean_pancancer[which(res_stat_with_mean_pancancer$p<=0.05),]
  print(paste("The number of rows in the pan-cancer analysis is:",nrow(stat.test_sig2_pancancer)))
  
  
  pancancerboxplot<-input_for_stat_pancancer[which(input_for_stat_pancancer$genes%in%stat.test_sig2_pancancer[,1]),]
  
  #
  # Compare the metastatic potential
  #
  
  #
  # Compare the metastatic potential
  #
  res_stat_df<-data.frame()
  
  for(tss in unique(tab_cna_with_ann_and_metpot$primary_tissue)){
    
    print(tss)
    
    input_for_stat<-tab_cna_with_ann_and_metpot[which(tab_cna_with_ann_and_metpot$primary_tissue%in%tss),] 
    freq_genes_status<- table(input_for_stat$genes,input_for_stat$value2)

    print("inside") 
      
    if(nrow(freq_genes_status)>=1){
      
    if(nrow(freq_genes_status)>=2){
        
    genes_to_consider<- rownames(freq_genes_status[freq_genes_status[,1]>1 & freq_genes_status[,2]>1,])
    
    }else{
    
    if(length(freq_genes_status[freq_genes_status[,1]>1 & freq_genes_status[,2]>1,])==2){
      
        genes_to_consider<- rownames(freq_genes_status)
      
    } else {
       
       genes_to_consider<- NULL
  
    }
      
    }
      
    }
      
    if(length(genes_to_consider)!=0){
      

      stat.test <-input_for_stat[which(input_for_stat$genes %in% genes_to_consider),] %>%
        group_by(primary_tissue, genes) %>%
        wilcox_test(mean ~ value2)
      
      stat.test<-stat.test[which(stat.test$n1 >= 10 & stat.test$n2 >= 10),]
      
      print(nrow(stat.test))
      
      if(nrow(stat.test)>=1){
        
        print("true inside")
        
        stats_df_temp<-data.frame(stat.test)
        
        input_for_mean<-input_for_stat[which(input_for_stat$genes %in%genes_to_consider),]
        
        input_for_mean$mean_unlog<-10^input_for_mean$mean
        
        stat_means<-aggregate(mean_unlog ~ primary_tissue + genes + value2 , data = input_for_mean, FUN = median, na.rm = TRUE)
        
        # stat_means$value2<-factor(ifelse(stat_means$value2>=1,"up","down"),levels=c("up","down"))
        
        if(length(unique(stat_means$value2))==2){
        
          print("TRUE")
          
        logFCmut_no_mut<-stat_means[stat_means$value2%in%"withCNA","mean_unlog"]/stat_means[stat_means$value2%in%"withoutCNA","mean_unlog"]

        logFCmut_no_mut_df<-data.frame(genes=as.character(stat_means[stat_means$value2%in%"withCNA","genes"]),logFCmut_no_mut=logFCmut_no_mut)
        
        res_stat_with_mean <- merge(stats_df_temp,logFCmut_no_mut_df,by.x="genes",by.y="genes")
        
        res_stat_df<-rbind(res_stat_df,data.frame(res_stat_with_mean))
        }
        
      }
      
    }
    
  }
  
  if(nrow(res_stat_df)>=1){
    
  stat.test_sig<-res_stat_df[which(res_stat_df$p<=0.05),]
  
  pdf(paste(name_comparison,".HEATMAP.metastatic_potential.CNV.July.pdf",sep=""))
  stat.test_sig_heat<-stat.test_sig
  stat.test_sig_heat$logFCmut_no_mut<-log(stat.test_sig_heat$logFCmut_no_mut+1,2)
  ph<-ggplot(stat.test_sig_heat, aes(primary_tissue, genes, fill= logFCmut_no_mut)) +
    geom_tile(color = "black")+scale_fill_gradient2(low = "white", mid = "orange2", high = "blue2", midpoint = .25)+theme_bw()
  print(ph)
  
  dev.off()
  
  res_stat_df$statistic<-as.numeric(res_stat_df$statistic)
  
  stat.test_sig<-stat.test_sig[order(stat.test_sig[,"p"],decreasing=F),]
  
  stat.test_sig$rank<-1:nrow(stat.test_sig)
  
  output_txt<-paste(name_comparison,".metastatic_potential.CNV.pval=0.05txt",sep="")
  
  write.table(stat.test_sig,file=output_txt,sep="\t",row.names=F,quote=F)
  
  if(nrow(stat.test_sig)!=0){
    
  input_boxplot<-tab_cna_with_ann_and_metpot[which(tab_cna_with_ann_and_metpot$primary_tissue %in% stat.test_sig[,2] & tab_cna_with_ann_and_metpot$genes %in% stat.test_sig[,1]),]

  pdf(paste(name_comparison,".metastatic_potential.CNV.July.pdf",sep=""))
  
  stat.test_sig$fold_change_status<-ifelse(stat.test_sig$logFCmut_no_mut>1,"CNA","NO_CNA")
  
  # ps<-ggplot(stat.test_sig, aes(x=rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
  #   geom_point(size=3) + 
  #   geom_text_repel(data=stat.test_sig, mapping=aes(x=rank, y=-log10(p), label=genes), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")
  # 
  # print(ps)
  
  ps2<-ggplot(stat.test_sig, aes(x=rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
    geom_point(size=3) + 
    geom_text_repel(data=stat.test_sig, mapping=aes(x=rank, y=-log10(p), label=genes), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")+coord_polar()+theme_bw()
  
  print(ps2)
  
  stat.test_sig_temp<-stat.test_sig[which(stat.test_sig$fold_change_status%in%"CNA"),]
  stat.test_sig_temp<-stat.test_sig_temp[order(stat.test_sig_temp$p),]

  if(nrow(stat.test_sig_temp)>=1){

    stat.test_sig_temp$new_rank<-1:nrow(stat.test_sig_temp)
    
    # ps<-ggplot(stat.test_sig_temp, aes(x=new_rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
    #   geom_point(size=3) + 
    #   geom_text_repel(data=stat.test_sig_temp, mapping=aes(x=new_rank, y=-log10(p), label=genes), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")  
    # 
    # 
    # print(ps)
    
    ps2<-ggplot(stat.test_sig_temp, aes(x=new_rank, y=-log10(p), color=primary_tissue,shape=factor(fold_change_status))) + 
      geom_point(size=3) +  scale_x_continuous(breaks=seq(1,max(stat.test_sig_temp$new_rank,by=1)))+
      geom_text_repel(data=stat.test_sig_temp, mapping=aes(x=new_rank, y=-log10(p), label=genes), size=2, box.padding = 0.3, force = 0.5,segment.size = 0.1, min.segment.length = 0.05, show.legend = F, color = "black")+coord_polar()+theme_bw()
    
    print(ps2)
    
  }
  
  
  print("inside")
  
  for(tis in unique(stat.test_sig$primary_tissue)){
  
      print(tis)
        
      sig_genes_tissue<-as.character(unique(stat.test_sig[which(stat.test_sig[,2]%in%tis),1]))
      
      input_boxplot2<-input_boxplot[which(input_boxplot$genes %in% sig_genes_tissue),]
    
      sig_genes_tissue2<-names(which(rowSums(table(input_boxplot2$genes,input_boxplot2$value2)) >= 2))
      
      print(length(unique(input_boxplot2$value2)))
      
      input_boxplot3<-input_boxplot2[which(input_boxplot2$primary_tissue%in%tis & input_boxplot2$genes%in% sig_genes_tissue2),]
      # input_boxplot3$value2<-ifelse(input_boxplot3$value2==1,"withCNA","withoutCNA")
      input_boxplot3$value2<-factor(input_boxplot3$value2,levels=c("withoutCNA","withCNA"))
      
  # for(g in sig_genes_tissue2){
  #   
  #     input_cont<-input_boxplot3[which(input_boxplot3$genes%in%g),]
  #     input_cont$met_potential_discr<-rep("fill",nrow(input_cont))
  #     input_cont[which(input_cont$mean<=-4),"met_potential_discr"]<-"non_metastatic"
  #     input_cont[which(input_cont$mean>-4 & input_cont$mean < -2),"met_potential_discr"]<-"weakly_metastatic"
  #     input_cont[which(input_cont$mean>=-2),"met_potential_discr"]<-"metastatic"
  #     
  #     input_cont$value2<-factor(input_cont$value2,levels=c("withoutCNA","withCNA"))
  #       
  #     if(length(unique(input_cont$met_potential_discr))==3){
  #       
  #       contingency_table<-table(input_cont$met_potential_discr,input_cont$value2)[c("metastatic","weakly_metastatic","non_metastatic"),c("withCNA","withoutCNA")]
  #       
  #       contingency_table2<-rbind(colSums(contingency_table[c("metastatic","weakly_metastatic"),]),contingency_table["non_metastatic",])
  #       rownames(contingency_table2)<-c("metastatic","non_metastatic")
  #       
  #     }else{
  #       
  #       contingency_table<-table(input_cont$met_potential_discr,input_cont$value2)[unique(input_cont$met_potential_discr),c("withCNA","withoutCNA")]
  #       contingency_table2<-contingency_table
  #       
  #     }
  #     
  #   
  #     print(assoc(contingency_table, shade = TRUE, las=3,main=paste("not collapsed met.pot:",g)))
  #     
  #     # print(assoc(contingency_table2, shade = TRUE, las=3,main=paste("collapsed met.pot:",g)))
  # 
  # }
  
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
      
      stat.test4 <- input_boxplot3 %>%
        group_by(genes) %>%
        wilcox_test(value ~ value2)
      
      stat.test4 <- stat.test4 %>% add_xy_position(x = "genes")
      
      p3<-ggplot(input_boxplot3, aes(x=genes, y=value, color=value2)) +
        geom_boxplot()+
        geom_point(width=0.1,aes(color=value2), position = position_jitterdodge())+
        ggtitle(tis)+ylab("Copy Number Levels")
      
      p4<-p3+stat_pvalue_manual(stat.test4, xmin = "genes", xmax = NULL)
      print(p4)
  
  }
  
  dev.off()
  
  }
  
  }

  }
  
  
  

  if(nrow(stat.test_sig)>=1){
    
    stat.test_sig2<-data.frame(ID=name_comparison,stat.test_sig)
    stat.test_sig_all<-rbind(stat.test_sig_all,stat.test_sig2)
    
  }
  
  if(nrow(stat.test_sig2_pancancer)>=1){
    
    stat.test_sig_all_pancancer<-rbind(stat.test_sig_all_pancancer,data.frame(ID=name_comparison,stat.test_sig2_pancancer))
    
  }
  
  pancancerboxplot_tot<-rbind(pancancerboxplot_tot,pancancerboxplot)
  
}

print("end analysis")

pdf("HEATMAP.metastatic_potential.CNV.July.pdf")

stat.test_sig_all_heat<-stat.test_sig_all[-which(stat.test_sig_all$ID%in%c("luad","common","breast")),]

stat.test_sig_all_heat$logFCmut_no_mut<-log(stat.test_sig_all_heat$logFCmut_no_mut,2)
stat.test_sig_all_heat$p<- -log(stat.test_sig_all_heat$p,10)

ph<-ggplot(stat.test_sig_all_heat, aes(primary_tissue, genes, fill= logFCmut_no_mut)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "brown", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))

p_annotation<-ggplot(stat.test_sig_all_heat, aes(1, genes))+geom_tile(aes(fill=ID), colour="white")+
  coord_equal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6))

print(ph+p_annotation)

dev.off()

pdf("PANCANCER_HEATMAP.metastatic_potential.CNV.July.pdf")

if(length(which(stat.test_sig_all_pancancer$ID%in%c("luad","common","breast")))==0){
  
stat.test_sig_pancancer_all_heat<-stat.test_sig_all_pancancer

} else {

stat.test_sig_pancancer_all_heat<-stat.test_sig_all_pancancer[-which(stat.test_sig_all_pancancer$ID%in%c("luad","common","breast")),]
  
}

stat.test_sig_pancancer_all_heat$logFCmut_no_mut<-log(stat.test_sig_pancancer_all_heat$logFCmut_no_mut,2)
stat.test_sig_pancancer_all_heat$p<- -log(stat.test_sig_pancancer_all_heat$p,10)

ph<-ggplot(stat.test_sig_pancancer_all_heat, aes(ID, genes, fill= logFCmut_no_mut)) +
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue2", mid = "white", high = "brown", midpoint = 0)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))

p_annotation<-ggplot(stat.test_sig_pancancer_all_heat, aes(1, genes))+geom_tile(aes(fill=ID), colour="white")+
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



# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# There are not significant genes in BRCA AND LUAD  June
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pdf("HEATMAP.metastatic_potential.CNV.BRCA_LUAD.July.pdf")
 
 stat.test_sig_all_heat<-stat.test_sig_all[which(stat.test_sig_all$ID%in%c("luad","breast")),]

 stat.test_sig_all_heat$logFCmut_no_mut<-log(stat.test_sig_all_heat$logFCmut_no_mut,2)
 stat.test_sig_all_heat$p<- -log(stat.test_sig_all_heat$p,10)
 
 p_annotation<-ggplot(stat.test_sig_all_heat, aes(1, genes))+geom_tile(aes(fill=ID), colour="white")+
   coord_equal()+
   theme(axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.text.x = element_text(angle = 90, vjust = 0.6))
 
 ph<-ggplot(stat.test_sig_all_heat, aes(primary_tissue, genes, fill= logFCmut_no_mut)) +
   geom_tile(color = "black")+
   scale_fill_gradient2(low = "blue2", mid = "white", high = "brown", midpoint = 0)+theme_bw()+
   theme(axis.text.x = element_text(angle = 90))+
   geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))
 
 print(ph+p_annotation)
 
dev.off()

# 
# 

 pdf("PANCANCER_HEATMAP.metastatic_potential.CNV.BRCA_LUAD.July.pdf")
 
 stat.test_sig_pancancer_all_heat<-stat.test_sig_all_pancancer[which(stat.test_sig_all_pancancer$ID%in%c("luad","breast")),]

 stat.test_sig_pancancer_all_heat$logFCmut_no_mut<-log(stat.test_sig_pancancer_all_heat$logFCmut_no_mut,2)
 stat.test_sig_pancancer_all_heat$p<- -log(stat.test_sig_pancancer_all_heat$p,10)
 
 p_annotation<-ggplot(stat.test_sig_pancancer_all_heat, aes(1, genes))+geom_tile(aes(fill=ID), colour="white")+
   coord_equal()+
   theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.text.x = element_text(angle = 90, vjust = 0.6))
 
 ph<-ggplot(stat.test_sig_pancancer_all_heat, aes(ID, genes, fill= logFCmut_no_mut)) +
   geom_tile(color = "black")+
  scale_fill_gradient2(low = "brown", mid = "white", high = "blue2", midpoint = 0)+theme_bw()+
   theme(axis.text.x = element_text(angle = 90))+
   geom_point(aes(colour = p,size=p))+scale_colour_gradientn(colours = ((brewer.pal(11,"Greys"))))
 
 print(ph+p_annotation)
 
dev.off()
