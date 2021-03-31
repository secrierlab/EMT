setwd("/home/guidantoniomt/pseudospace/ml_for_ppt/cosmic_focal_broad_variants")

#
# Compare lasso with random forest
#

list_rdata_lasso<-c("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData","HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData","HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

list_rdata_rf<-c("RF_HMM_nstates3_mock.mes.vs.epi.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData","RF_HMM_nstates3_mock.mes.vs.mix.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData","RF_HMM_nstates3_mock.mix.vs.epi.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData")

list_class1<-c("mes","mes","mix")
list_class2<-c("epi","mix","epi")

list_lasso_res<-vector(mode="list",2)
list_common_genes_between_experiments<-vector(mode="list",2)

for(i in 1:length(list_rdata_lasso)){

	print(i)

	class1<-list_class1[i]
	class2<-list_class2[i]

	load(list_rdata_lasso[i])

	# get top 50 genes from lasso	
	df_coef2<-do.call(rbind,res_lasso)
	
	df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

        features<-data.frame(table(df_coef2[,1]))

        lasso_features_top50<-features[order(features[,2],decreasing=T),1][1:50]
	
	# get the most frequent (80%) genes from lasso
        df_coef2<-do.call(rbind,res_lasso)

        df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

        features<-data.frame(table(df_coef2[,1]))

        features_freq_df<-features[order(features[,2],decreasing=T),]

        lasso_features_80<-features_freq_df[features_freq_df[,2]>=800,1]

	list_lasso_res[[1]]<-lasso_features_top50
	list_lasso_res[[2]]<-lasso_features_80

	names(list_lasso_res)<-c("top50","perc_80")
        library(VennDetail)

	pdf(paste(class1,"_",class2,"_compare_lasso_specific.pdf",sep=""))
        venn<-venndetail(list_lasso_res)
        plot(venn)
        p<-plot(venn, type = "upset")
        print(p)
        dev.off()

	load(list_rdata_rf[i])

	rf_features<-features_selected_rf[,1]
	
	library(VennDiagram)
	library(RColorBrewer)
	library(VennDetail)

	myCol <- brewer.pal(3, "Pastel2")

	common_genes<-Reduce(intersect,list(
                         Lasso80_perc= as.character(lasso_features_80),
                         RF=as.character(rf_features)))

	list_common_genes_between_experiments[[i]]<-common_genes

	pdf(paste(class1,"_",class2,"_compare_lasso_with_rf.pdf",sep=""))
	venn<-venndetail(list(
       		        Lasso80_perc= as.character(lasso_features_80),
            		RF=as.character(rf_features)))
	plot(venn)
	
	p<-plot(venn, type = "upset")
	print(p)
	dev.off()

}
names(list_common_genes_between_experiments)<-c("mes_vs_epi",
						"mes_vs_mix",
						"mix_vs_epi")

save(x=list_common_genes_between_experiments,
	   file="common_genomic_features_lasso_random_forest.RData")

#
# Compare lasso results
#

list_rdata_lasso<-c("HMM_nstates3_mock.mes.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData","HMM_nstates3_mock.mes.vs.mix.tissue.TRUE.1000.cosmic_arms_focal.RData",
"HMM_nstates3_mock.mix.vs.epi.tissue.TRUE.1000.cosmic_arms_focal.RData")

results_lasso_50<-vector(mode="list",3)
names(results_lasso_50)<-c("mes_vs_epi",
			   "mes_vs_mix",
			   "mix_vs_epi")

results_lasso_80<-vector(mode="list",3)
names(results_lasso_80)<-c("mes_vs_epi",
			   "mes_vs_mix",
			   "mix_vs_epi")

for(i in 1:length(list_rdata_lasso)){

        print(i)

	load(list_rdata_lasso[i])

        # get top 50 genes from lasso   
        df_coef2<-do.call(rbind,res_lasso)

        df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

        features<-data.frame(table(df_coef2[,1]))

        lasso_features_top50<-features[order(features[,2],decreasing=T),1][1:50]

	results_lasso_50[[i]]<-as.character(lasso_features_top50)

        # get the most frequent (80%) genes from lasso
        df_coef2<-do.call(rbind,res_lasso)

        df_coef2<-df_coef2[grep(df_coef2[,1],pattern=paste(c("Intercept","tumors"),collapse="|"),invert=T),]

        features<-data.frame(table(df_coef2[,1]))

        features_freq_df<-features[order(features[,2],decreasing=T),]

        lasso_features_80<-features_freq_df[features_freq_df[,2]>=800,1]

	results_lasso_80[[i]]<-as.character(lasso_features_80)

}


        library(VennDiagram)
        library(VennDetail)
        
	pdf(paste("venny_lasso_analysis_top50.pdf",sep=""))
        venn<-venndetail(results_lasso_50)
        plot(venn)
        p<-plot(venn, type = "upset")
        print(p)
        dev.off()


        pdf(paste("venny_lasso_analysis_80perc.pdf",sep=""))
        venn<-venndetail(results_lasso_80)
        plot(venn)
        p<-plot(venn, type = "upset")
        print(p)
        dev.off()


#
# Compare random forest 
#

list_rdata_rf<-c("RF_HMM_nstates3_mock.mes.vs.epi.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData","RF_HMM_nstates3_mock.mes.vs.mix.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData","RF_HMM_nstates3_mock.mix.vs.epi.tissue.FALSE.1000.topgenes.800.cosmic_arms_focal_train_0.80.RData")

results_rf<-vector(mode="list",3)
names(results_rf)<-c("mes_vs_epi",
                           "mes_vs_mix",
                           "mix_vs_epi")

list_class1<-c("mes","mes","mix")
list_class2<-c("epi","mix","epi")

for(i in 1:length(list_rdata_lasso)){

        print(i)

	load(list_rdata_rf[i])

        rf_features<-features_selected_rf[,1]

	results_rf[[i]]<-as.character(rf_features)
}

        pdf(paste("venny_rf_analysis.pdf",sep=""))
        venn<-venndetail(results_rf)
        plot(venn)
        p<-plot(venn, type = "upset")
        print(p)
        dev.off()

