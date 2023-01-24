findComponents<-function(ps,output,xstring,kappa=3,manual=FALSE,cutoff_man=NULL){

        mixmdl = normalmixEM(ps$pseudospace,k=kappa,fast=TRUE,maxit=500000)

        p<-ggplot(ps, aes(x = pseudospace, fill = ..density..)) + geom_density(fill = "#E69F00")+xlab(xstring)+geom_vline(xintercept=mixmdl$mu,color='red')+ylab("Nearest Neighbor density across Pseudospace") +theme(text = element_text(size = 15),axis.text.y = element_text(size = 10)) + monocle:::monocle_theme_opts()

        pdf(paste(output,'.mixmdl.pdf',sep=''))
        print(p)
        plot(mixmdl,density=T,whichplots = 3)
        dev.off()

        print(mixmdl$mu)

        cutoff<-max(mixmdl$mu)
        print(cutoff)
	
	if(manual==T){
	cutoff<-cutoff_man
        ps$groups<-ifelse(ps$pseudospace>=cutoff,1,0)

	}

	ps$groups<-ifelse(ps$pseudospace>=cutoff,1,0)

        return(ps)
}

findDEGsgenes<-function(mtx_exp,ps,qvalue,fold_change,output,gs_presence=TRUE,gs=NULL,assay='TCGA'){

		if(assay=='TCGA'){
		ps<-ps[ps[,2]%in%colnames(mtx_exp),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
		mtx_exp2<-mtx_exp[,c(ps[,2]%in%colnames(mtx_exp))]
	        mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,2],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'


		} else{

		ps<-ps[ps[,1]%in%colnames(mtx_exp),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps[,1]%in%colnames(mtx_exp))]
		mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,1],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'
	
		}

		
		if(gs_presence==TRUE){
			
		mtx_exp3<-mtx_exp3[mtx_exp3[,1]%in%gs,]		
		print(dim(mtx_exp3))

		}

                # collapse the expression of the genes
                # mtx_exp4<-aggregate(.~genesID,mtx_exp3,mean)
		mtx_exp4<-as.data.frame(setDT(mtx_exp3)[, lapply(.SD, mean), by = .(genesID)])
                rownames(mtx_exp4)<-mtx_exp4[,1]

                mtx_for_degs<-data.frame(groups=as.factor(ps$groups),t(mtx_exp4[,-1]))

                print('find degs')

                baseformula<-"~groups"

                resPVAL<-NULL
                resFC<-NULL
                resTRENDPVAL<-NULL
                resTRENDFC<-NULL

                for(i in 2:ncol(mtx_for_degs)){
                print(i)
                # https://blog.minitab.com/blog/adventures-in-statistics-2/choosing-between-a-nonparametric-test-and-a-parametric-test

                formula<-paste(colnames(mtx_for_degs)[i],baseformula,sep='')
                # p <- summary(aov(as.formula(formula), data=mtx_for_degs))[[1]][["Pr(>F)"]][1]
		p<-wilcox.test(as.formula(formula), data=mtx_for_degs)$p.value

                fc<-aggregate(as.formula(formula),mtx_for_degs,mean)
                #High minus low metastatic
                fc2<-fc[2,2]-fc[1,2]

                # p <- coef(summary(glm(as.formula(formula),data=mtx_for_degs,family=binomial('logit'))))[2,4]          

                resPVAL<-c(resPVAL,p)
                resFC<-c(resFC,fc2)


                # test for trend 
                current_gene<-data.frame(mtx_for_degs[,i])
                current_gene$patients<-rownames(mtx_for_degs)

                psy<-ps[order(ps$pseudospace),]

                current_gene3<-current_gene[match(psy$patients,current_gene$patients),]

                rm<-smooth.spline(current_gene3[!is.na(current_gene3[,1]),1])$y

                restrend<-cox.stuart.test(rm)$p.value

                fold_change_trend<-max(rm)-min(rm)


                resTRENDPVAL<-c(resTRENDPVAL,restrend)
                resTRENDFC<-c(resTRENDFC,fold_change_trend)

                }

                print('end degs')

                names(resPVAL)<-colnames(mtx_for_degs)[-1]
                names(resTRENDPVAL)<-colnames(mtx_for_degs)[-1]

                statistics<-data.frame(genes=colnames(mtx_for_degs)[-1],pvalue=resPVAL,pvalue_corrected=p.adjust(resPVAL,'BH'),fold_change=resFC,cox_stuart_test=resTRENDPVAL,cox_stuart_test_corrected=p.adjust(resTRENDPVAL,'BH'),fc_trend=resTRENDFC)

                statistics2<-statistics[statistics$pvalue_corrected<=qvalue,]

                statistics3<-statistics2[which(statistics2$fold_change>=fold_change | statistics2$fold_change<=-fold_change),]

                write.table(statistics3,file=paste(output,'.txt',sep=''),sep='\t',row.names=F,quote=F)

                statistic_trend_sign<-statistics[statistics$pvalue_corrected<=qvalue & statistics$fc_trend>=fold_change,]

                write.table(statistic_trend_sign,file=paste(output,'.coxstuart.txt',sep=''),sep='\t',row.names=F,quote=F)

                list_results<-vector(mode='list',length=3)

                list_results[[1]]<-statistics
                list_results[[2]]<-statistics3
                list_results[[3]]<-statistic_trend_sign

                return(list_results)
}





findDEGsgenesForCancer<-function(mtx_exp,ps,qvalue,fold_change,output,gs_presence=TRUE,gs=NULL,assay='TCGA'){

                if(assay=='TCGA'){

                ps<-ps[ps[,2]%in%colnames(mtx_exp),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps[,2]%in%colnames(mtx_exp))]
                mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,2],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'

                } else{

                ps<-ps[ps[,1]%in%colnames(mtx_exp),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps[,1]%in%colnames(mtx_exp))]
                mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,1],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'

                }


                if(gs_presence==TRUE){

                mtx_exp3<-mtx_exp3[mtx_exp3[,1]%in%gs,]
                print(dim(mtx_exp3))

                }

                # collapse the expression of the genes
                # mtx_exp4<-aggregate(.~genesID,mtx_exp3,mean)
                mtx_exp4<-as.data.frame(setDT(mtx_exp3)[, lapply(.SD, mean), by = .(genesID)])
                rownames(mtx_exp4)<-mtx_exp4[,1]

                mtx_for_degs_init<-data.frame(groups=as.factor(ps$groups),t(mtx_exp4[,-1]))


		list_results<-data.frame()

		for(tum in unique(ps[,1])){
		
		print(tum)

		samples_for_current_tumors<-ps[ps[,1]%in%tum,'patients']

		mtx_for_degs<-mtx_for_degs_init[which(rownames(mtx_for_degs_init)%in%samples_for_current_tumors),]


		if(length(unique(mtx_for_degs[,1]))>1){

                print('find degs')

                baseformula<-"~groups"

                resPVAL<-NULL
                resFC<-NULL
                resTRENDPVAL<-NULL
                resTRENDFC<-NULL

                for(i in 2:ncol(mtx_for_degs)){

                # https://blog.minitab.com/blog/adventures-in-statistics-2/choosing-between-a-nonparametric-test-and-a-parametric-test

                formula<-paste(colnames(mtx_for_degs)[i],baseformula,sep='')

                p<-wilcox.test(as.formula(formula), data=mtx_for_degs)$p.value

                fc<-aggregate(as.formula(formula),mtx_for_degs,mean)
                #High minus low metastatic
                fc2<-fc[2,2]-fc[1,2]

                # p <- coef(summary(glm(as.formula(formula),data=mtx_for_degs,family=binomial('logit'))))[2,4]          

                resPVAL<-c(resPVAL,p)
                resFC<-c(resFC,fc2)

                # test for trend 
                current_gene<-data.frame(mtx_for_degs[,i])
                current_gene$patients<-rownames(mtx_for_degs)

                psy<-ps[order(ps$pseudospace),]

                current_gene3<-current_gene[match(psy$patients,current_gene$patients),]

                rm<-smooth.spline(current_gene3[!is.na(current_gene3[,1]),1])$y

                restrend<-cox.stuart.test(rm)$p.value

                fold_change_trend<-max(rm)-min(rm)


                resTRENDPVAL<-c(resTRENDPVAL,restrend)
                resTRENDFC<-c(resTRENDFC,fold_change_trend)

                }

                print('end degs')

                names(resPVAL)<-colnames(mtx_for_degs)[-1]
                names(resTRENDPVAL)<-colnames(mtx_for_degs)[-1]

                statistics<-data.frame(tumor=tum,genes=colnames(mtx_for_degs)[-1],pvalue=resPVAL,pvalue_corrected=p.adjust(resPVAL,'BH'),fold_change=resFC,cox_stuart_test=resTRENDPVAL,cox_stuart_test_corrected=p.adjust(resTRENDPVAL,'BH'),fc_trend=resTRENDFC)

                statistics2<-statistics[statistics$pvalue_corrected<=qvalue,]

                statistics3<-statistics2[which(statistics2$fold_change>=fold_change | statistics2$fold_change<=-fold_change),]

		list_results<-rbind(list_results,statistics3)

		}

		}

                return(list_results)
}


plotGenesWithPseudospace<-function(mtx_exp,mtx_cnv,gs,ps,tab_clinical,clinical_feature,clusters,trim=0.1,trim_status=FALSE,output){

                mtx_exp2<-data.frame(mtx_exp[,1],mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])])
                colnames(mtx_exp2)[1]<-'genes'

                mtx_exp2_select<-melt(mtx_exp2[as.character(mtx_exp2[,1])%in%gs,])
		
                subjects<-sapply(strsplit(as.character(mtx_exp2_select[,2]),'\\.'),'[[',4)
                tumors<-sapply(strsplit(as.character(mtx_exp2_select[,2]),'\\.'),'[[',1)

                mtx_exp2_select$tumors<-tumors
                colnames(mtx_exp2_select)[2]<-'patients'

                mtx_exp2_select3<-merge(mtx_exp2_select,ps[,c(2:3)],by='patients')
                colnames(mtx_exp2_select3)[3]<-'expression'

                mtx_exp2_select3$patients<-gsub(make.names(subjects,unique=F),pattern='\\.',replacement='_')

                mtx_exp2_select4<-merge(mtx_exp2_select3,tab_clinical[,c('samples',clinical_feature)],by.x='patients',by.y='samples',all.x=T)

		resKmean<-data.frame()

		for(i in unique(mtx_exp2_select4$tumors)){

		input_for_kmeans<-mtx_exp2_select4[mtx_exp2_select4$tumors==i,c('patients','genes','expression')]

		input_for_kmeans_wide<-reshape(input_for_kmeans,idvar='genes',timevar='patients',direction='wide')
		rownames(input_for_kmeans_wide)<-input_for_kmeans_wide[,1]

		kmdf<-kmeans(input_for_kmeans_wide[,-1],clusters)

                kmdf<-data.frame(tumors=i,clusters=kmdf$cluster,genes=input_for_kmeans_wide[,1])

		resKmean<-rbind(resKmean,kmdf)

		}

		mtx_exp2_select5<-merge(mtx_exp2_select4,resKmean,by.x=c('genes','tumors'),by.y=c('genes','tumors'),all.x=T)

		pdf(paste(output,'.pdf',sep=''))

                for(cn in unique(mtx_exp2_select4$tumors)) {

		print(cn)

                input_for_plot<-mtx_exp2_select5[mtx_exp2_select5$tumors==cn,]

		if(trim_status==TRUE){
                input_for_plot<-input_for_plot[which(input_for_plot$expression>trim),]
		} else {
                input_for_plot[,clinical_feature]<-as.factor(input_for_plot[,clinical_feature])
		}

                q <- ggplot(aes(pseudospace, expression), data = input_for_plot)

                q2 <- q + geom_point(aes(pseudospace,expression,color=tumor_stage),data=input_for_plot)+geom_smooth(aes(color = NULL),method = "loess")+ggtitle(cn)

		print(q2)	
		
		}

		dev.off()
	
		write.table(resKmean,file=paste(output,'.kmeans.txt',sep=''),sep='\t',row.names=F,quote=F)
		
		pdf(paste(output,'.kmeans.pdf',sep=''))

                for(cn in unique(mtx_exp2_select5$tumors)) {

                print(cn)

                input_for_plot<-mtx_exp2_select5[mtx_exp2_select5$tumors==cn,]
		
		if(trim_status==TRUE){
                input_for_plot<-input_for_plot[which(input_for_plot$expression>trim),]
		}else{
                input_for_plot[,clinical_feature]<-as.factor(input_for_plot[,clinical_feature])
		}

                q <- ggplot(aes(pseudospace, expression), data = input_for_plot)

                q2 <- q + geom_point(aes(pseudospace,expression,color=tumor_stage),data=input_for_plot)+geom_smooth(aes(color = NULL),method = "loess")+ggtitle(cn) + facet_wrap(~clusters)
                print(q2)

                }

                dev.off()	

}

createMatrixFromSegFile<-function(dictionary_copy_number,annotation_file,current_genes,ps){
	
	input_annotation<-annotation_file[annotation_file[,5]%in%current_genes,]
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)

        colnames(dictionary_copy_number)[c(1:3)]<-c('chrom','start','end')
        RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))

        qh<-res_overlaps$queryHits
        sh<-res_overlaps$subjectHits

        res_query<-dictionary_copy_number[qh,]
        res_subject<-input_annotation[sh,]

        RES_ANNOTATION_CURRENT_CANCER<-data.frame(res_query,res_subject)

	dfcan_signal<-RES_ANNOTATION_CURRENT_CANCER[,c('samples','GENE_NAME','Segment_Mean')]
        dfcan_signal[,1]<-sapply(strsplit(dfcan_signal[,1],split='\\-'),'[',3)

	dfcan2<-as.data.frame(setDT(dfcan_signal)[, lapply(.SD, mean), by = .(samples, GENE_NAME)])
	
	dfcan3<-dcast(dfcan2, GENE_NAME ~ samples)

	return(dfcan3)
}


createMatrixFromSegFileForEnrichment<-function(dictionary_copy_number,annotation_file,current_genes,ps,qvalue){

	print('Pre-processing, this step is slow')
	
	dictionary_copy_number2<-dictionary_copy_number[,c(1,2,3,27,24)]

	annotation_file<-annotation_file[annotation_file[,5]%in%current_genes,]
	current_genes<-intersect(current_genes,annotation_file[,5])

	input_annotation<-annotation_file
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)
	
	colnames(dictionary_copy_number2)[c(1:3)]<-c('chrom','start','end')
 
	RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number2)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))

        qh<-res_overlaps$queryHits
        sh<-res_overlaps$subjectHits

	dictionary_copy_number2<-as.data.table(dictionary_copy_number2)
	input_annotation<-as.data.table(input_annotation)

        res_query<-dictionary_copy_number2[qh,]
        res_subject<-input_annotation[sh,]

        RES_ANNOTATION_CURRENT_CANCER<-as.data.table(data.frame(res_query,res_subject)[,c('samples','GENE_NAME','Type_CNVs_canonical')])
	
	# Crucial, order the data.frame by genes. 
	RES_ANNOTATION_CURRENT_CANCER2<-RES_ANNOTATION_CURRENT_CANCER[order(RES_ANNOTATION_CURRENT_CANCER$GENE_NAME)]

	# Get a matrix in which for each gene I known the position in the data.frame start and end
	print('I am creting a matrix GENE: START - END')

	rl<-rle(RES_ANNOTATION_CURRENT_CANCER2$GENE_NAME)
	lst<-lapply(split(seq_along(RES_ANNOTATION_CURRENT_CANCER2$GENE_NAME),rep(seq_along(rl$values),rl$lengths)),range)
	names(lst)<-rl$values
	
	rle_df<-do.call(rbind,lst)	

	print('end processing')

	print('start analysis')

	RESPVAL<-data.frame()
        RH<-NULL
        RL<-NULL

        RHAN<-NULL
        RLAN<-NULL

        RHDN<-NULL
        RLDN<-NULL

	
	for(i in 1:nrow(rle_df)){

	current_gene<-rownames(rle_df)[i]
	start<-rle_df[i,1]
	end<-rle_df[i,2]
		
	dfcan_signal<-as.data.frame(RES_ANNOTATION_CURRENT_CANCER2[start:end,])
	print(dim(dfcan_signal))
	
        dfcan_signal[,1]<-sapply(strsplit(dfcan_signal[,1],split='\\-'),'[',3)
	colnames(dfcan_signal)[1]<-'patients'	

	# The number of samples in the two datasets could be different because it is different the number of patients with CNV and GE
	dfcan_signal2<-merge(dfcan_signal,ps,'patients')
	
	dfcan_signal3<-dfcan_signal2
	dfcan_signal3$Type_CNVs_canonical<-as.character(dfcan_signal3$Type_CNVs_canonical)

	patients_with_low_emt<-as.numeric(table(dfcan_signal3[which(dfcan_signal3$groups==0),'Type_CNVs_canonical']))
	patients_with_high_emt<-as.numeric(table(dfcan_signal3[which(dfcan_signal3$groups==1),'Type_CNVs_canonical']))

	mtxstats<-matrix(c(patients_with_high_emt,patients_with_low_emt),ncol=2)
	colnames(mtxstats)<-c('high','low')
	rownames(mtxstats)<-c('depletion','amplification','neutral')

	p<-PT<-pairwiseNominalIndependence(mtxstats)
	p_matrix<-data.frame(t(data.frame(PT$p.Fisher)))
	colnames(p_matrix)<-gsub(PT$Comparison,pattern='\\ : ',replacement='_')

	#look only for amplification and depletion events 	
	mtxstats2<-mtxstats[c('amplification','depletion'),]
        ratios<-mtxstats2[1,]/mtxstats2[2,]
        
        ratios_high<-ratios[1]
        ratios_low<-ratios[2]
	
	#look only for amplification and neutral events       
        mtxstats2<-mtxstats[c('amplification','neutral'),]
        ratios<-mtxstats2[1,]/mtxstats2[2,]

	ratios_high_amp_neutral<-ratios[1]
        ratios_low_amp_neutral<-ratios[2]

	#look only for depletion and neutral events       
        mtxstats2<-mtxstats[c('depletion','neutral'),]
        ratios<-mtxstats2[1,]/mtxstats2[2,]

        ratios_high_depletion_neutral<-ratios[1]
        ratios_low_depletion_neutral<-ratios[2]

	RESPVAL<-rbind(RESPVAL,p_matrix)

	RH<-c(RH,ratios_high)
	RL<-c(RL,ratios_low)

	#HIGH: amplified and neutral ratios	
	RHAN<-c(RHAN,ratios_high_amp_neutral)
	#LOW: amplified and neutral ratios     
	RLAN<-c(RLAN,ratios_low_amp_neutral)

	#HIGH: depleted and neutral ratios     
	RHDN<-c(RLDN,ratios_high_depletion_neutral)
	#LOW: depleted and neutral ratios     
	RLDN<-c(RLDN,ratios_low_depletion_neutral)

	}
	
	RH2<-ifelse(RH>1,'amplification','depletion')
        RL2<-ifelse(RL>1,'amplification','depletion')
	

	RHAN2<-ifelse(RHAN>1,'amplification','neutral')
	RLAN2<-ifelse(RLAN>1,'amplification','neutral')

	RHDN2<-ifelse(RHDN>1,'depletion','neutral')
	RLDN2<-ifelse(RLDN>1,'depletion','neutral')

	QPVAL<-apply(RESPVAL,2,p.adjust)	
	
	statistics<-data.frame(genesID=rownames(rle_df),
	p_value=RESPVAL,
	pvalue_corrected=QPVAL,

	class_high_amp_dep=RH,
	class_low_amp_dep=RL,
	class_high_amp_neu=RHAN,
	class_high_dep_neu=RHDN,
	class_low_amp_neu=RLAN,
	class_low_dep_neu=RLDN,

	prevalence_high_amp_dep=RH2,
	prevalence_low_amp_dep=RL2,
	prevalence_high_amp_neu=RHAN2,
	prevalence_high_dep_neu=RHDN2,
	prevalence_low_amp_neu=RLAN2,
	prevalence_low_dep_new=RLDN2

	)
	
	# Select interesting biology events
	# 1) Depletion and Amplified regions in high 
	
	stats_amp_dep<-data.frame(statistics[statistics$pvalue_corrected.depletion_amplification<=qvalue,],type='Amplification_depletion')
	stats_depletion_neutral<-data.frame(statistics[statistics$pvalue_corrected.depletion_neutral<=qvalue,],type='Depletion_neutral')
	stats_amplification_neutral<-data.frame(statistics[statistics$pvalue_corrected.amplification_neutral<=qvalue,],type='Amplification_neutral')

	statistics_significant<-rbind(stats_amp_dep,stats_depletion_neutral,stats_amplification_neutral)

	list_results_stats<-vector(mode='list',2)
	
	list_results_stats[[1]]<-statistics
	list_results_stats[[2]]<-statistics_significant

        return(list_results_stats)
}


runGaia<-function(dictionary_copry_number,annotation_file,current_genes,ps){

        dictionary_copy_number2<-dictionary_copy_number[,c(1,2,3,27,24,5)]

        annotation_file<-annotation_file[annotation_file[,5]%in%current_genes,]
        current_genes<-intersect(current_genes,annotation_file[,5])

        input_annotation<-annotation_file
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)

        colnames(dictionary_copy_number2)[c(1:3)]<-c('chrom','start','end')

        RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number2)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))

        qh<-res_overlaps$queryHits
        sh<-res_overlaps$subjectHits

        dictionary_copy_number2<-as.data.table(dictionary_copy_number2)
        input_annotation<-as.data.table(input_annotation)

        res_query<-dictionary_copy_number2[qh,]
        res_subject<-input_annotation[sh,]

        RES_ANNOTATION_CURRENT_CANCER<-as.data.table(data.frame(res_query,res_subject)[,c('samples','GENE_NAME','Type_CNVs_canonical','Segment_Mean')])



}



createMatrix_CNV_tumor<-function(dictionary_copy_number,annotation_file,current_genes){

        print('Pre-processing, this step is slow')

        dictionary_copy_number2<-dictionary_copy_number[,c(1,2,3,27,24,5)]

        annotation_file<-annotation_file[annotation_file[,5]%in%current_genes,]
        current_genes<-intersect(current_genes,annotation_file[,5])

        input_annotation<-annotation_file
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)

        colnames(dictionary_copy_number2)[c(1:3)]<-c('chrom','start','end')

        RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number2)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))

        res_overlaps<-findOverlaps(query=annotation_GR,subject=RES_ALL_CANCERS_GR,type='within')

        qh<-queryHits(res_overlaps)
        sh<-subjectHits(res_overlaps)

        dictionary_copy_number2<-as.data.table(dictionary_copy_number2)
        input_annotation<-as.data.table(input_annotation)

        res_subject<-dictionary_copy_number2[sh,]
        res_query<-input_annotation[qh,]

        RES_ANNOTATION_CURRENT_CANCER<-as.data.table(data.frame(res_subject,res_query)[,c('samples','GENE_NAME','Type_CNVs_canonical','Segment_Mean')])

        # Crucial, order the data.frame by genes. 
        RES_ANNOTATION_CURRENT_CANCER2<-RES_ANNOTATION_CURRENT_CANCER[order(RES_ANNOTATION_CURRENT_CANCER$GENE_NAME)]

	# Crucial, order the data.frame by genes. 
        RES_ANNOTATION_CURRENT_CANCER2<-RES_ANNOTATION_CURRENT_CANCER[order(RES_ANNOTATION_CURRENT_CANCER$GENE_NAME)]

	return(RES_ANNOTATION_CURRENT_CANCER2)	
}


createMatrixFromSegFileForEnrichment_one_vs_rest<-function(dictionary_copy_number,annotation_file,current_genes,ps,already_matrix=T){
	
	if(already_matrix!=T){
	
        print('Pre-processing, this step is slow')

        dictionary_copy_number2<-dictionary_copy_number[,c(1,2,3,27,24,5)]

        annotation_file<-annotation_file[annotation_file[,5]%in%current_genes,]
        current_genes<-intersect(current_genes,annotation_file[,5])

        input_annotation<-annotation_file
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)

        colnames(dictionary_copy_number2)[c(1:3)]<-c('chrom','start','end')

        RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number2)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))
	  
	res_overlaps<-findOverlaps(query=annotation_GR,subject=RES_ALL_CANCERS_GR,type='within')

        qh<-queryHits(res_overlaps)
        sh<-subjectHits(res_overlaps)

        dictionary_copy_number2<-as.data.table(dictionary_copy_number2)
        input_annotation<-as.data.table(input_annotation)

        res_subject<-dictionary_copy_number2[sh,]
        res_query<-input_annotation[qh,]

	print('Slow process...wait...')

        RES_ANNOTATION_CURRENT_CANCER<-as.data.table(data.frame(res_query,res_subject)[,c('samples','GENE_NAME','Type_CNVs_canonical','Segment_Mean')])

        # Crucial, order the data.frame by genes. 
        RES_ANNOTATION_CURRENT_CANCER2<-RES_ANNOTATION_CURRENT_CANCER[order(RES_ANNOTATION_CURRENT_CANCER$GENE_NAME)]

        # Get a matrix in which for each gene I known the position in the data.frame start and end
        print('I am creting a matrix GENE: START - END indices, this is step is slow put it is to speed up the rest of the analysis')
        rl<-rle(RES_ANNOTATION_CURRENT_CANCER2$GENE_NAME)
        lst<-lapply(split(seq_along(RES_ANNOTATION_CURRENT_CANCER2$GENE_NAME),rep(seq_along(rl$values),rl$lengths)),range)
        names(lst)<-rl$values

        rle_df<-do.call(rbind,lst)

	} else {

	RES_ANNOTATION_CURRENT_CANCER = dictionary_copy_number
		
	RES_ANNOTATION_CURRENT_CANCER2<-RES_ANNOTATION_CURRENT_CANCER[order(RES_ANNOTATION_CURRENT_CANCER$GENE_NAME)]

	rl<-rle(RES_ANNOTATION_CURRENT_CANCER2$GENE_NAME)
        lst<-lapply(split(seq_along(RES_ANNOTATION_CURRENT_CANCER2$GENE_NAME),rep(seq_along(rl$values),rl$lengths)),range)
        names(lst)<-rl$values

        rle_df<-do.call(rbind,lst)
	

	} 

       	print('end processing')

        print('start analysis')

        RESPVAL_AMP<-NULL
	RESPVAL_DEL<-NULL
	RESOR_AMP<-NULL
	RESOR_DEL<-NULL

	GENES_ANALYZED<-NULL

        for(i in 1:nrow(rle_df)){
	
	print(i)

        current_gene<-rownames(rle_df)[i]

        dfcan_signal<-as.data.frame(RES_ANNOTATION_CURRENT_CANCER2[rle_df[i,1]:rle_df[i,2],])
	change_value<-sapply(strsplit(dfcan_signal[,1],split='\\-'),'[',3)
	dfcan_signal<-data.frame(change_value,dfcan_signal[,-1],stringsAsFactors=F)
        colnames(dfcan_signal)[1]<-'patients'
	
	dfcan_signal<-data.table(dfcan_signal,key='patients')

	ps<-data.table(ps,key='patients')
	
        # The number of samples in the two datasets could be different because it is different the number of patients with CNV and GE
        dfcan_signal2<-merge(dfcan_signal,ps,by='patients')

	# It the gene has only one row, does not consider for the analysis
	if(nrow(dfcan_signal2)!=1 & length(unique(dfcan_signal2$Type_CNVs_canonical))!=1){

	GENES_ANALYZED<-c(GENES_ANALYZED,current_gene)

        dfcan_signal3<-unique(dfcan_signal2[,c(1,2,3,7,9)])

        dfcan_signal3$Type_CNVs_canonical<-as.character(dfcan_signal3$Type_CNVs_canonical)

	tbl_low<-table(dfcan_signal3[which(dfcan_signal3$groups==0),'Type_CNVs_canonical'])
	tbl_high<-table(dfcan_signal3[which(dfcan_signal3$groups==1),'Type_CNVs_canonical'])

        patients_with_low_emt<-as.numeric(table(dfcan_signal3[which(dfcan_signal3$groups==0),'Type_CNVs_canonical']))
        patients_with_high_emt<-as.numeric(table(dfcan_signal3[which(dfcan_signal3$groups==1),'Type_CNVs_canonical']))
	
	length_high<-length(patients_with_high_emt)
	length_low<-length(patients_with_low_emt)

	names_high<-names(table(dfcan_signal3[which(dfcan_signal3$groups==1),'Type_CNVs_canonical']))
	names_low<-names(table(dfcan_signal3[which(dfcan_signal3$groups==0),'Type_CNVs_canonical']))

	if(length_high!=length_low){

	if(length_high<length_low){	
	
	template<-c(0,0,0)
	names(template)<-c('depletion','gain','LOH_Neutral')

	data_available<-intersect(names_high,names_low)

	template[data_available]<-tbl_high[data_available]

	patients_with_high_emt<-template
	
	} else {

		
	template<-c(0,0,0)
        names(template)<-c('depletion','gain','LOH_Neutral')

        data_available<-intersect(names_low,names_high)

        template[data_available]<-tbl_low[data_available]

        patients_with_low_emt<-template


	}

	}

	if((length_high==length_low) & (length_high<3 & length_low<3)){

	template_high<-c(0,0,0)
	names(template_high)<-c('depletion','gain','LOH_Neutral')

	data_available_high<-intersect(names(template_high),names_high)
        template_high[data_available_high]<-tbl_high[data_available_high]

        patients_with_high_emt<-template_high

	template_low<-c(0,0,0)
	names(template_low)<-c('depletion','gain','LOH_Neutral')

	data_available_low<-intersect(names(template_low),names_low)
        template_low[data_available_low]<-tbl_low[data_available_low]

        patients_with_low_emt<-template_low


	}

        mtxstats<-matrix(c(patients_with_high_emt,patients_with_low_emt),ncol=2)
        colnames(mtxstats)<-c('high','low')
        rownames(mtxstats)<-c('depletion','amplification','neutral')

	# p-value enrichment amplification vs the rest
	mtxstats_amp_vs_all<-rbind(mtxstats['amplification',],colSums(mtxstats[c('depletion','neutral'),]))
	rownames(mtxstats_amp_vs_all)<-c('amplification','not_amp')

	fish_amp<-fisher.test(mtxstats_amp_vs_all)
        p_amp_rest<-fish_amp$p.value
	or_amp_rest<-fish_amp$estimate

        #look only for amplification and not amplification  
        mtxstats2<-mtxstats_amp_vs_all[c('amplification','not_amp'),]
        ratios<-mtxstats2[1,]/mtxstats2[2,]

        ratios_high<-ratios[1]
        ratios_low<-ratios[2]

	# p-value enrichment depletion vs the rest
	mtxstats_dep_vs_all<-rbind(mtxstats['depletion',],colSums(mtxstats[c('amplification','neutral'),]))
        rownames(mtxstats_dep_vs_all)<-c('depletion','not_dep')

	fish_del<-fisher.test(mtxstats_dep_vs_all)
        p_dep_rest<-fish_del$p.value
	or_dep_rest<-fish_del$estimate

	mtxstats2<-mtxstats_dep_vs_all[c('depletion','not_dep'),]
        ratios_dep<-mtxstats2[1,]/mtxstats2[2,]

        ratios_high_dep<-ratios_dep[1]
        ratios_low_dep<-ratios_dep[2]

	RESPVAL_AMP<-c(RESPVAL_AMP,p_amp_rest)
	RESPVAL_DEL<-c(RESPVAL_DEL,p_dep_rest)
	
	RESOR_AMP<-c(RESOR_AMP,or_amp_rest)
	RESOR_DEL<-c(RESOR_DEL,or_dep_rest)
	}

        }
	
	
	stats<-data.frame(genes=GENES_ANALYZED,
			  p_value_amp=RESPVAL_AMP,
			  q_value_amp=p.adjust(RESPVAL_AMP,'BH'),
			  p_value_del=RESPVAL_DEL,
			  q_value_del=p.adjust(RESPVAL_DEL,'BH'),
			  OR_amp_high=RESOR_AMP,
			  OR_del=RESOR_DEL)

	
        return(stats)

}

heatmapResCNV<-function(dictionary_copy_number, annotation_file,current_genes,sig_cnv){

	dictionary_copy_number2<-dictionary_copy_number[,c(1,2,3,27,24,5,28)]

        annotation_file<-annotation_file[annotation_file[,5]%in%current_genes,]
        current_genes<-intersect(current_genes,annotation_file[,5])

        input_annotation<-annotation_file
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)

        colnames(dictionary_copy_number2)[c(1:3)]<-c('chrom','start','end')

        RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number2)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))

        qh<-res_overlaps$queryHits
        sh<-res_overlaps$subjectHits

        dictionary_copy_number2<-as.data.table(dictionary_copy_number2)
        input_annotation<-as.data.table(input_annotation)

        res_query<-dictionary_copy_number2[qh,]
        res_subject<-input_annotation[sh,]

        RES_ANNOTATION_CURRENT_CANCER<-data.frame(res_query,res_subject)
	RES_ANNOTATION_CURRENT_CANCER<-RES_ANNOTATION_CURRENT_CANCER[,which(colnames(RES_ANNOTATION_CURRENT_CANCER)%in%c('samples','GENE_NAME','Segment_Mean','cancer'))]	

	sig_cnv2<-sig_cnv[,c(1,ncol(sig_cnv))]
	colnames(sig_cnv2)[1]<-'GENE_NAME'

	rasig<-merge(RES_ANNOTATION_CURRENT_CANCER,sig_cnv2)




	rasig_amp<-as.data.frame(setDT(rasig[which(rasig$status%in%'amplification'),c(1:3)])[, lapply(.SD, mean), by = .(samples, GENE_NAME)])[,c(1:3)]
	colnames(rasig_amp)[3]<-'copy_number'

	rasig_depletion<-as.data.frame(setDT(rasig[which(rasig$status%in%'deletion'),c(1:3)])[, lapply(.SD, mean), by = .(samples, GENE_NAME)])[,c(1:3)]
	colnames(rasig_depletion)[3]<-'copy_number'
	

	amp_matrix<-as.data.frame.matrix(xtabs(copy_number~GENE_NAME+samples, data=rasig_amp))
	dep_matrix<-as.data.frame.matrix(xtabs(copy_number~GENE_NAME+samples, data=rasig_depletion))

	colnames(amp_matrix)<-gsub(make.names(sapply(strsplit(colnames(amp_matrix),split='-'),'[[',3),unique=T),pattern='\\.',replacement='_')
	colnames(dep_matrix)<-gsub(make.names(sapply(strsplit(colnames(dep_matrix),split='-'),'[[',3),unique=T),pattern='\\.',replacement='_')

	pseudospace_patients_amp<-intersect(gsub(make.names(ps2[order(ps2$pseudospace),2],unique=T),pattern='\\.',replacement='_'),colnames(amp_matrix))
	pseudospace_patients_dep<-intersect(gsub(make.names(ps2[order(ps2$pseudospace),2],unique=T),pattern='\\.',replacement='_'),colnames(dep_matrix))


	amp_matrix2<-amp_matrix[,pseudospace_patients_amp]
	dep_matrix2<-dep_matrix[,pseudospace_patients_dep]

	amp_matrix2<-amp_matrix2[,match(pseudospace_patients_amp,colnames(amp_matrix2))]

	pdf('test.pdf')
	   quantile_breaks <- function(xs, n = nbreaks-1) {
	   breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
   	   breaks[!duplicated(breaks)]
  	   }

        mat_breaks <- quantile_breaks(as.matrix(amp_matrix), n = 100) 
        hmcols <- inferno(length(mat_breaks) - 1)

	ph<-pheatmap(amp_matrix2, useRaster = T, cluster_cols = T,
        cluster_rows = T, show_rownames = F, show_colnames = F,color=hmcols,breaks=mat_breaks)
	ph
	dev.off()

}


heatmapPseudospaceSelectedGenes<-function(mtx_exp,ps,gs){
	       

		ps<-ps[ps[,2]%in%colnames(mtx_exp),]

		ps<-ps[order(ps$pseudospace),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps[,2]%in%colnames(mtx_exp))]
                mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,2],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'

		mtx_exp3<-mtx_exp3[mtx_exp3[,1]%in%gs,]
		rownames(mtx_exp3)<-mtx_exp3[,1]

		spline<-function(x,ps){
		res<-smooth.spline(x,ps[,3])$y
		return(res)
		}

		

		m=mtx_exp3[,-1]		
		m2<-apply(m,1,scale)

	
		bks <- seq(-3, 3, by = 0.2)
 	        hmcols <- inferno(length(bks) - 1)

		 ph<-pheatmap(t(m2), useRaster = T, cluster_cols = F,cluster_rows = T, show_rownames = T, show_colnames = F,color=hmcols,breaks=bks) 


        }

compareProfilesGenes<-function(mtx_exp1,ps1,mtx_exp2,ps2,gs,label1='inner',label2='outer',trim=TRUE,trim_value=0.5){


		mtx_exp1_genes<-data.frame(melt(mtx_exp1[mtx_exp1[,1]%in%gs,]),group=label1)
	       	colnames(mtx_exp1_genes)[2]<-'patients'
		mtx_exp1_genes_ps<-merge(mtx_exp1_genes,ps1,by='patients')
 
		mtx_exp2_genes<-data.frame(melt(mtx_exp2[mtx_exp2[,1]%in%gs,]),group=label2)
                colnames(mtx_exp2_genes)[2]<-'patients'
                mtx_exp2_genes_ps<-merge(mtx_exp2_genes,ps2,by='patients')
		
		input_plot_tot<-rbind(mtx_exp1_genes_ps,mtx_exp2_genes_ps)
	
		input_plot_tot<-input_plot_tot[input_plot_tot[,2]%in%'VIM',]

		input_plot_tot$new_class<-ifelse(input_plot_tot$value>=10,'greater','lower')

		if(trim==TRUE){

		input_plot_tot<-input_plot_tot[which(input_plot_tot$value>=trim_value),]
		
		bp<-ggplot(input_plot_tot,aes(x=pseudospace,y=value,color=group,shape=group))+geom_point(aes(alpha=0.9))+scale_color_manual(values=c('#E69F00','#56B4E9'))+facet_wrap(~genes)+geom_smooth(method='loess',fullrange=T,aes(linetype=group),color='black')+geom_smooth(method='lm',formula=y~1)

		}


		ggboxplot(input_plot_tot,'new_class','value',color='group')
}

plotProfileSingleGenes<-function(mtx_exp,ps,gs,output){

	        ps<-ps[ps[,2]%in%colnames(mtx_exp),]
                ps<-ps[order(ps$pseudospace),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps[,2]%in%colnames(mtx_exp))]
                mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,2],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'

                mtx_exp3<-mtx_exp3[mtx_exp3[,1]%in%gs,]
                rownames(mtx_exp3)<-mtx_exp3[,1]
		
		melt_exp<-melt(mtx_exp3)
		colnames(melt_exp)[2]<-'patients'

		melt_exp2<-merge(melt_exp,ps,by='patients')	
	
		plot_output <- paste(output,"_exp_pseudotime",".png",sep='')

		if(length(unique(melt_exp2$ratio))==2){
		bp<-ggplot(melt_exp2,aes(x=pseudospace,y=value,color=ratio,shape=ratio))+geom_point(aes(alpha=0.9))+scale_color_manual(values=c('#E69F00','#56B4E9'))+facet_wrap(~genesID)+geom_smooth(method='loess',fullrange=T,aes(linetype=ratio),color='black')}else{
		 bp<-ggplot(melt_exp2,aes(x=pseudospace,y=value,color=ratio,shape=ratio))+geom_point(aes(alpha=0.9))+facet_wrap(~genesID)+geom_smooth(method='loess',fullrange=T,aes(linetype=ratio),color='black')
		}
		
		png(plot_output, height=6,width=10,units='in',res=300)
		print(bp)
		dev.off()
		
		if(length(unique(melt_exp2$ratio))==2){
		bp2<-ggboxplot(melt_exp2,x='ratio',y='value',col='ratio',palette=c('#E69F00','#56B4E9'))+stat_compare_means(comparisons=list(c('outer','inner')),method='wilcox.test')+facet_wrap(~genesID)
		} else {
		bp2<-ggboxplot(melt_exp2,x='ratio',y='value',col='ratio')+stat_compare_means(comparisons=list(c('outer','inner')),method='wilcox.test')+facet_wrap(~genesID)
		}
		
		plot_output2 <- paste(output,"_boxplot_exp_pseudotime",".png",sep='')
		
		png(plot_output2,width=8,height=10,units='in',res=300)
		print(bp2)
		dev.off()
}

plotProfileSingleGenesHighEMT<-function(mtx_exp,ps,gs,output){

                ps<-ps[ps[,2]%in%colnames(mtx_exp),]
                ps<-ps[order(ps$pseudospace),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps[,2]%in%colnames(mtx_exp))]
                mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps[,2],colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'

                mtx_exp3<-mtx_exp3[mtx_exp3[,1]%in%gs,]
                rownames(mtx_exp3)<-mtx_exp3[,1]

                melt_exp<-melt(mtx_exp3)
                colnames(melt_exp)[2]<-'patients'

                melt_exp2<-merge(melt_exp,ps,by='patients')
		print(head(melt_exp2))

                plot_output <- paste(output,"_exp_pseudotime",".png",sep='')

                if(length(unique(melt_exp2$ratio))==2){
                bp<-ggplot(melt_exp2,aes(x=pseudospace,y=value,color=ratio,shape=ratio))+geom_point(aes(alpha=0.9))+scale_color_manual(values=c('#E69F00','#56B4E9'))+facet_wrap(~genesID)+geom_smooth(method='loess',fullrange=T,aes(linetype=ratio),color='black')}else{
                 bp<-ggplot(melt_exp2,aes(x=pseudospace,y=value,color=ratio,shape=ratio))+geom_point(aes(alpha=0.9))+facet_wrap(~genesID)+geom_smooth(method='loess',fullrange=T,aes(linetype=ratio),color='black')
                }
                
                png(plot_output, height=6,width=10,units='in',res=300)
                print(bp)
                dev.off()
                
                if(length(unique(melt_exp2$ratio))==2){
                bp2<-ggboxplot(melt_exp2,x='ratio',y='value',col='ratio',palette=c('#E69F00','#56B4E9'))+stat_compare_means(comparisons=list(c('High-EMT','Low-EMT')),method='wilcox.test')+facet_wrap(~genesID)
                } else {
                bp2<-ggboxplot(melt_exp2,x='ratio',y='value',col='ratio')+stat_compare_means(comparisons=list(c('High-EMT','Low-EMT')),method='wilcox.test')+facet_wrap(~genesID)
                }
                
                plot_output2 <- paste(output,"_boxplot_exp_pseudotime",".png",sep='')
                
                png(plot_output2,width=8,height=10,units='in',res=300)
                print(bp2)
                dev.off()
}


plotProfileSingleGenesStage<-function(mtx_exp,ps1,gs,output){

		colnames(mtx_exp)<-gsub(make.names(sapply(strsplit(colnames(mtx_exp),split='\\.'),'[',4),unique=T),pattern='\\.',replacement='_')

                ps<-ps[ps$patients%in%colnames(mtx_exp),]
                ps<-ps[order(ps$pseudospace),]

                # mtx_exp2<-mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])]
                mtx_exp2<-mtx_exp[,c(ps$patients%in%colnames(mtx_exp))]
                mtx_exp3<-data.frame(mtx_exp[,1],mtx_exp2[,match(ps$patients,colnames(mtx_exp2))])
                colnames(mtx_exp3)[1]<-'genesID'

                mtx_exp3<-mtx_exp3[mtx_exp3[,1]%in%gs,]
                rownames(mtx_exp3)<-mtx_exp3[,1]

                melt_exp<-melt(mtx_exp3)
                colnames(melt_exp)[2]<-'patients'

                melt_exp2<-merge(melt_exp,ps,by='patients')

                plot_output <- paste(output,"_exp_stage_pseudotime",".png",sep='')

                bp<-ggplot(melt_exp2,aes(x=pseudospace,y=value,color=ratio_mock,shape=ratio_mock))+geom_point(aes(alpha=0.9))+scale_color_manual(values=c('#E69F00','#56B4E9','#999999','#0DB275'))+facet_wrap(~genesID)+geom_smooth(method='loess',fullrange=T,aes(linetype=ratio),color='black')


                png(plot_output, height=6,width=10,units='in',res=300)
                print(bp)
                dev.off()

                bp2<-ggboxplot(melt_exp2,x='ratio',y='value',col='ratio',palette=c('#E69F00','#56B4E9','#999999'))+stat_compare_means(comparisons=list(c('outer','inner')),method='wilcox.test')+facet_wrap(~genesID)


                plot_output2 <- paste(output,"_boxplot_exp_stage_pseudotime",".png",sep='')

                png(plot_output2,width=8,height=10,units='in',res=300)
                print(bp2)
                dev.off()
}




saveOutputDegs<-function(resDegsAllcancer,string_output,thr){

frequencies_of_degs_along_cancer<-data.frame(rowSums(table(resDegsAllcancer$tumor,resDegsAllcancer$genes))/ncol(table(resDegsAllcancer$tumor,resDegsAllcancer$genes)))

write.table(frequencies_of_degs_along_cancer,file=paste(paste(string_output,'most_frequent_cancers_genes',qvalue,sep='_'),'.txt',sep=''),sep='\t',row.names=T,quote=F)

most_frequent_degs_across_cancers<-data.frame(colSums(table(resDegsAllcancer$tumor,resDegsAllcancer$genes))/nrow(table(resDegsAllcancer$tumor,resDegsAllcancer$genes)))
most_frequent_degs_across_cancers<-data.frame(rownames(most_frequent_degs_across_cancers),most_frequent_degs_across_cancers)
most_frequent_degs_across_cancers<-most_frequent_degs_across_cancers[as.numeric(most_frequent_degs_across_cancers[,2])>=thr,]

write.table(most_frequent_degs_across_cancers,file=paste(paste(string_output,'most_frequent_genes',qvalue,sep='_'),'.txt',sep=''),sep='\t',row.names=F,quote=F)

}


saveOutputDegsCNV<-function(resDegsAllcancer,string_output,thr=0.20){

frequencies_of_degs_along_cancer<-data.frame(rowSums(table(resDegsAllcancer$tumors,resDegsAllcancer$genes))/ncol(table(resDegsAllcancer$tumors,resDegsAllcancer$genes)))

write.table(frequencies_of_degs_along_cancer,file=paste(paste(string_output,'most_frequent_cancers_genes',qvalue,sep='_'),'.cnv.txt',sep=''),sep='\t',row.names=T,quote=F)

most_frequent_degs_across_cancers<-data.frame(colSums(table(resDegsAllcancer$tumors,resDegsAllcancer$genes))/nrow(table(resDegsAllcancer$tumors,resDegsAllcancer$genes)))
most_frequent_degs_across_cancers<-data.frame(rownames(most_frequent_degs_across_cancers),most_frequent_degs_across_cancers)
most_frequent_degs_across_cancers<-most_frequent_degs_across_cancers[as.numeric(most_frequent_degs_across_cancers[,2])>=thr,]

write.table(most_frequent_degs_across_cancers,file=paste(paste(string_output,'most_frequent_genes',qvalue,sep='_'),'.cnv.txt',sep=''),sep='\t',row.names=F,quote=F)

genes_for_amp_depletion<-table(resDegsAllcancer$genes,resDegsAllcancer$status)

genes_for_amp_depletion2<-genes_for_amp_depletion[which(rownames(genes_for_amp_depletion)%in%most_frequent_degs_across_cancers[,1]),]

write.table(most_frequent_degs_across_cancers,file=paste(paste(string_output,'most_frequent_genes',qvalue,sep='_'),'.cnv.txt',sep=''),sep='\t',row.names=T,quote=F)


}



# 










































createMatrixFromSegFilePurity<-function(dictionary_copy_number,annotation_file,current_genes,ps){

        input_annotation<-annotation_file[annotation_file[,5]%in%current_genes,]
        colnames(input_annotation)[c(1:3)]<-c('chrom','start','end')

        annotation_GR<-makeGRangesFromDataFrame(input_annotation)

        colnames(dictionary_copy_number)[c(1:3)]<-c('chrom','start','end')
        RES_ALL_CANCERS_GR<-makeGRangesFromDataFrame(dictionary_copy_number)

        res_overlaps<-as.data.frame(findOverlaps(query=RES_ALL_CANCERS_GR,subject=annotation_GR))

        qh<-res_overlaps$queryHits
        sh<-res_overlaps$subjectHits

        res_query<-dictionary_copy_number[qh,]
        res_subject<-input_annotation[sh,]

        RES_ANNOTATION_CURRENT_CANCER<-data.frame(res_query,res_subject)

        dfcan_signal<-RES_ANNOTATION_CURRENT_CANCER[,c('samples','GENE_NAME','R_abs')]
        dfcan_signal[,1]<-sapply(strsplit(dfcan_signal[,1],split='\\-'),'[',3)

        dfcan2<-as.data.frame(setDT(dfcan_signal)[, lapply(.SD, mean), by = .(samples, GENE_NAME)])

        dfcan3<-dcast(dfcan2, GENE_NAME ~ samples)

        return(dfcan3)
}

plotGenesWithPseudospaceCNV<-function(mtx_exp,mtx_cnv,gs,ps,tab_clinical,clinical_feature,clusters,trim=0.1,output){

                mtx_exp2<-data.frame(mtx_exp[,1],mtx_exp[,c(colnames(mtx_exp)%in%ps[,2])])
                colnames(mtx_exp2)[1]<-'genes'

                mtx_exp2_select<-melt(mtx_exp2[as.character(mtx_exp2[,1])%in%gs,])

                colnames(mtx_exp2_select)[2]<-'patients'

                mtx_exp2_select3<-merge(mtx_exp2_select,ps,by='patients')
                colnames(mtx_exp2_select3)[3]<-'expression'

                mtx_exp2_select4<-merge(mtx_exp2_select3,tab_clinical[,c('samples',clinical_feature)],by.x='patients',by.y='samples',all.x=T)

                resKmean<-data.frame()

                for(i in unique(mtx_exp2_select4$tumors)){
                
		input_for_kmeans<-mtx_exp2_select4[mtx_exp2_select4$tumors==i,c('patients','genes','expression')]

                input_for_kmeans_wide<-reshape(input_for_kmeans,idvar='genes',timevar='patients',direction='wide')
                rownames(input_for_kmeans_wide)<-input_for_kmeans_wide[,1]

		input_for_kmeans_wide[is.na(input_for_kmeans_wide)]<-0

                kmdf<-kmeans(input_for_kmeans_wide[,-1],clusters)

                kmdf<-data.frame(tumors=i,clusters=kmdf$cluster,genes=input_for_kmeans_wide[,1])

                resKmean<-rbind(resKmean,kmdf)

                }

                mtx_exp2_select5<-merge(mtx_exp2_select4,resKmean,by.x=c('genes','tumors'),by.y=c('genes','tumors'),all.x=T)

                pdf(paste(output,'.pdf',sep=''))

                for(cn in unique(mtx_exp2_select4$tumors)) {

                print(cn)

                input_for_plot<-mtx_exp2_select5[mtx_exp2_select5$tumors==cn,]

                input_for_plot<-input_for_plot[which(input_for_plot$expression>trim),]

                input_for_plot[,clinical_feature]<-as.factor(input_for_plot[,clinical_feature])

                q <- ggplot(aes(pseudospace, expression), data = input_for_plot)


                q2 <- q + geom_point(aes(pseudospace,expression,color=tumor_stage),data=input_for_plot)+geom_smooth(aes(color = NULL),method = "loess")+ggtitle(cn)

                print(q2)

                }

                dev.off()

                write.table(resKmean,file=paste(output,'.kmeans.txt',sep=''),sep='\t',row.names=F,quote=F)

                pdf(paste(output,'.kmeans.pdf',sep=''))

                for(cn in unique(mtx_exp2_select5$tumors)) {

                print(cn)

                input_for_plot<-mtx_exp2_select5[mtx_exp2_select5$tumors==cn,]

                input_for_plot<-input_for_plot[which(input_for_plot$expression>trim),]

                input_for_plot[,clinical_feature]<-as.factor(input_for_plot[,clinical_feature])

                q <- ggplot(aes(pseudospace, expression), data = input_for_plot)

                q2 <- q + geom_point(aes(pseudospace,expression,color=tumor_stage),data=input_for_plot)+geom_smooth(aes(color = NULL),method = "loess")+ggtitle(cn) + facet_wrap(~clusters)
                print(q2)

                }

                dev.off()

}


plot_pseudotime_heatmap_binSize <- function (mtx_exp,ps,gs, bin_size = 5,num_clusters=6,nbreaks=20,main)
{
    
    psx<-ps[order(ps$pseudospace),]

    mtx_exp2<-data.frame(mtx_exp[,1],mtx_exp[,c(colnames(mtx_exp)%in%psx[,2])])
    colnames(mtx_exp2)[1]<-'genes'

    mtx_exp2_select<-melt(mtx_exp2[as.character(mtx_exp2[,1])%in%gs,])

    colnames(mtx_exp2_select)[2]<-'patients'
   
    mtx_exp2_select<-data.table(mtx_exp2_select,key='patients')
    ps<-data.table(ps,key='patients')

    mtx_exp2_select3<-merge(mtx_exp2_select,ps,by='patients')
    colnames(mtx_exp2_select3)[3]<-'expression'

    m=mtx_exp2_select3[order(mtx_exp2_select3$pseudospace),c('genes','patients','expression')]
    m2=reshape2::dcast(m, genes~patients, value.var = "expression")[,-1]
 
    m3<-m2[,match(psx[,2],colnames(m2))]
	
    smoothing_function<-function(X,bin){sapply(split(X, (seq_along(X) - 1) %/% bin), mean)}

    test<-apply(m3,1,smoothing_function,bin_size)
    test<-apply(test,1,scale)

   #https://slowkow.com/notes/pheatmap-tutorial/

   pdf('test.pdf')
   quantile_breaks <- function(xs, n = nbreaks-1) {
   breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
   breaks[!duplicated(breaks)]
   }

   test[test>2]<-2
   mat_breaks <- quantile_breaks(test, n = nbreaks)
  

   hmcols <- inferno(length(mat_breaks) - 1)

    pheatmap(test,cluster_rows=F,cluster_cols=F, useRaster = T, show_rownames = F,
        show_colnames = F, cutree_rows = num_clusters, treeheight_row = 10,
        breaks = mat_breaks, fontsize = 6, color = hmcols,
        main = main, annotation_colors = F)


    dev.off() 
}

