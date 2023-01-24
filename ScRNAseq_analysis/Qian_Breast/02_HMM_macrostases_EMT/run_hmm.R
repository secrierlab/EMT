library(data.table)
library(monocle)
library(irlba)
library(FNN)
library(TSCAN)
library(igraph)
library(slingshot)
library(RColorBrewer)
library(depmixS4)
library(pheatmap)
library(ComplexHeatmap)

setwd("~/Documents/EMT_paper_revision/02_HMM_macrostases_EMT")
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)
load("TCGA_exp_rid.RData") #This is a toy example of TCGA expression data

input_file<-c("proj_pseudospace_mcf_mock_mock_Mclust3_with_pemt.txt")

setwd(input_dir)
pseudospace_input<-read.delim(file="proj_pseudospace_mcf_tgfb_mock_Mclust3_with_pemt.txt")

pseudospace_input_sort<-pseudospace_input[order(pseudospace_input$mock,decreasing=F),]

setwd(input_dir)
markers_genes_read<-read.table(file="EMT_and_pEMT_markers.txt",header=T)
markers_genes<-markers_genes_read[,1]

TCGA_GEXP_ALL_rid<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,]

pseudospace_input_sort<-pseudospace_input_sort[which(pseudospace_input_sort$patients2%in%colnames(TCGA_GEXP_ALL_rid)),]

TCGA_GEXP_ALL_sort<-log(t(TCGA_GEXP_ALL_rid[,match(pseudospace_input_sort$patients2,colnames(TCGA_GEXP_ALL_rid))])+1,2)
colnames(TCGA_GEXP_ALL_sort)<-TCGA_GEXP_ALL[TCGA_GEXP_ALL[,1]%in%markers_genes,1]

#TCGA_GEXP_ALL_sort<-log(TCGA_GEXP_ALL_sort+1,2)
#zscore<-function(x){(x-mean(x))/sd(x)}
#TCGA_GEXP_ALL_sort<-apply(TCGA_GEXP_ALL_sort,1,zscore)

setwd(output_dir)

set.seed(12345)

library(glmnet)    

cvfit <-cv.glmnet(x = TCGA_GEXP_ALL_sort, y = pseudospace_input_sort$mock, alpha = 1)

pdf("mock_glmnet.pdf")
plot(cvfit)
dev.off()

dfcoef<-as.matrix(coef(cvfit,s= "lambda.min"))
genes_to_use<- names(dfcoef[dfcoef[,1]!=0,])[-1]

list_HMM<- lapply(as.list(paste(colnames(TCGA_GEXP_ALL_sort[,genes_to_use]),"~1")),as.formula)

list_families<-vector(mode="list",length(genes_to_use))

for(lf in 1:length(genes_to_use)){

	list_families[[lf]]<-gaussian()
}


for(state in 3){

print(state)

nstates=state

input_HMM<-data.frame(time=pseudospace_input_sort$mock,TCGA_GEXP_ALL_sort[,genes_to_use])

HMM<-depmix(list_HMM,input_HMM[,-1],nstates=nstates,family=list_families)

HMMfit<-fit(HMM, verbose = FALSE) #fit our model to the data set

#use also summary summary(HMMfit, which = "transition") to get the transitions

idx_trans<-which(names(getpars(HMMfit))=="")

transition_matrix<-matrix(getpars(HMMfit)[idx_trans],byrow=T,nrow=nstates)
rownames(transition_matrix)<-paste("from",1:nstates,sep="")
colnames(transition_matrix)<-paste("to",1:nstates,sep="")

gp <- graph.adjacency(transition_matrix, mode = "directed", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)

        
pdf(paste("HMM_graph_transition_mock",nstates,"pdf",sep="."),width=12)
plot(gp,edge.label=round(E(gp)$weight, 3))
dev.off()

library(diagram)

pdf(paste("HMM_graph_transition_mock_v2",nstates,"pdf",sep="."))
plotmat(round(transition_matrix,3),
	lwd = 1, 
	box.lwd = 2, 
        cex.txt = 0.7, 
        box.size=0.1, 
	box.prop=0.5,
	box.type = "circle", 
        arr.length=.1,
        arr.width=.1,
        self.cex = .3,
        self.shifty = -.01,
        self.shiftx = .15,
        main = "",
	relsize=0.9,
	box.col=c('#e69f00','#56b4e9','#ee0c0c'))
dev.off()


write.table(transition_matrix,file=paste("HMM_graph_transition_mock_v3_transition_matrix",nstates,"txt",sep="."),sep="\t",quote=F)

pdf(paste("HMM_graph_transition_mock_v3",nstates,"pdf",sep="."))
plotmat(t(transition_matrix),
        lwd = 1,
        box.lwd = 2,
        cex.txt = 0.7,
        box.size=0.1,
        box.prop=0.5,
        box.type = "circle",
        arr.length=.1,
        arr.width=.1,
        self.cex = .3,
        self.shifty = -.01,
        self.shiftx = .15,
        main = "",
        relsize=0.9,
        box.col=c('#e69f00','#56b4e9','#ee0c0c'))
dev.off()


pdf(paste("HMM_graph_transition_mock_mst",nstates,"pdf",sep="."),width=12)
plot(dp_mst,edge.label=round(E(dp_mst)$weight, 3))
dev.off()



HMMpost<-posterior(HMMfit)

output_HMM<-data.frame(samples=rownames(TCGA_GEXP_ALL_sort),HMM_states=as.character(HMMpost[,1]),pseudospace=pseudospace_input_sort$mock,stage=pseudospace_input_sort$tumor_stage)
output_HMM$HMM_states<-as.character(output_HMM$HMM_states)
output_HMM$biological_states<-output_HMM$HMM_states
output_HMM$biological_states<-gsub(output_HMM$biological_states,pattern="3",replacement="mes")
output_HMM$biological_states<-gsub(output_HMM$biological_states,pattern="2",replacement="epi")
output_HMM$biological_states<-gsub(output_HMM$biological_states,pattern="1",replacement="mix")
write.table(output_HMM,file=paste("HMM_results_nstates_",nstates,".txt",sep=""),row.names=F,quote=F,sep="\t")

TCGA_GEXP_ALL_sort2<-TCGA_GEXP_ALL_sort[,genes_to_use]
zscore<-function(x){(x-mean(x))/sd(x)}
TCGA_GEXP_ALL_sort2<-t(apply(TCGA_GEXP_ALL_sort2,2,zscore))


res_HMM<-data.frame(cbind(t(TCGA_GEXP_ALL_sort2),time=pseudospace_input_sort$mock,state=HMMpost[,1]))
res_HMM$state<-as.factor(res_HMM$state)

library(ggplot2)


resHMM_for_heatmap<-cbind(samples=rownames(TCGA_GEXP_ALL_sort),HMMpost)

#
# Create heatmaps to check the expression of the genes along the states
#

	annotation_colors=list(tumors=c("ACC"="#996b2f","BLCA"="#724fb3","BRCA"="#ba449c",
                                  "CESC"="#adba3d","CHOL"="#8f2e43","COAD"="#191970",
                                  "ESCA"="#5060ac","HNSC"="#65b06f","KICH"="#B266FF",
                                  "KIRC"="#969841","KIRP"="#bf6d29","LIHC"="#d1972c",
                                  "LUAD"="#36dee6","LUSC"="#cc4b3e","OV"="#d97eb8",
                                  "PAAD"="#c43d7e","PRAD"="#d077d8","READ"="#FFB266","SKCM"="#452c7c",
                                  "STAD"="#43c29e","THYM"="#de6991","THCA"="#757ee8",
                                  "UCS"="#88321a","UCEC"="#75b550","UVM"="#8f2931"),
                         tumor_stage=c("low_stage"="green3","high_stage"="#e69f00","metastatic"="#ee0c0c","not_reported"="gainsboro"))


	row_colors=list(type=c("Mesenchymal_marker"="firebrick4","Epithelial_marker"="dodgerblue3","pEMT"="orange2"))
  
	input_cp_scale2<-pseudospace_input[order(pseudospace_input$mock,decreasing=F),]
  
	ann_samples<-input_cp_scale2[input_cp_scale2$patients2%in%colnames(TCGA_GEXP_ALL_sort2),c(3,11)]
	ann_samples[,1]<-as.character(ann_samples[,1])
	ann_samples[,2]<-as.character(ann_samples[,2])
	colnames(ann_samples)[1:2]<-c("tumors","tumor_stage")

	ann_samples$tumor_stage<-gsub(ann_samples$tumor_stage,pattern=" ",replacement="_")

	row_colors=list(type=c("Mesenchymal_marker"="firebrick4","Epithelial_marker"="dodgerblue3","pEMT"="#e69f00"))

	max_col<-ncol(input_cp_scale2)

	ann_genes<-markers_genes_read[markers_genes_read[,1]%in%rownames(TCGA_GEXP_ALL_sort2),]
	colnames(ann_genes)<-c("genes","type")

	order_samples_colors_genes<-match(rownames(TCGA_GEXP_ALL_sort2),ann_genes[,1])

	ann_genes2<-data.frame(ann_genes[order_samples_colors_genes,2])
	colnames(ann_genes2)<-"type"

	row_ann<-rowAnnotation(df=ann_genes2,col=row_colors)

	epi_genes<-intersect(genes_to_use, markers_genes_read[markers_genes_read$status=="Epithelial_marker",1])
	mes_genes<-intersect(genes_to_use, markers_genes_read[markers_genes_read$status=="Mesenchymal_marker",1])

	# Compute an EMT score of the patients
	# VIM + CDH2 + FOXC2 + SNAI1 + SNAI2 + TWIST1 + FN1 + ITGB6 + MMP2 + MMP3 + MMP9 + SOX10 + GCS − CDH1 − DSP − OCLN
	# https://www.nature.com/articles/s41598-018-21061-1

	all_scores_emt<-NULL

	for(iemt in 1:ncol(TCGA_GEXP_ALL_sort2)){

        mean_epi<-mean(TCGA_GEXP_ALL_sort2[epi_genes,iemt])
        mes_epi<-mean(TCGA_GEXP_ALL_sort2[mes_genes,iemt])

        score_emt<-mes_epi-mean_epi

        all_scores_emt<-c(all_scores_emt,score_emt)

	}


	library(circlize)

	col_fun = colorRamp2(c(min(all_scores_emt), 0, max(all_scores_emt)), c("deepskyblue3", "gainsboro", "darkorange2"))

	col_ann<-HeatmapAnnotation(emt_score=all_scores_emt,
				   tumors=ann_samples[,1],
				   stage=ann_samples[,2],
				   col=list(emt_score=col_fun,
           tumors=annotation_colors[[1]],
           stage=annotation_colors[[2]]))

	library(ComplexHeatmap)
	library(forcats)

	pdf(paste("HMM_heatmap_mock",nstates,"pdf",sep="."),width=12)
	p<-Heatmap(TCGA_GEXP_ALL_sort2,
    column_split=as.factor(HMMpost$state),
    show_column_names=F,
		left_annotation = row_ann,
		top_annotation = col_ann,
		row_km=3)
	print(p)
	dev.off()

	scores_EMT<-data.frame(tumors= sapply(strsplit(colnames(TCGA_GEXP_ALL_sort2),split="\\."),"[[",1),
			       samples=colnames(TCGA_GEXP_ALL_sort2),
			       score_emt=all_scores_emt,
			       states=HMMpost$state)

	write.table(scores_EMT,file=paste("HMM_results_nstates_tumors_for_states",nstates,".withEMT.txt",sep=""),row.names=T,quote=F,sep="\t")

	tumors_for_states<-table(scores_EMT$tumors,scores_EMT$states)
	write.table(tumors_for_states,file=paste("HMM_results_nstates_tumors_for_states",nstates,".txt",sep=""),row.names=T,quote=F,sep="\t")

        pdf(paste("HMM_boxplot_EMT",nstates,"pdf",sep="."),width=12)
        p<-ggplot(scores_EMT, aes(x=reorder(tumors,score_emt,FUN=mean), y=score_emt, fill=tumors)) + geom_boxplot()
        print(p+scale_fill_manual(values=annotation_colors[[1]])+geom_hline(yintercept=0)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
        dev.off()

	pdf(paste("HMM_boxplot_EMT_v2",nstates,"pdf",sep="."),width=12)
	scores_EMT$states<-factor(scores_EMT$states,levels=c("2","1","3"))
	p<-ggplot(scores_EMT, aes(x=reorder(tumors,score_emt,FUN=mean), y=score_emt, fill=tumors)) + geom_boxplot()+theme_classic()+facet_wrap(~states)
	print(p+scale_fill_manual(values=annotation_colors[[1]])+geom_hline(yintercept=0)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()

        #
        # Plot scores EMT along pseudotime 
        #
	pseudospace_with_score<-merge(pseudospace_input_sort,scores_EMT,by.x="patients2",by.y="samples")

	pseudospace_with_score$states<-as.factor(pseudospace_with_score$states)
	                
	pdf(paste("HMM_pseudotime_vs_scores_emt_mock",nstates,"pdf",sep="."),width=12,pointsize=12)
	p<-ggplot(pseudospace_with_score, aes(x=mock, y=score_emt,color=states)) +
	geom_point(shape=18,size=2)+scale_x_reverse()+ scale_color_manual(values=c('#56b4e9','#e69f00','#ee0c0c'))+geom_smooth(method = "lm", formula = y ~ poly(x, 10),se=TRUE, linetype="dashed",color="blue")+geom_hline(yintercept=0, linetype="dashed", color = "black")+ theme_classic()
	print(p)
	dev.off()

	#
	# heatmap expression of specifics markers along the pseudotime
	#
	list_mark_emt<-colnames(TCGA_GEXP_ALL_sort)

	all_markers<-data.frame()

	for(mark_emt in 1:length(list_mark_emt)){
	
  exp_current_mark<-TCGA_GEXP_ALL_sort[,mark_emt]

  df_exp_current_mark<-data.frame(samples=names(exp_current_mark),genes=list_mark_emt[mark_emt],exp=exp_current_mark)

  df_exp_current_mark_pseudo<-merge(pseudospace_with_score,df_exp_current_mark,by.x="patients2",by.y="samples")
	
  all_markers<-rbind(all_markers,df_exp_current_mark_pseudo)

	}


	pdf(paste("HMM_pseudotime_vs_scores_emt_mock_MARKERS",nstates,"pdf",sep="."),width=12,height=8,pointsize=12)
	
	df_exp_current_mark_pseudo_selected<-all_markers[all_markers$genes%in%c("CDH1","CRB3","DSP","CDH2","FN1","VIM"),]
	df_exp_current_mark_pseudo_selected$genes<-as.character(df_exp_current_mark_pseudo_selected$genes)
  df_exp_current_mark_pseudo_selected$genes<-factor(df_exp_current_mark_pseudo_selected$genes,levels=c("CDH1","CRB3","DSP","CDH2","FN1","VIM"))


  p<-ggplot(df_exp_current_mark_pseudo_selected, 
  aes(x=mock, y=exp,color=states))+geom_point(shape=17,size=2)+ylim(0,max(TCGA_GEXP_ALL_sort))+facet_wrap(.~genes,nrow=3,ncol=3)+scale_x_reverse()+ scale_color_manual(values=c('#e69f00','#56b4e9','#ee0c0c'))+geom_smooth(aes(group=states,linetype=states),method = "lm", formula = y ~ poly(x,4),se=TRUE,color="black")+ scale_linetype_manual(values=c("twodash","dotted","solid"))+theme_classic()
	print(p)	
	
	dev.off()




}




