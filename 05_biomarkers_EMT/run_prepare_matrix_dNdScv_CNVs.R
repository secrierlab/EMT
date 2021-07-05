library(data.table)
library(plyr)

setwd("/data/pseudospace/dndscv/hmm_3_states")

list_files<-c("TCGA_projection_mock_HMM_3_withpemt.dndscv_siggenes.0.1.epi.txt",
	      "TCGA_projection_mock_HMM_3_withpemt.dndscv_siggenes.0.1.mes.txt",
	      "TCGA_projection_mock_HMM_3_withpemt.dndscv_siggenes.0.1.mix.txt")

names(list_files)<-c("epithelial","Mesenchymal","mixed")

list_anno<-c("TCGA_projection_mock_HMM_3_withpemt.dndscv_anno.0.1.epi.txt",
              "TCGA_projection_mock_HMM_3_withpemt.dndscv_anno.0.1.mes.txt",
              "TCGA_projection_mock_HMM_3_withpemt.dndscv_anno.0.1.mix.txt")
names(list_anno)<-c("epi","mes","mix")

setwd('/data/pseudospace/HMM')
input_file<-c("HMM_results_nstates_3.txt")

knn_df_tcga<-read.delim(file=input_file)
knn_df_tcga$tumors<-sapply(strsplit(as.character(knn_df_tcga[,1]),split="\\."),"[[",1)
colnames(knn_df_tcga)[1]<-"patients2"

classes_patients<-split(knn_df_tcga[,c("patients2","biological_states","tumors")],f=knn_df_tcga$biological_states)


#
# Prepare the variants data
#

setwd("/data/pseudospace/dndscv/hmm_3_states")

list_results_var<-vector(mode="list",length=length(list_files))

for(var_file in 1:length(list_files)){

	#
	# Load the significant variants
	#
	
	#group_current_var_file<-names(list_files)[var_file]
	#sig_var_input<-fread(list_files[var_file],data.table=F)
	#sig_genes<-unique(sig_var_input[,2])

	#
	# Load the annotation related with the significant variant
	#
	
	group_current_anno_file<-names(list_anno)[var_file]

        anno_var_input<-fread(list_anno[var_file],data.table=F)

	matrix_genes_pts<-table(anno_var_input$sampleID,anno_var_input$gene)
	attributes(matrix_genes_pts)$class <- "matrix" 

	list_results_var[[var_file]]<-data.frame(matrix_genes_pts)
}

res_tot_var<-rbind.fill(list_results_var)
rownames(res_tot_var)<-gsub(make.names(unlist(lapply(list_results_var,rownames)),unique=T),pattern="\\.",replacement="-")
colnames(res_tot_var)<-paste(colnames(res_tot_var),"_dndscv",sep="")

# change NA values with zero
res_tot_var[is.na(res_tot_var)]<-0
res_tot_var[res_tot_var>1]<-1

#
# Prepare the copy-number data
#
setwd("/data/pseudospace/gistic/hmm_3_states")

folders_cnv<-c("epi","mes","mix")

ALL_CNA<-vector(mode="list",length=length(folders_cnv))

for(fc in 1:length(folders_cnv)){
	
	fold_cnv<-folders_cnv[fc]

	print(fold_cnv)

	setwd(paste("/data/pseudospace/gistic/hmm_3_states",fold_cnv,sep="/"))

	# for each group there are multiple directories for cancer type
	sub_folders<-dir()

        # for each group there are multiple directories for cancer type
        res_by_tum<-vector(mode="list",length(sub_folders))

	for(sb in 1:length(sub_folders)){
	       
	sb2<-sub_folders[sb]

	setwd(paste(paste("/data/pseudospace/gistic/hmm_3_states",fold_cnv,sep="/"),sb2,sep="/"))
	
	# get all the lesions
	
	#check the existence of all_lesions.conf_99, if not skip the analysis
	if(length( grep(dir(),pattern="all_lesions.conf_99.txt"))!=0){

	all_lesions_cnv<-fread("all_lesions.conf_99.txt",data.table=F)
	all_lesions_cnv<-all_lesions_cnv[all_lesions_cnv[,9]=="Actual Copy Change Given",]
	regions_to_use<-all_lesions_cnv$Descriptor

	all_lesions_cnv_peak<-all_lesions_cnv[grep(all_lesions_cnv[,1],pattern="Peak"),c(1,2,10:ncol(all_lesions_cnv))]
	all_lesions_cnv_peak[,1]<-sapply(strsplit(all_lesions_cnv_peak[,1],split=" "),"[[",1)
	
	# take the significant samples
	# can exists multiple arms regions with the same bands so i use make.names =T 
	rownames(all_lesions_cnv_peak)<-gsub(make.names(gsub(paste(all_lesions_cnv_peak[,1],all_lesions_cnv_peak[,2],sep="_"),pattern="\\.",replacement="_"),unique=T),pattern="\\.",replacement="")

	all_lesions_cnv_peak<-all_lesions_cnv_peak[,grep(colnames(all_lesions_cnv_peak),pattern="TCGA")]
        
	# first output gistic to use
	colnames(all_lesions_cnv_peak)<-unlist( lapply(strsplit(colnames(all_lesions_cnv_peak),split="-"),FUN=function(x){paste(x[1:4],collapse="-")}))

	# get focal events genes-levels
        focal_cnv<-fread("focal_data_by_genes.txt",data.table=F)
	rownames(focal_cnv)<-paste(focal_cnv[,1],"focal",sep="_")

	# take only the regions in the lesions files because are those one significant
	focal_cnv2<-focal_cnv[which(focal_cnv$Cytoband %in% regions_to_use),-c(1:3)]
	colnames(focal_cnv2)<-unlist( lapply(strsplit(colnames(focal_cnv2),split="-"),FUN=function(x){paste(x[1:4],collapse="-")}))
	
	broad_sig<-fread("broad_significance_results.txt",data.table=F)
	broad_values_by_arm<-fread("broad_values_by_arm.txt",data.table=F)
	colnames(broad_values_by_arm)<-unlist(lapply(strsplit(colnames(broad_values_by_arm),split="-"),FUN=function(x){paste(x[1:4],collapse="-")}))

	arm_amp_significant_regions<-broad_sig[which(broad_sig[,"Amp q-value"]<=0.10),"Arm"]
        arm_del_significant_regions<-broad_sig[which(broad_sig[,"Del q-value"]<=0.10),"Arm"]

	# output arms_events
	arm_amp_final<-broad_values_by_arm[which(broad_values_by_arm[,1] %in% arm_amp_significant_regions),]
	if(nrow(arm_amp_final)!=0){
	rownames(arm_amp_final)<-paste(arm_amp_final[,1],"arm_amplification",sep="_")
	}
	
	arm_del_final<-broad_values_by_arm[which(broad_values_by_arm[,1] %in% arm_del_significant_regions),]
	if(nrow(arm_del_final)!=0){
	rownames(arm_del_final)<-paste(arm_del_final[,1],"arm_deletion",sep="_")
	}

	if(nrow(arm_amp_final)!=0 & nrow(arm_del_final)!=0){
	arms_global2<-rbind(arm_del_final[,-1],arm_amp_final[,-1])
	}
	        
	if(nrow(arm_amp_final)!=0 & nrow(arm_del_final)==0){
	arms_global2<-rbind(arm_amp_final[,-1])
	}

	if(nrow(arm_del_final)!=0 & nrow(arm_amp_final)==0){
        arms_global2<-rbind(arm_del_final[,-1])
        }


	#uniform the columns names for safe
	colnames(focal_cnv2)<-gsub(make.names(colnames(focal_cnv2),unique=T,allow = TRUE),pattern="\\.",replacement="_")
	colnames(all_lesions_cnv_peak)<-gsub(make.names(colnames(all_lesions_cnv_peak),unique=T),pattern="\\.",replacement="_")
	colnames(arms_global2)<-gsub(make.names(colnames(arms_global2),unique=T),pattern="\\.",replacement="_")

	cna_events<-rbind(focal_cnv2,all_lesions_cnv_peak,arms_global2)

	print(dim(cna_events))

	}# close if existence file all_config_99.conf

	# change the columns names cna_events: recreate it

	colnames(cna_events)<- lapply(strsplit(colnames(cna_events),split="_"),FUN=function(x){paste(paste(x[1:4],collapse="-"),sep="_")})

	# create a data.frame in which the first column contains the genomic features
	df_cna<-data.frame(genomic_features= rownames(cna_events),
			   cna_events)
	df_cna[,1]<-as.character(df_cna[,1])
	
	res_by_tum[[sb]]<-df_cna

	}

	#merge the results from all tumors for the current state
	all_genomic_features<-data.frame(unique(unlist(lapply(res_by_tum,FUN=function(x){x[,1]}))))
	colnames(all_genomic_features)[1]<-"genomic_features"
	all_genomic_features[,1]<-as.character(all_genomic_features[,1])

	for(rtm in 1:length(res_by_tum)){
	
	print(rtm)	
	print(dim(res_by_tum[[rtm]]))
	all_genomic_features<-merge(all_genomic_features,res_by_tum[[rtm]],all.x=T,by="genomic_features")
	
	}
	rownames(all_genomic_features)<-all_genomic_features[,1]

	ALL_CNA[[fc]]<-all_genomic_features

	setwd(paste("/data/pseudospace/gistic/hmm_3_states"))

}


# different groups of samples can have the same genomic alterations

all_genomic_features<-data.frame(unique(unlist(lapply(res_by_tum,FUN=function(x){x[,1]}))))      
colnames(all_genomic_features)[1]<-"genomic_features"

merge1<-merge(all_genomic_features,ALL_CNA[[1]],all.x=T,by="genomic_features")
merge2<-merge(all_genomic_features,ALL_CNA[[2]],all.x=T,by="genomic_features")
merge3<-merge(all_genomic_features,ALL_CNA[[3]],all.x=T,by="genomic_features")

#There are patients that could be classified in one group or another if they are replicates
#Therefore the data can be repeated multiple times

#MIX: grep(a,pattern="TCGA.BK.A26L.01A",value=T)
#[1] "UCEC.TCGA.BK.A26L.01A.11R.A277.07"
# grep(b,pattern="TCGA.BK.A26L.01A",value=T)
#character(0)
#MES: grep(c,pattern="TCGA.BK.A26L.01A",value=T)
#[1] "UCEC.TCGA.BK.A26L.01A.11R.A16F.07"
common1<-intersect(colnames(merge1),colnames(merge2))[-1]
common2<-intersect(colnames(merge1),colnames(merge3))[-1]
ucommon<-c(common1,common2)

ALL_CNA2<-cbind(merge1[,-which(colnames(merge1)%in%ucommon)],
		merge2[,-which(colnames(merge1)%in%c("genomic_features",ucommon))],	       
		merge3[,-which(colnames(merge1)%in%c("genomic_features",ucommon))])

rownames(ALL_CNA2)<-ALL_CNA2[,1]
colnames(ALL_CNA2)<-gsub(colnames(ALL_CNA2),pattern="\\.",replacement="-")

#
# Now merge the results of the variants with the results of the copy number by samples
#

input_cna<-data.frame(ID=as.character(colnames(ALL_CNA2[,-1])),t(ALL_CNA2[,-1]))
input_mut<-data.frame(ID=as.character(rownames(res_tot_var)),res_tot_var)

#input_cna[,1]<-sapply(strsplit(as.character(input_cna[,1]),split="-"),"[[",3)
#input_mut[,1]<-sapply(strsplit(as.character(input_mut[,1]),split="-"),"[[",3)

input_tot<-merge(input_cna,input_mut,by="ID",all.x=T,all.y=T)
input_tot[,1]<-unlist(lapply(strsplit(as.character(input_tot[,1]),split="\\-"),FUN=function(x){paste(x[1:4],collapse="-")}))

sub_pseudotime<-knn_df_tcga
sub_pseudotime[,1]<-as.character(sub_pseudotime[,1])
sub_pseudotime[,1]<-unlist(lapply(strsplit(sub_pseudotime[,1],split="\\."),FUN=function(x){paste(x[2:5],collapse="-")}))

colnames(sub_pseudotime)[1]<-"patients"

input_tot_pseudotime<-merge(x=sub_pseudotime,y=input_tot,by.y="ID",by.x="patients")

setwd("/data/pseudospace/ml_for_ppt")

write.table(input_tot_pseudotime,file="input_for_ml_hmm_states_3_mock.txt",sep="\t",row.names=F,quote=F)

input_tot_pseudotime[is.na(input_tot_pseudotime)]<-0

write.table(input_tot_pseudotime,file="input_for_ml_hmm_states_3_mock_v2.txt",sep="\t",row.names=F,quote=F)

setwd("/data/datasets/TCGA/TCGA_supp/aneuploidy_score")
aneuploidy_score<-read.delim(file="aneuploidy_score_taylor_et_al2018.legacy.txt",stringsAsFactors=F)
aneuploidy_score<-aneuploidy_score[,c(1,3,6)]
colnames(aneuploidy_score)[c(2:3)]<-c("aneuploidy","genome_doubling")

input_tot_pseudotime[,1]<-substring(input_tot_pseudotime[,1],1,15)

input_tot_pseudotime_aneuploidy<-merge(x = input_tot_pseudotime,
				       y = aneuploidy_score,
				       by.x = "patients",
				       by.y = "Sample",
				       all.x=T)

input_tot_pseudotime_aneuploidy[is.na(input_tot_pseudotime_aneuploidy)]<-0

setwd("/data/pseudospace/ml_for_ppt")

write.table(input_tot_pseudotime_aneuploidy,file="input_for_ml_hmm_states_3_mock_as_gd.txt",sep="\t",row.names=F,quote=F)
