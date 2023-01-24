#####################
#CellChat - GSE75688
####################


#Load required libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)


############################################################################
#PART I: Data input and processing and initialization of the CellChat object
############################################################################

#Load data:
setwd("~/Documents/Dormancy_TME_project_data/scRNAseq/GSE75688/Expression/")
annotation <- read.table("GSE75688_final_sample_information.txt", header = TRUE,sep = "\t")
TPM.data <- read.table("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",header = TRUE, sep = "\t")
##Reorganise the datafame:
TPM.data$gene_name <- as.character(TPM.data$gene_name)
TPM.data$gene_id <- NULL
TPM.data$gene_type <- NULL
TPM.data <- TPM.data %>% group_by(gene_name) %>% mutate_each(funs(mean)) %>% distinct
TPM.data <- data.frame(TPM.data)
rownames(TPM.data) <- TPM.data$gene_name
TPM.data$gene_name <- NULL
TPM.data <- as.matrix(TPM.data)
TPM.data <- log2(TPM.data + 1) #Data needs to be log-transformed 
TPM.data <- data.frame(TPM.data)



#Select only primary tumour samples and samples which have not undergone chemotherapy:
annotation$sample <- as.character(annotation$sample)
annotation$patient.sample <- sapply(annotation$sample, function(x)
  strsplit(x,"_")[[1]][1])
annotation$type <- as.character(annotation$type)
annotation <- annotation[annotation$type %in% "SC",]
annotation$Label <- annotation$index3


#Add quiescence score annotation:
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/01_Reconstruction_EMT_bulk/")
pseudotime <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE, sep = "\t")
pseudotime <- pseudotime[order(pseudotime$mock),]
split_pseudotime <- split(pseudotime, rep(1:3, length.out = nrow(pseudotime), each = ceiling(nrow(pseudotime)/3)))
pseudotime_low <- split_pseudotime[[1]]
pseudotime_mid <- split_pseudotime[[2]]
pseudotime_high <- split_pseudotime[[3]]
pseudotime_low <- as.character(pseudotime_low$patients)
pseudotime_mid <- as.character(pseudotime_mid$patients)
pseudotime_high <- as.character(pseudotime_high$patients)
annotation_normal <- annotation[annotation$index %in% "nonTumor",]
annotation_tumour <- annotation[annotation$index %in% "Tumor",]
#Identify quiescent samples
annotation_tumour$Label <- sapply(annotation_tumour$sample, function(x)
  ifelse(x %in% pseudotime_low, "LowPseudotime_tumour",
         ifelse(x %in% pseudotime_high, "HighPseudotime_tumour",
                ifelse(x %in% pseudotime_mid, "MidPseudotime_tumour","Other"))))
annotation <- rbind(annotation_normal, annotation_tumour)
table(annotation$Label)
TPM.data <- TPM.data[,colnames(TPM.data) %in% annotation$sample]
rownames(annotation) <- annotation$sample
all(rownames(annotation) == colnames(TPM.data))
TPM.data <- TPM.data[,order(colnames(TPM.data))]
annotation <- annotation[order(rownames(annotation)),]
all(rownames(annotation) == colnames(TPM.data))
TPM.data <- as.matrix(TPM.data) #Ultimate format is a matrix with samples as colnames and genes as rownames


###Format for Cellphonedb:
Gene <- rownames(TPM.data)
Gene <- data.frame(Gene)
TPM.data <- cbind(Gene, TPM.data)
rownames(TPM.data) <- NULL

annotation <- annotation[,colnames(annotation) %in% c("sample","Label")]
colnames(annotation) <- c("Cell","cell_type")
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/CellphoneDb/Input/")
#Save:
write.table(annotation, file = "test_meta.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(TPM.data, file = "test_counts.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
