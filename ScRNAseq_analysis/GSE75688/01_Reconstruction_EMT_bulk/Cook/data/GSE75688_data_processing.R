######################################
#GSE75688 - select tumour cells
#####################################



##Load the data 
setwd("~/Documents/EMT_paper_revision/Data")
annotation <- read.table("GSE75688_final_sample_information.txt", header = TRUE,sep = "\t")
TPM.data <- read.table("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",header = TRUE, sep = "\t")


##Reorganise the datafame:
TPM.data$gene_name <- as.character(TPM.data$gene_name)
TPM.data$gene_id <- NULL
TPM.data$gene_type <- NULL
colnames(TPM.data)[1] <- "X"


#Select only tumour cells
annotation$sample <- as.character(annotation$sample)
annotation$type <- as.character(annotation$type)
annotation <- annotation[annotation$index %in% "Tumor",]
samples <- annotation$sample
TPM.data <- TPM.data[,colnames(TPM.data) %in% c(samples,"X")]

#Save refined expr and anno
setwd("~/Documents/EMT_paper_revision/Data")
save(TPM.data, file = "TPM_data.RData")
save(annotation, file = "anno.RData")
