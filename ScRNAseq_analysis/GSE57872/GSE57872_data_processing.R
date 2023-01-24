######################################
#GSE57872 - select tumour cells
#####################################



##Load the data 
setwd("~/Documents/EMT_paper_revision/GSE57872/Data/")
TPM.data <- read.table("GSE57872_GBM_data_matrix.txt",header = TRUE, sep = "\t")

#Select cells of interes
sample <- colnames(TPM.data)
sample <- sample[-1]
annotation <- data.frame(sample)
annotation$Patient <- sapply(annotation$sample, function(x)
  strsplit(x,"_")[[1]][1])
table(annotation$Patient)
annotation <- annotation[annotation$Patient %in% c("MGH26","MGH28","MGH29","MGH30","MGH31","MGH26Tumor","MGH28Tumor","MGH29Tumor","MGH30Tumor","MGH31Tumor"),]
annotation$sample <- as.character(annotation$sample)
samples <- annotation$sample
TPM.data <- TPM.data[,colnames(TPM.data) %in% c(samples,"X")]

#Save refined expr and anno
setwd("~/Documents/EMT_paper_revision/GSE57872/Data/")
save(TPM.data, file = "TPM_data.RData")
save(annotation, file = "anno.RData")
