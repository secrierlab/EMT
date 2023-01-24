#####################
#GSE75688 - consensusTME
####################

library(dplyr)

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
annotation <- annotation[annotation$type %in% "Bulk",]
annotation$Label <- annotation$index3
annotation <- annotation[annotation$Label %in% "Tumor",]
TPM.data <- TPM.data[,colnames(TPM.data) %in% annotation$sample]
TPM.data <- as.matrix(TPM.data)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/consensusTME")
save(TPM.data, file = "TPM_expr_bulk_samples.RData")

#Consensus TME
library(ConsensusTME)
immune_infiltration <- ConsensusTME::consensusTMEAnalysis(TPM.data, cancer = "BRCA", statMethod = "ssgsea")
immune_infiltration <- data.frame(t(immune_infiltration))
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/consensusTME")
save(immune_infiltration, file = "immune_infiltration.RData")
#Add pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/01_Reconstruction_EMT_bulk/")
pseudospace_input<-read.delim(file="AveragedPseudotime_midconf_withzeros.txt")
pseudospace_input <- pseudospace_input[pseudospace_input$patients %in% rownames(immune_infiltration),]
immune_infiltration$Sample <- rownames(immune_infiltration)
pseudospace_input <- pseudospace_input[,colnames(pseudospace_input) %in% c("patients","mock")]
colnames(pseudospace_input) <- c("Sample","Pseudotime")
immune_infiltration <- merge(immune_infiltration, pseudospace_input,
                             by.x = "Sample", by.y = "Sample")
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/consensusTME")
save(immune_infiltration, file = "immune_infiltration.RData")




