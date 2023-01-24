########################################
#CellphoneDB input - Breast canccer
########################################


#Load data:
setwd("~/Documents/Dormancy_TME_project_data/scRNAseq/Qian2020/Breast/")
load("expr_data.RData")
#Load annotation:
load("QS_all_cell.RData")
all(colnames(expr.data) == merged.data$Cell)
#Load EMT annotation
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/QianBreast/01_Reconstruction_EMT_bulk/")
pseudotime <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE, sep = "\t")

#Add quiescence score annotation:
table(merged.data$CellType)
B_cell <- merged.data[merged.data$CellType %in% "B_cell",]
B_cell <- as.character(B_cell$Cell)
DC <- merged.data[merged.data$CellType %in% "DC",]
DC <- as.character(DC$Cell)
EC <- merged.data[merged.data$CellType %in% "EC",]
EC <- as.character(EC$Cell)
Fibroblast <- merged.data[merged.data$CellType %in% "Fibroblast",]
Fibroblast <- as.character(Fibroblast$Cell)
Mast <- merged.data[merged.data$CellType %in% "Mast",]
Mast <- as.character(Mast$Cell)
Myeloid <- merged.data[merged.data$CellType %in% "Myeloid",]
Myeloid <- as.character(Myeloid$Cell)
T_cell <- merged.data[merged.data$CellType %in% "T_cell",]
T_cell <- as.character(T_cell$Cell)
#Split cancer cells according to pseudotime categories
pseudotime <- pseudotime[order(pseudotime$mock),]
split_pseudotime <- split(pseudotime, rep(1:3, length.out = nrow(pseudotime), each = ceiling(nrow(pseudotime)/3)))
pseudotime_low <- split_pseudotime[[1]]
pseudotime_mid <- split_pseudotime[[2]]
pseudotime_high <- split_pseudotime[[3]]
pseudotime_low <- as.character(pseudotime_low$patients)
pseudotime_mid <- as.character(pseudotime_mid$patients)
pseudotime_high <- as.character(pseudotime_high$patients)
merged.data$cell_type <- sapply(merged.data$Cell, function(x)
  ifelse(x %in% pseudotime_low,"LowPseudotime_tumour",
         ifelse(x %in% B_cell, "B_cell",
                ifelse(x %in% DC,"DC",
                       ifelse(x %in% EC,"EC",
                              ifelse(x %in% Fibroblast, "Fibroblast",
                                     ifelse(x %in% Mast, "Mast",
                                            ifelse(x %in% Myeloid, "Myeloid",
                                                   ifelse(x %in% T_cell, "T_cell",
                                                          ifelse(x %in% pseudotime_mid, "MidPseudotime_tumour","HighPseudotime_tumour"))))))))))
table(merged.data$cell_type)
merged.data <- merged.data[,colnames(merged.data) %in% c("Cell","cell_type")]
all(merged.data$Cell == colnames(expr.data))



########################################################
##Analysis on a patient-by-patient basis
########################################################

Gene <- rownames(expr.data)
Gene <- data.frame(Gene)
expr.data <- cbind(Gene, expr.data)
rownames(expr.data) <- NULL


##################

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/QianBreast/CellphoneDb/Input/")

#Save:
write.table(merged.data, file = "test_meta.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(expr.data, file = "test_counts.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


