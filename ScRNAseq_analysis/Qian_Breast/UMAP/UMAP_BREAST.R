#########
###UMAP:
#########

#Identification of highly variable features:
bc.data <- FindVariableFeatures(bc.data, selection.method = "vst", nfeatures = 2000)

#Scaling the data:
all.genes <- rownames(bc.data)
bc.data <- ScaleData(bc.data, features = all.genes)

#Perform linear dimensionality reduction
bc.data <- RunPCA(bc.data, features = VariableFeatures(object = bc.data))

#UMAP:
bc.data <- RunUMAP(bc.data, dims = 1:50,seed.use = 123)
UMAP_coordinates <- bc.data@reductions$umap
UMAP_coordinates <- UMAP_coordinates@cell.embeddings

#Save UMAP coordinates:
setwd("~/Documents/Dormancy_TME_project/scRNAseq_datasets/UMAP/Qian_2021_variable_genes")
#setwd("/home/annawiecek/Qian_UMAP")

save(UMAP_coordinates, file = "UMAP_breast_50.RData")





######################################
#Do this on a patient by patient basis 
Patients <- c(41,42,43,44,45,46,47,48,49,50,51,52,53,54)


for (i in Patients) {
  
  print(i)
  
  bc.data.patient <- subset(x = bc.data.save, subset = Patient == i)
  
  
  #Identification of highly variable features:
  bc.data.patient <- FindVariableFeatures(bc.data.patient, selection.method = "vst", nfeatures = 2000)
  
  #Scaling the data:
  all.genes <- rownames(bc.data.patient)
  bc.data.patient <- ScaleData(bc.data.patient, features = all.genes)
  
  #Perform linear dimensionality reduction
  bc.data.patient <- RunPCA(bc.data.patient, features = VariableFeatures(object = bc.data.patient))
  
  
  #UMAP:
  bc.data.patient <- RunUMAP(bc.data.patient, dims = 1:50,seed.use = 123)
  UMAP_coordinates <- bc.data.patient@reductions$umap
  UMAP_coordinates <- UMAP_coordinates@cell.embeddings
  
  #Save UMAP coordinates:
  save(UMAP_coordinates, file = paste("UMAP_breast_patient_",i,"_50.RData",sep = ""))
  
  
  
}





##########################################
#Plots:

#Main:
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/QianBreast/UMAP/")
load('UMAP_breast_50.RData')
#Load annotation
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/QianBreast/CellphoneDb/Input/")
UMAP_coordinates <- data.frame(UMAP_coordinates)
UMAP_coordinates$Cell <- rownames(UMAP_coordinates)
anno <- read.table("test_meta.txt", header = TRUE, sep = "")
merge.data <- merge(UMAP_coordinates, anno,
                    by.x = "Cell", by.y = "Cell")
setwd("~/Documents/Dormancy_TME_project_data/scRNAseq/Qian2020/Breast/")
anno <- read.csv("BC_metadata.csv", header = TRUE)
merge.data <- merge(merge.data, anno,
                    by.x = "Cell", by.y = "Cell")
merge.data$PatientNumber <- factor(merge.data$PatientNumber)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/QianBreast/UMAP/")
library(ggplot2)
library(RColorBrewer)
pdf("BreastCancer_PatientNumber.pdf",height = 10,width = 10)
ggplot(merge.data, aes(x=UMAP_1, y=UMAP_2, colour = PatientNumber)) +
  geom_point(size = 1.5, alpha = 0.6) + theme_classic() 
dev.off()
table(merge.data$cell_type)

cols <- c("B_cell" = "#A6CEE3", "DC" = "#1F78B4", "EC" = "#B2DF8A", "Fibroblast" = "#33A02C","HighPseudotime_tumour" = "#E31A1C", 
          "LowPseudotime_tumour" = "#FDBF6F", "Mast" = "#FB9A99","MidPseudotime_tumour" = "#FF7F00", "Myeloid" = "#CAB2D6", "T_cell" = "#6A3D9A")
pdf("BreastCancer_CellType.pdf",height = 10,width = 10)
p <- ggplot(merge.data, aes(x=UMAP_1, y=UMAP_2, colour = cell_type)) +
  geom_point(size = 1.5, alpha = 0.6) + theme_classic() 
p + scale_color_manual(values = cols)
dev.off()

#Patient-by-patient basis


for (i in Patients) {
  
  print(i)
  
  setwd("~/Documents/Dormancy_TME_project/scRNAseq_datasets/UMAP/Qian_2021_variable_genes")
  
  load(paste("UMAP_breast_patient_",i,"_50.RData",sep = ""))
  UMAP_coordinates <- data.frame(UMAP_coordinates)
  UMAP_coordinates$Cell <- rownames(UMAP_coordinates)
  setwd("~/Documents/EMT_paper_revision/Qian_Breast/CellphoneDb/")
  anno <- read.table("test_meta.txt", header = TRUE, sep = "")
  merge.data <- merge(UMAP_coordinates, anno,
                      by.x = "Cell", by.y = "Cell")
  setwd("~/Documents/Dormancy_TME_project_data/scRNAseq/Qian_2020/Breast/")
  anno <- read.csv("BC_metadata.csv", header = TRUE)
  merge.data <- merge(merge.data, anno,
                      by.x = "Cell", by.y = "Cell")
  setwd("~/Documents/EMT_paper_revision/Qian_Breast/UMAP/Figures/")
  
  pdf(paste("BreastCancer_CellType_patient_",i,".pdf",sep = ""),height = 10,width = 10)
  p <- ggplot(merge.data, aes(x=UMAP_1, y=UMAP_2, colour = cell_type)) +
    geom_point(size = 1.5, alpha = 0.6) + theme_classic() 
  print(p + scale_color_manual(values = cols))
  dev.off()
  
  
}