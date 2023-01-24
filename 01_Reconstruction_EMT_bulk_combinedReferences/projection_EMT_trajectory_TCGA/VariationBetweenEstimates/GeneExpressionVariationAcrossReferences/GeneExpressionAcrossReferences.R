library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)
library(reshape)
library(ggpubr)

#Load list of EMT marker genes:
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/Data")
emt_markers <- read.table("EMT_and_pEMT_markers.txt", header = TRUE,sep = "\t")
emt_markers <- emt_markers$genes
expr.data <- data.frame(emt_markers)

#Load and combine data:
references <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF","MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")
i <- references[1]
#Load data:
input_string<-paste(i,".rds",sep="")
input <-readRDS(input_string)
no_sti_emt<-rownames(input@meta.data[grep(input@meta.data$Time,pattern="_rm",invert=T),])
matrix_temp<-input@assays$RNA@scale.data
print(dim(matrix_temp))
mock_norm2 <-data.frame(gene_symbol=rownames(input@assays$RNA@scale.data),matrix_temp)
colnames(mock_norm2)[1]<-"genes.gene_short_name"
mock_pData <-data.frame(cell=rownames(input@meta.data),sample=i,Pseudotime=as.numeric(input@meta.data$Pseudotime))
mock_norm2 <- mock_norm2[,colnames(mock_norm2) %in% c("genes.gene_short_name",no_sti_emt)]
mock_pData <- mock_pData[mock_pData$cell %in% no_sti_emt,]
mock_pData$Pseudotime<-(mock_pData$Pseudotime/max(mock_pData$Pseudotime))*100
mock_norm2 <- mock_norm2[mock_norm2$genes.gene_short_name %in% emt_markers,]
expr.data <- merge(expr.data, mock_norm2,
                   by.x = "emt_markers", by.y = "genes.gene_short_name", all = TRUE)
anno <- mock_pData

references <- references[-1]
for (i in references) {
  
  input_string<-paste(i,".rds",sep="")
  input <-readRDS(input_string)
  no_sti_emt<-rownames(input@meta.data[grep(input@meta.data$Time,pattern="_rm",invert=T),])
  matrix_temp<-input@assays$RNA@scale.data
  print(dim(matrix_temp))
  mock_norm2 <-data.frame(gene_symbol=rownames(input@assays$RNA@scale.data),matrix_temp)
  colnames(mock_norm2)[1]<-"genes.gene_short_name"
  mock_pData <-data.frame(cell=rownames(input@meta.data),sample=i,Pseudotime=as.numeric(input@meta.data$Pseudotime))
  mock_norm2 <- mock_norm2[,colnames(mock_norm2) %in% c("genes.gene_short_name",no_sti_emt)]
  mock_pData <- mock_pData[mock_pData$cell %in% no_sti_emt,]
  mock_pData$Pseudotime<-(mock_pData$Pseudotime/max(mock_pData$Pseudotime))*100
  mock_norm2 <- mock_norm2[mock_norm2$genes.gene_short_name %in% emt_markers,]
  expr.data <- merge(expr.data, mock_norm2,
                     by.x = "emt_markers", by.y = "genes.gene_short_name", all = TRUE)
  anno <- rbind(anno,mock_pData)
  
  
}
#Add MCF10 data
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/MCF10/01_Reconstruction_EMT_bulk/data")
load("mock_norm2.RData")
load("mock_pData.RData")
mock_pData$sample <- "MCF10_mock"
anno <- rbind(anno,mock_pData)
mock_norm2$genes.id <- NULL
mock_norm2 <- mock_norm2[mock_norm2$genes.gene_short_name %in% emt_markers,]
expr.data <- merge(expr.data, mock_norm2,
                   by.x = "emt_markers", by.y = "genes.gene_short_name", all = TRUE)
load("tgfb_norm2.RData")
load("tgfb_pData.RData")
tgfb_pData$sample <- "MCF10_TGFB"
anno <- rbind(anno,tgfb_pData)
tgfb_norm2$genes.id <- NULL
tgfb_norm2 <- tgfb_norm2[tgfb_norm2$genes.gene_short_name %in% emt_markers,]
expr.data <- merge(expr.data, tgfb_norm2,
                   by.x = "emt_markers", by.y = "genes.gene_short_name", all = TRUE)
rownames(expr.data) <- expr.data$emt_markers
expr.data$emt_markers <- NULL
expr.data <- data.frame(t(expr.data))
expr.data$Sample <-  gsub('\\.', '-', row.names(expr.data))
anno$cell <- gsub('\\.', '-', anno$cell)
all(anno$cell == expr.data$Sample)
total.data <- merge(anno, expr.data,
                    by.x = "cell", by.y = "Sample")
total.data[is.na(total.data)] <- 0



##################################################################
#Summarise results
Summary <- data.frame(emt_markers)
#Summarise dropout:
Dropout <- NULL
for (i in emt_markers) {
  
  print(i)
  selected.data <- total.data[total.data[[i]] > 0,]
  references <- unique(selected.data$sample)
  dropout <- ((length(references))/10) * 100
  Dropout <- c(Dropout, dropout)
  
}
Summary$PercentageReportedDatasets <- Dropout
Summary <- rbind(Summary, Summary)

######################
#Cook referneces:
Cook_median <- NULL
Cook_sd <- NULL
for (i in emt_markers) {
  
  print(i)
  selected.data <- total.data[total.data[[i]] > 0,]
  references <- unique(selected.data$sample)
  selected.data <- total.data[!(total.data$sample %in% c("MCF10_TGFB","MCF10_mock")),]
  selected.data <- selected.data[selected.data$sample %in% references,]
  median <- median(selected.data[[i]])
  sd <- sd(selected.data[[i]])
  Cook_median <- c(Cook_median, median)
  Cook_sd <-c(Cook_sd,sd)
}




######################
#McFaline referneces:
McFaline_median <- NULL
McFaline_sd <- NULL
for (i in emt_markers) {
  
  print(i)
  selected.data <- total.data[total.data[[i]] > 0,]
  references <- unique(selected.data$sample)
  selected.data <- total.data[total.data$sample %in% c("MCF10_TGFB","MCF10_mock"),]
  selected.data <- selected.data[selected.data$sample %in% references,]
  median <- median(selected.data[[i]])
  sd <- sd(selected.data[[i]])
  McFaline_median <- c(McFaline_median, median)
  McFaline_sd <-c(McFaline_sd,sd)
}
Summary$Dataset <- c(rep("Cook",32),rep("McFalineFigueroa",32))
Summary$MedianExpression <- c(Cook_median, McFaline_median)
Summary$ExpressionSD <- c(Cook_sd, McFaline_sd)
Summary[is.na(Summary)] <- 0

#######################
#PLOT
ggplot(Summary, aes(x=MedianExpression, y=ExpressionSD, size = PercentageReportedDatasets, label = emt_markers)) +
  geom_point(alpha=0.7) + facet_wrap(~Dataset,scales = "free")

library(ggrepel)
dodge <- position_dodge(1)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/VariationBetweenEstimates/GeneExpressionVariationAcrossReferences")
pdf("EMT_GeneExpressionVariation.pdf",height = 10, width = 10)
ggplot(Summary, aes(x = MedianExpression, y = ExpressionSD, label = emt_markers, size = PercentageReportedDatasets)) +
  geom_point(
    position = dodge,
    alpha = 0.5
  ) +
  geom_text_repel(
    position = dodge,
    size = 4,
    alpha = 0.9,
    segment.size = .25,
    segment.alpha = .8,
    force = 1
  ) + theme_classic() + facet_wrap(~Dataset, scales = "free")
dev.off()
