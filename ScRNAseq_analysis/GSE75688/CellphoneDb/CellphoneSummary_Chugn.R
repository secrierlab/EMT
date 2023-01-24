####################################################
#Number of interactions with differnet TME compnents
####################################################


#Load required packages:
library(reshape)
library(ggplot2)

#Load metadata to see how many cells are in the differnet categories:
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/CellphoneDb/Input/")
meta.data <- read.table("test_meta.txt", header = TRUE, sep = "\t")
table(meta.data$cell_type)


#Load interaction data
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/CellphoneDb/Input/out")
df.net <- read.table("significant_means.txt", header = TRUE, sep = "\t")
df.net <- melt(df.net, id.vars = c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"))


#Check number of interactions
cell_types <- unique(as.character(meta.data$cell_type))
cell_types <- cell_types[!(cell_types %in% c("LowPseudotime_tumour","MidPseudotime_tumour","HighPseudotime_tumour"))]

Interaction_number <- NULL

for (i in cell_types) {
  
  print(i)
  df.net.low <- df.net[df.net$variable %in% c(paste(i,".LowPseudotime_tumour",sep = ""),paste("LowPseudotime_tumour.",i,sep = "")),]
  df.net.low <- na.omit(df.net.low)
  df.net.low <- dim(df.net.low)[1]
  Interaction_number <- c(Interaction_number, df.net.low)
  
  
  df.net.mid <- df.net[df.net$variable %in% c(paste(i,".MidPseudotime_tumour",sep = ""),paste("MidPseudotime_tumour.",i,sep = "")),]
  df.net.mid <- na.omit(df.net.mid)
  df.net.mid <- dim(df.net.mid)[1]
  Interaction_number <- c(Interaction_number, df.net.mid)
  
  
  df.net.high <- df.net[df.net$variable %in% c(paste(i,".HighPseudotime_tumour",sep = ""),paste("HighPseudotime_tumour.",i,sep = "")),]
  df.net.high <- na.omit(df.net.high)
  df.net.high <- dim(df.net.high)[1]
  Interaction_number <- c(Interaction_number, df.net.high)
  
  
  
}

CellType <- NULL

for (i in cell_types) {
  
  print(i)
  
  CellType <- c(CellType, rep(i,3))
  
}

l <- length(cell_types)
PseudotimeStatus <- rep(c("Low","Mid","High"),l)
summary <- data.frame(CellType, PseudotimeStatus, Interaction_number)
summary$PseudotimeStatus <- factor(summary$PseudotimeStatus, levels = c("Low","Mid","High"))


#Boxplot summary
pdf("Chung2017_InteractionNumberSummary.pdf",height = 5,width = 10)
ggplot(summary, aes(x=PseudotimeStatus, y=Interaction_number)) + 
  geom_bar(stat = "identity") +  facet_wrap(.~CellType, nrow = 1) + theme_classic() + ylab("Total Number of Interactions with Cancer Cells") + xlab("Pseudotime Status")
dev.off()
