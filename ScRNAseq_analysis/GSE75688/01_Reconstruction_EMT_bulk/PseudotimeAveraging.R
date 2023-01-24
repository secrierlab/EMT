#######################################
#Combining pseudotime estimates
#######################################



####Pseudotime mid confidence
#Load MCF10 pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/01_Reconstruction_EMT_bulk/MCF10/output_dir/")
pseudotime <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
pseudotime$mock <- pseudotime$mock + pseudotime$tgfb

conditions <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF",
                "MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/01_Reconstruction_EMT_bulk/Cook/output_dir/")
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ''), header = TRUE,sep = "\t")
  print(all(pseudotime$patients == pseudotime_add$patients))
  pseudotime$mock <- pseudotime$mock + pseudotime_add$pseudospace
  
  
}

summary(pseudotime$mock)
pseudotime$mock <- pseudotime$mock / 10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/01_Reconstruction_EMT_bulk/")
write.table(pseudotime, "AveragedPseudotime_midconf_withzeros.txt",row.names=F,quote=F,sep="\t")






########################
#Extra analysis
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE75688/01_Reconstruction_EMT_bulk/")
pseudospace_input<-read.delim(file="AveragedPseudotime_midconf_withzeros.txt")
pseudospace_input$Patient <- sapply(pseudospace_input$patients, function(x)
  strsplit(x,"_")[[1]][1])
pseudospace_input$SampleType <- sapply(pseudospace_input$patients, function(x)
  strsplit(x,"_")[[1]][2])
table(pseudospace_input$Patient)
table(pseudospace_input$SampleType)
pseudospace_input$Stage <- sapply(pseudospace_input$Patient, function(x)
  ifelse(x %in% c("BC03LN","BC07LN"),"Metastases","Primary"))
pseudospace_input$patients2 <- pseudospace_input$patients
pseudospace_input_sort<-pseudospace_input[order(pseudospace_input$mock,decreasing=F),]


#sc vs bulk
table(pseudospace_input_sort$SampleType)
bulk_patients <- pseudospace_input_sort[pseudospace_input_sort$SampleType %in% "Pooled",]
bulk_patients <- unique(as.character(bulk_patients$Patient))
mean_sc <- NULL
bulk <- NULL 
for (i in bulk_patients) {
  
  print(i)
  test <- pseudospace_input_sort[pseudospace_input_sort$Patient %in% i,]
  test1 <- test[!(test$SampleType %in% "Pooled"),]
  test1 <- mean(test1$mock)
  test2 <- test[test$SampleType %in% "Pooled",]
  test2 <- mean(test2$mock)
  mean_sc <- c(mean_sc, test1)
  bulk <- c(bulk, test2)
  
}
summary <- data.frame(mean_sc, bulk)
library(ggpubr)
pdf("GSE75688_bulk_vs_scRNAseq_pseudotime.pdf", width = 5, height = 5)
p <- ggscatter(summary, x = "bulk", y = "mean_sc",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               cor.coef = TRUE, alpha = 0.1, size = 2)
p + labs(x = "Bulk sample peseudotime", y = "mean scRNAseq sample pseudotime") 
dev.off()




