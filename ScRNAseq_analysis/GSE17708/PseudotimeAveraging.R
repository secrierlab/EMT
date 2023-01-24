#######################################
#Combining pseudotime estimates
#######################################



####Pseudotime mid confidence
#Load MCF10 pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/MCF10/output_dir/")
pseudotime <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
pseudotime$mock <- pseudotime$mock + pseudotime$tgfb

conditions <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF",
                "MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/Cook/output_dir/")
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ''), header = TRUE,sep = "\t")
  print(all(pseudotime$patients == pseudotime_add$patients))
  pseudotime$mock <- pseudotime$mock + pseudotime_add$pseudospace
  
  
}

summary(pseudotime$mock)
pseudotime$mock <- pseudotime$mock / 10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/")
write.table(pseudotime, "AveragedPseudotime_midconf_withzeros.txt",row.names=F,quote=F,sep="\t")





##############################
#Pseudotime across timepoints
library(ggplot2)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/ScRNAseqAnalysis/GSE17708/")
pseudotime <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE, sep = "\t")
pseudotime$patients
pseudotime$Timepoint <- c(rep(0,3),rep(0.5,3), rep(1,3),rep(2,2),rep(4,3),rep(8,3),rep(16,3),rep(24,3),rep(72,3))
colnames(pseudotime)
pseudotime$Timepoint <- factor(pseudotime$Timepoint, levels = c(0,0.5,1,2,4,8,16,24,72))
pseudotime$Rep <- c(1,2,3,1,2,3,1,2,3,1,2,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
pseudotime$Rep <- factor(pseudotime$Rep, levels = c(1,2,3))

pdf("GSE17708_timepoint_vs_pseudotime.pdf",height = 4,width = 8)
ggplot(pseudotime, aes(x=Timepoint, y=mock, shape = Rep, color = Rep)) + 
  geom_point() + xlab("Timepoint (hours)") + ylab("Pseudotime") + theme_classic()
dev.off()

