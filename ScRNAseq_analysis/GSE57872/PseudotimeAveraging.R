#######################################
#Combining pseudotime estimates
#######################################



####Pseudotime mid confidence
#Load MCF10 pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_POG/MCF10/output_dir/")
pseudotime <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
pseudotime$mock <- pseudotime$mock + pseudotime$tgfb

conditions <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF",
                "MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_POG/Cook/output_dir/")
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ''), header = TRUE,sep = "\t")
  print(all(pseudotime$patients == pseudotime_add$patients))
  pseudotime$mock <- pseudotime$mock + pseudotime_add$pseudospace
  
  
}

summary(pseudotime$mock)
pseudotime$mock <- pseudotime$mock / 10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_POG/")
write.table(pseudotime, "AveragedPseudotime_midconf_withzeros.txt",row.names=F,quote=F,sep="\t")


