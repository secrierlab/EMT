#######################################
#Combining pseudotime estimates
#######################################


####Pseudotime low confidence
#Load MCF10 pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/MCF10/01_Reconstruction_EMT_bulk/output_dir")
pseudotime <- read.table("proj_pseudospace.txt", header = TRUE,sep = "\t")
pseudotime$mock <- pseudotime$mock + pseudotime$tgfb

conditions <- c("A549_EGF","A549_TGFB1","A549_TNF","DU145_EGF","DU145_TGFB1","DU145_TNF",
                "MCF7_EGF","MCF7_TGFB1","MCF7_TNF","OVCA420_EGF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/",i,"/01_Reconstruction_EMT_bulk/output_dir",sep = ""))
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ''), header = TRUE,sep = "\t")
  print(all(pseudotime$patients == pseudotime_add$patients))
  pseudotime$mock <- pseudotime$mock + pseudotime_add$pseudospace
  
  
}

summary(pseudotime$mock)
pseudotime$mock <- pseudotime$mock / 14
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged")
write.table(pseudotime, "AveragedPseudotime_lowconf.txt",row.names=F,quote=F,sep="\t")







####Pseudotime mid confidence
#Load MCF10 pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/MCF10/01_Reconstruction_EMT_bulk/output_dir")
pseudotime <- read.table("proj_pseudospace.txt", header = TRUE,sep = "\t")
pseudotime$mock <- pseudotime$mock + pseudotime$tgfb

conditions <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF",
                "MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/",i,"/01_Reconstruction_EMT_bulk/output_dir",sep = ""))
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ''), header = TRUE,sep = "\t")
  print(all(pseudotime$patients == pseudotime_add$patients))
  pseudotime$mock <- pseudotime$mock + pseudotime_add$pseudospace
  
  
}

summary(pseudotime$mock)
pseudotime$mock <- pseudotime$mock / 10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged")
write.table(pseudotime, "AveragedPseudotime_midconf.txt",row.names=F,quote=F,sep="\t")




####Pseudotime high confidence
#Load MCF10 pseudotime estimates
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/MCF10/01_Reconstruction_EMT_bulk/output_dir")
pseudotime <- read.table("proj_pseudospace.txt", header = TRUE,sep = "\t")
pseudotime$mock <- pseudotime$mock + pseudotime$tgfb

conditions <- c("A549_TGFB1","DU145_TGFB1",
                "MCF7_TGFB1","OVCA420_TGFB1")

for (i in conditions) {
  
  print(i)
  
  setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/",i,"/01_Reconstruction_EMT_bulk/output_dir",sep = ""))
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor.txt",sep = ''), header = TRUE,sep = "\t")
  print(all(pseudotime$patients == pseudotime_add$patients))
  pseudotime$mock <- pseudotime$mock + pseudotime_add$pseudospace
  
  
}

summary(pseudotime$mock)
pseudotime$mock <- pseudotime$mock / 6
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged")
write.table(pseudotime, "AveragedPseudotime_highconf.txt",row.names=F,quote=F,sep="\t")


