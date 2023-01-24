#######################################
#CCLE pseudotime averaging
#######################################

library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/MCF10/output_dir/")
load("MCF10_mock_pseudotime.RData")
final_pseudotime <- pseudotimeCCLEderivedMCF
load("MCF10_tgfb_pseudotime.RData")
final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime + pseudotimeCCLEderivedMCF$all_pseudotime

conditions <- c("A549_EGF","A549_TGFB1","A549_TNF","DU145_EGF","DU145_TGFB1","DU145_TNF",
                "MCF7_EGF","MCF7_TGFB1","MCF7_TNF","OVCA420_EGF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/CookReference/output_dir/")
  
  load(paste("Pseudotime_",i,".RData",sep = ""))
  final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime + pseudotimeCCLEderivedMCF$all_pseudotime
  
  
}
summary(final_pseudotime$all_pseudotime)
final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime / 14
pseudotimeCCLEderivedMCF <- final_pseudotime
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/")
save(pseudotimeCCLEderivedMCF, file = "pseudotime_lowconf.RData")








setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/MCF10/output_dir/")
load("MCF10_mock_pseudotime.RData")
final_pseudotime <- pseudotimeCCLEderivedMCF
load("MCF10_tgfb_pseudotime.RData")
final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime + pseudotimeCCLEderivedMCF$all_pseudotime

conditions <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF",
                "MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/CookReference/output_dir/")
  
  load(paste("Pseudotime_",i,".RData",sep = ""))
  final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime + pseudotimeCCLEderivedMCF$all_pseudotime
  
  
}
summary(final_pseudotime$all_pseudotime)
final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime / 10
pseudotimeCCLEderivedMCF <- final_pseudotime
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/")
save(pseudotimeCCLEderivedMCF, file = "pseudotime_midconf.RData")







setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/MCF10/output_dir/")
load("MCF10_mock_pseudotime.RData")
final_pseudotime <- pseudotimeCCLEderivedMCF
load("MCF10_tgfb_pseudotime.RData")
final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime + pseudotimeCCLEderivedMCF$all_pseudotime

conditions <- c("A549_TGFB1","DU145_TGFB1",
                "MCF7_TGFB1","OVCA420_TGFB1")


for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/CookReference/output_dir/")
  
  load(paste("Pseudotime_",i,".RData",sep = ""))
  final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime + pseudotimeCCLEderivedMCF$all_pseudotime
  
  
}
summary(final_pseudotime$all_pseudotime)
final_pseudotime$all_pseudotime <- final_pseudotime$all_pseudotime / 6
pseudotimeCCLEderivedMCF <- final_pseudotime
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_tarjectory_CCLE/")
save(pseudotimeCCLEderivedMCF, file = "pseudotime_highconf.RData")


