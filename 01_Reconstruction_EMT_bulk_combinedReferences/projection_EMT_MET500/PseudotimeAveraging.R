#######################################
#MET500 EMT analysis
#######################################

library(data.table)
library(monocle)
library(DropletUtils)
library(DESeq2)
library(irlba)
library(FNN)

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/MCF10/output_dir/")
input_cp_final <- read.table("KNN_projection_TCGA_MET500_to_MCF10A_treated_cells.txt", header = TRUE,sep = "\t")
input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace + input_cp_final$tgfb_pseudospace

conditions <- c("A549_EGF","A549_TGFB1","A549_TNF","DU145_EGF","DU145_TGFB1","DU145_TNF",
                "MCF7_EGF","MCF7_TGFB1","MCF7_TNF","OVCA420_EGF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/CookReference/output_dir/")
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_MET500_to_",i,"_treated_cells.txt",sep = ""), header = TRUE, sep = "\t")
  print(all(input_cp_final$patients == pseudotime_add$patients))
  input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace + pseudotime_add$mock_pseudospace
  
  
}
summary(input_cp_final$mock_pseudospace)
input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace / 14
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/")
write.table(input_cp_final, "AveragedPseudotime_lowconf.txt",row.names=F,quote=F,sep="\t")









setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/MCF10/output_dir/")
input_cp_final <- read.table("KNN_projection_TCGA_MET500_to_MCF10A_treated_cells.txt", header = TRUE,sep = "\t")
input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace + input_cp_final$tgfb_pseudospace

conditions <- c("A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF",
                "MCF7_TGFB1","MCF7_TNF","OVCA420_TGFB1","OVCA420_TNF")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/CookReference/output_dir/")
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_MET500_to_",i,"_treated_cells.txt",sep = ""), header = TRUE, sep = "\t")
  print(all(input_cp_final$patients == pseudotime_add$patients))
  input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace + pseudotime_add$mock_pseudospace
  
  
}
summary(input_cp_final$mock_pseudospace)
input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace / 10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/")
write.table(input_cp_final, "AveragedPseudotime_midconf.txt",row.names=F,quote=F,sep="\t")








setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/MCF10/output_dir/")
input_cp_final <- read.table("KNN_projection_TCGA_MET500_to_MCF10A_treated_cells.txt", header = TRUE,sep = "\t")
input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace + input_cp_final$tgfb_pseudospace

conditions <- c("A549_TGFB1","DU145_TGFB1",
                "MCF7_TGFB1","OVCA420_TGFB1")

for (i in conditions) {
  
  print(i)
  
  setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/CookReference/output_dir/")
  
  pseudotime_add <- read.table(paste("KNN_projection_TCGA_MET500_to_",i,"_treated_cells.txt",sep = ""), header = TRUE, sep = "\t")
  print(all(input_cp_final$patients == pseudotime_add$patients))
  input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace + pseudotime_add$mock_pseudospace
  
  
}
summary(input_cp_final$mock_pseudospace)
input_cp_final$mock_pseudospace <- input_cp_final$mock_pseudospace / 6
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_MET500/")
write.table(input_cp_final, "AveragedPseudotime_highconf.txt",row.names=F,quote=F,sep="\t")








