# Ridgeline plot pseudotime scores

setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/MCF10/01_Reconstruction_EMT_bulk/output_dir/")
pseudospace_input <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
pseudospace_input <- pseudospace_input[,colnames(pseudospace_input) %in% c("patients","mock","tgfb")]
library(reshape2)
pseudospace_summary <- melt(pseudospace_input, id.vars = c("patients"))
colnames(pseudospace_summary) <- c("Patients","Reference","Pseudotime")
pseudospace_summary <- pseudospace_summary[,order(colnames(pseudospace_summary))]
pseudospace_summary$Reference <- sapply(pseudospace_summary$Reference, function(x)
  ifelse(x %in% "mock","MCF10_mock",
         ifelse(x %in% "tgfb","MCF10_tgfb","other")))

References <- c("A549_TNF","A549_TGFB1","DU145_TNF","DU145_TGFB1",
                "MCF7_TNF","MCF7_TGFB1","OVCA420_TNF","OVCA420_TGFB1")

for (i in References) {

  print(i)  
  
  setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/",i,"/01_Reconstruction_EMT_bulk/output_dir/",sep = ""))
  pseudospace_input <- read.table(paste("KNN_projection_TCGA_to_",i,"_treated_cells_no_correction_primarytumor_withzeros.txt",sep = ''), header = TRUE,sep = "\t")
  pseudospace_input$Reference <- i
  colnames(pseudospace_input) <- c("Patients","Pseudotime","Reference")
  pseudospace_input <- pseudospace_input[,order(colnames(pseudospace_input))]
  pseudospace_summary <- rbind(pseudospace_summary, pseudospace_input)
  
}
#Add combined reference:
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/")
full.ref <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE,sep = "\t")
full.ref <- full.ref[,colnames(full.ref) %in% c("patients","mock")]
full.ref$Reference <- "Combined"
colnames(full.ref) <- c("Patients","Pseudotime","Reference")
pseudospace_summary <- rbind(pseudospace_summary, full.ref)

table(pseudospace_summary$Reference)
pseudospace_summary$Reference <- factor(pseudospace_summary$Reference, levels = c("Combined","A549_TGFB1","A549_TNF","DU145_TGFB1","DU145_TNF","MCF7_TGFB1","MCF7_TNF","MCF10_mock","MCF10_tgfb","OVCA420_TGFB1","OVCA420_TNF"))
library(ggridges)
library(ggplot2)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/VariationBetweenEstimates")
pdf("PseudotimeReferneceRidgelinePlot.pdf", height = 5,width = 10)
ggplot(pseudospace_summary, aes(x = Pseudotime, y = Reference, fill = Reference)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()

