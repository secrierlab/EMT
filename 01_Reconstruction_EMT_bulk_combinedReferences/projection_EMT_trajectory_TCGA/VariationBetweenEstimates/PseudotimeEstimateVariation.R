####################################################################
#Correlation between pseudotime estimates from different references
####################################################################

library(ggpubr)
library(corrplot)


#Correlation between MCF10 reference estimate and estimates from all references 
#Load estimates from full reference
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/")
full.ref <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE,sep = "\t")
#Load estimates from MCF10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/MCF10/01_Reconstruction_EMT_bulk/output_dir/")
mcf10.ref <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
all(full.ref$patients == mcf10.ref$patients)
full.ref$MCF10_mock <- mcf10.ref$mock
full.ref$MCF10_tgfb <- mcf10.ref$tgfb
full.ref$MCF10 <- (full.ref$MCF10_mock + full.ref$MCF10_tgfb) /2
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/VariationBetweenEstimates")
pdf("Full_vs_MCF10_reference_pseudotime_correlations.pdf", height = 8,width = 8)
print(ggscatter(full.ref, x = "mock", y="MCF10_mock", cor.coef = TRUE, add = "reg.line") + ylab("MCF10 mock pseudotime reference") + xlab("Full refernece"))
print(ggscatter(full.ref, x = "mock", y="MCF10_tgfb", cor.coef = TRUE, add = "reg.line") + ylab("MCF10 tgfb pseudotime reference") + xlab("Full refernece"))
print(ggscatter(full.ref, x = "mock", y="MCF10", cor.coef = TRUE, add = "reg.line") + ylab("MCF10 pseudotime references (mock+tgfb)") + xlab("Full refernece"))
dev.off()



#########
#Correlation between current pseudotime and all other estimates 
#Create results table
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/")
full.ref <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE,sep = "\t")
Induction_method <- c("TNF","TGFB","mock")
Results_corr <- data.frame(Induction_method)
Results_pval <- data.frame(Induction_method)
#MCF10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/MCF10/01_Reconstruction_EMT_bulk/output_dir/")
mcf10.ref <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
Correlation <- NA
Pvalue <- NA
cor <- cor.test(full.ref$mock, mcf10.ref$tgfb)
Correlation <- c(Correlation, cor$estimate)
Pvalue <- c(Pvalue, cor$p.value)
cor <- cor.test(full.ref$mock, mcf10.ref$mock)
Correlation <- c(Correlation, cor$estimate)
Pvalue <- c(Pvalue, cor$p.value)
Results_corr$MCF10 <- Correlation
Results_pval$MCF10 <- Pvalue
#Other references
references <- c("A549","DU145","MCF7","OVCA420")
conditions <- c("TNF","TGFB1")
for (i in references) {
  
  print(i)
  Correlation <- NULL
  Pvalue <- NULL
  for (a in conditions) {
    setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/",i,"_",a,"/01_Reconstruction_EMT_bulk/output_dir/",sep = ""))
    pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_",a,"_treated_cells_no_correction_primarytumor_withzeros.txt",sep = ''), header = TRUE,sep = "\t")
    cor <- cor.test(full.ref$mock, pseudotime_add$pseudospace)
    Correlation <- c(Correlation, cor$estimate)
    Pvalue <- c(Pvalue, cor$p.value)
    
  }
  
  Correlation <- c(Correlation,NA)
  Pvalue <- c(Pvalue, NA)
  Results_corr[[i]] <- Correlation
  Results_pval[[i]] <- Pvalue

}
#Plot
rownames(Results_corr) <- Results_corr$Induction_method
Results_corr$Induction_method <- NULL
rownames(Results_pval) <- Results_pval$Induction_method
Results_pval$Induction_method <- NULL
Results_corr <- as.matrix(Results_corr)
Results_pval <- as.matrix(Results_pval)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/VariationBetweenEstimates")
pdf("Full_vs_individual_reference_correlations.pdf",height = 5, width = 5)
corrplot(Results_corr, method = "circle", na.label.col = "white", tl.col = "black", tl.srt = 45, p.mat = Results_pval, sig.level = 0.05, insig = "blank",addCoef.col = "black")
dev.off()








###################################################################
#Correlations acorss individual tissues 
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/")
full.ref <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE,sep = "\t")
full.ref$MCF7 <- 0
full.ref$A549 <- 0
full.ref$DU145 <- 0
full.ref$OVCA420 <- 0
references <- c("A549","DU145","MCF7","OVCA420")
conditions <- c("TNF","TGFB1")
for (i in references) {
  
  for (a in conditions) {
    setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/",i,"_",a,"/01_Reconstruction_EMT_bulk/output_dir/",sep = ""))
    pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_",a,"_treated_cells_no_correction_primarytumor_withzeros.txt",sep = ''), header = TRUE,sep = "\t")
    full.ref[[i]] <- full.ref[[i]] + pseudotime_add$pseudospace

  }
  
}
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/MCF10/01_Reconstruction_EMT_bulk/output_dir/")
mcf10.ref <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
full.ref$MCF10 <- mcf10.ref$mock + mcf10.ref$tgfb
full.ref$MCF10_MCF7 <- full.ref$MCF10 + full.ref$MCF7
full.ref$MCF10_MCF7 <- full.ref$MCF10_MCF7 /4
full.ref$A549 <- full.ref$A549 /2
full.ref$DU145 <- full.ref$DU145 /2
full.ref$OVCA420 <- full.ref$OVCA420 /2
#Correaltions:
full.ref$patients <- gsub('\\.', '-', full.ref$patients)
full.ref$CancerType <- sapply(full.ref$patients, function(x)
  strsplit(x,"-")[[1]][1])
full.ref <- full.ref[order(full.ref$CancerType),]
CancerType <- unique(full.ref$CancerType)
CancerType <- CancerType[CancerType %in% c("BRCA","LUAD","PRAD","OV")]
Results_corr <- data.frame(c("Pancancer",CancerType))
Results_pval <- data.frame(c("Pancancer",CancerType))
References <- c("MCF10_MCF7","A549","DU145","OVCA420")
colnames(Results_corr)[1] <- "CancerType" 
colnames(Results_pval)[1] <- "CancerType"
for (a in References) {
  
  print(a)
  Correlation <- NULL
  Pvalue <- NULL
  cor <- cor.test(full.ref$mock, full.ref[[a]])
  Correlation <- c(Correlation, cor$estimate)
  Pvalue <- c(Pvalue, cor$p.value)
  
  for (i in CancerType) {
    
    selected.data <- full.ref[full.ref$CancerType %in% i,]
    cor <- cor.test(selected.data$mock, selected.data[[a]])
    Correlation <- c(Correlation, cor$estimate)
    Pvalue <- c(Pvalue, cor$p.value)
    
  }
  Results_corr[[a]] <- Correlation
  Results_pval[[a]] <- Pvalue
  
}
rownames(Results_corr) <- Results_corr$CancerType
Results_corr$CancerType <- NULL
rownames(Results_pval) <- Results_pval$CancerType
Results_pval$CancerType <- NULL
Results_corr <- as.matrix(Results_corr)
Results_pval <- as.matrix(Results_pval)
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/VariationBetweenEstimates")
pdf("Full_vs_cancerwise_reference_correlations.pdf",height = 5, width = 5)
corrplot(Results_corr, method = "circle", na.label.col = "white", tl.col = "black", tl.srt = 45, p.mat = Results_pval, sig.level = 0.05, insig = "blank",addCoef.col = "black")
dev.off()




###################################################################
#Mean change in pseudotime values:
###################################################################
#Create results table
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/")
full.ref <- read.table("AveragedPseudotime_midconf_withzeros.txt", header = TRUE,sep = "\t")
#MCF10
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/MCF10/01_Reconstruction_EMT_bulk/output_dir/")
mcf10.ref <- read.table("proj_pseudospace_withzeros.txt", header = TRUE,sep = "\t")
full.ref$MCF10_mock <- full.ref$mock - mcf10.ref$mock
full.ref$MCF10_tgfb <- full.ref$mock - mcf10.ref$tgfb
#Other references
references <- c("A549","DU145","MCF7","OVCA420")
conditions <- c("TNF","TGFB1")
for (i in references) {
  
  print(i)
  for (a in conditions) {
    setwd(paste("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/",i,"_",a,"/01_Reconstruction_EMT_bulk/output_dir/",sep = ""))
    pseudotime_add <- read.table(paste("KNN_projection_TCGA_to_",i,"_",a,"_treated_cells_no_correction_primarytumor_withzeros.txt",sep = ''), header = TRUE,sep = "\t")
    full.ref[[paste(i,a,sep = "_")]] <- full.ref$mock - pseudotime_add$pseudospace

  }
  
}
full.ref <- full.ref[,colnames(full.ref) %in% c("patients","MCF10_mock","MCF10_tgfb","A549_TNF","A549_TGFB1","DU145_TNF",
                                                "DU145_TGFB1","MCF7_TNF","MCF7_TGFB1,OVCA420_TNF","OVCA420_TGFB1")]
library(reshape)
full.ref <- melt(full.ref, id.vars = "patients")
colnames(full.ref) <- c("Patients","Reference","PseudotimeDifference")
setwd("~/Documents/EMT_paper_revision/TCGA_PurityScaledAnalysis/MCF10_Cook_averaged/01_Reconstruction_EMT_bulk/projection_EMT_trajectory_TCGA/VariationBetweenEstimates")
pdf("BoxplotPseudotimeDifferences.pdf",height = 5, width = 8)
ggplot(full.ref, aes(x=Reference, y=PseudotimeDifference)) +
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme_classic()
dev.off()

