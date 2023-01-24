###############################
#GSE17708 Expression datasetup:
###############################

library(dplyr)

setwd("~/Documents/EMT_paper_revision/GSE17708/")
tcga <- read.csv("GSE17708_Keshamouni_TGFB1_logs.csv", header = TRUE)
tcga[1:10,1:10]
colnames(tcga)[1] <- "gene_symbol"
tcga <- tcga %>% group_by(gene_symbol) %>% mutate_each(funs(mean)) %>% distinct
tcga[1:10,1:10]
tcga <- data.frame(tcga)
save(tcga, file = "GSE17708_expr_data.RData")
