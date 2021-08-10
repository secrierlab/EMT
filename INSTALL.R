install.packages(c("mclust","clustvarsel","phenopath","vcd","rootSolve","cluster",
                   "factoextra","mltest","irlba","FNN","igraph","depmixS4","pheatmap",
                   "diagram","tidyverse","FactoMineR","dotwhisker","interactions",
                   "circlize","parallel","FSelector","nnet","gridExtra","reshape2",
                   "jtools","ggstance","broom.mixed","MASS","ggfortify","lme4","MuMIn",
                   "biglasso","doParallel","glmnetUtils","ROCR","MLeval","caret","ranger",
                   "fastshap","plotrix","data.table","VennDiagram","forcats","ggforce","gtools",
                   "ggpmisc","dplyr","grid","finalfit","knitr","kableExtra","survival","survminer","readxl",
                   "tidyr","glmnet","reshape","viridis","ggridges","cowplot","plyr","corrplot","data.table",
                   "rstatix","ggpubr","vcd","RColorBrewer","openxlsx","patchwork","ggplot2"),dependencies=T)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install(c("DropletUtils","DESeq2",
                       "monocle","TSCAN","slingshot",
                       "ComplexHeatmap","GSVA","VennDetail",
                       "biomaRt","depmap","sva","TCGAbiolinks"))


