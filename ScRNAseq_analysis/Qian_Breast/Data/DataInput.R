######################################
#Qian Breast cancer - select tumour cells
#####################################



##Load the data 
setwd("~/Documents/Dormancy_TME_project_data/scRNAseq/Qian_2020/Breast/")
load("expr_data.RData")
load("QS_all_cell.RData")
table(merged.data$CellType)
merged.data <- merged.data[merged.data$CellType %in% "Cancer",]
expr.data <- expr.data[,colnames(expr.data) %in% merged.data$Cell]
dim(expr.data)


##Reorganise the datafame:
X <- rownames(expr.data)
expr.data <- data.frame(expr.data)
expr.data <- cbind(X, expr.data)


#Save refined expr and anno
setwd("~/Documents/EMT_paper_revision/Qian_Breast/Data/")
save(expr.data, file = "expr_data.RData")
annotation <- merged.data
save(annotation, file = "anno.RData")
