source("prepare_data.R")
source("Linear_regression.R")
source("plotting.R")
library(readr)

## Kidney Data ##
kid_CNA <- read.delim("Project_Data/kirc_CNv_core.txt", header = TRUE, sep = "\t", dec = ".")
Significant_kid <- read.csv("Project_Data/Genes.txt", sep="") # Significant Genes 
kid_gene_exp <- read.delim("Project_Data/Cancer.txt")

# prepare Kidney Data
Final_compined_Kidney_df <- new.Prepare_data( kid_gene_exp , kid_CNA , Significant_kid)

# Drop columns with 50% zeros or NAN
Cleaned_Kidney_df <- new.Drop_zeros_NA(Final_compined_Kidney_df,0.50)

# linear regression on each gene without feature selection
model_1 <- new.MultipleLinearRegression(Cleaned_Kidney_df,"NAV2")
model_2 <- new.MultipleLinearRegression(Cleaned_Kidney_df,"AEN")
model_3 <- new.MultipleLinearRegression(Cleaned_Kidney_df,"PFKFB2")
model_4 <- new.MultipleLinearRegression(Cleaned_Kidney_df,"BSPRY")
model_5 <- new.MultipleLinearRegression(Cleaned_Kidney_df,"PLEKHB1")

summary(model_5_FS)
# linear regression on each gene with feature selection
model_1_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_Kidney_df,Cleaned_Kidney_df$NAV2,"NAV2")
model_2_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_Kidney_df,Cleaned_Kidney_df$AEN,"AEN")
model_3_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_Kidney_df,Cleaned_Kidney_df$PFKFB2,"PFKFB2")
model_4_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_Kidney_df,Cleaned_Kidney_df$BSPRY,"BSPRY")
model_5_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_Kidney_df,Cleaned_Kidney_df$PLEKHB1,"PLEKHB1")


# Plot kidney Data R2 and R-Adj
gene_names_list <- c("NAV2","AEN","PFKFB2","BSPRY","PLEKHB1")
MLR_FS_models <- list(model_1_FS,model_2_FS,model_3_FS,model_4_FS,model_5_FS)
MLR_models <- list(model_1,model_2,model_3,model_4,model_5)

MLR_result <- new.Prepare_For_plot(gene_names_list,MLR_models)
MLR__FS_result <- new.Prepare_For_plot(gene_names_list,MLR_FS_models) 


## MLR Plot ##
windows(width = 8, height = 5)
par(mar=c(5, 5, 1, 10), xpd=TRUE)
colors = c("blue", "yellow")
b <- barplot(height = as.matrix(MLR_result), xlab = "Genes",ylab = "Value"                      # Grouped barplot using Base R
             ,main = "Model's Performance",beside = TRUE,col = colors,ylim=c(0,0.7),)
legend("topright", inset=c(-0.2, 0.2), 
       legend=c("R2", "R-Adj"), pch=c(19),col = c("blue","yellow"),cex=0.7)


## MLR_FS plot ##
windows(width = 8, height = 5)
par(mar=c(5, 5, 1, 10), xpd=TRUE)
colors = c("blue", "yellow")
b <- barplot(height = as.matrix(MLR__FS_result), xlab = "Genes",ylab = "Value"                      # Grouped barplot using Base R
             ,main = "Model's Performance",beside = TRUE,col = colors,ylim=c(0,0.7),)
legend("topright", inset=c(-0.2, 0.2), 
       legend=c("R2", "R-Adj"), pch=c(19),col = c("blue","yellow"),cex=0.7)





##### lung Data #####

lung_CNA <- read.delim("Project_Data/lusc_CNV_core.txt", header = TRUE, sep = "\t", dec = ".")
Significant_lung <- read_csv("Project_Data/sig_genes_lung.csv") # Significant Genes 
lung_gene_exp <- read.csv("Project_Data/Lung_Cancer.txt", sep="")


# prepare lung Data
Final_compined_lung_df <- new.Prepare_data(lung_gene_exp , lung_CNA , Significant_lung)

# Drop columns with 50% zeros or NAN
Cleaned_lung_df <- new.Drop_zeros_NA(Final_compined_lung_df,0.50)

# linear regression on each gene without feature selection
model_6 <- new.MultipleLinearRegression(Cleaned_lung_df,"TNFAIP8L2",type = "lung")
model_7 <- new.MultipleLinearRegression(Cleaned_lung_df,"PFDN2",type = "lung")
model_8 <- new.MultipleLinearRegression(Cleaned_lung_df,"GINS1",type = "lung")
model_9 <- new.MultipleLinearRegression(Cleaned_lung_df,"GPX3",type = "lung")
model_10 <- new.MultipleLinearRegression(Cleaned_lung_df,"LRP2BP",type ="lung")


# linear regression on each gene with feature selection
model_6_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_lung_df,Cleaned_lung_df$TNFAIP8L2,"TNFAIP8L2",type ="lung")
model_7_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_lung_df,Cleaned_lung_df$PFDN2,"PFDN2",type ="lung")
model_8_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_lung_df,Cleaned_lung_df$GINS1,"GINS1",type ="lung")
model_9_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_lung_df,Cleaned_lung_df$GPX3,"GPX3",type ="lung")
model_10_FS <- new.MultipleLinearRegression_WithFeatureSelection(Cleaned_lung_df,Cleaned_lung_df$LRP2BP,"LRP2BP",type ="lung")



## plot ## 

# Plot Lung Data R2 and R-Adj
gene_names_list_lung <- c("TNFAIP8L2","PFDN2","GINS1","GPX3","LRP2BP")
MLR_FS_models_lung <- list(model_6_FS,model_7_FS,model_8_FS,model_9_FS,model_10_FS)
MLR_models_lung <- list(model_6,model_7,model_8,model_9,model_10)

MLR_result_lung <- new.Prepare_For_plot(gene_names_list_lung,MLR_models_lung)
MLR__FS_result_lung <- new.Prepare_For_plot(gene_names_list_lung,MLR_FS_models_lung) 



## MLR Plot ##
windows(width = 8, height = 5)
par(mar=c(5, 5, 1, 10), xpd=TRUE)
colors = c("blue", "yellow")
b <- barplot(height = as.matrix(MLR_result_lung), xlab = "Genes",ylab = "Value"                      # Grouped barplot using Base R
             ,main = "Model's Performance MLR Lung",beside = TRUE,col = colors,ylim=c(0,1.5),)
legend("topright", inset=c(-0.2, 0.2), 
       legend=c("R2", "R-Adj"), pch=c(19),col = c("blue","yellow"),cex=0.7)


## MLR_FS plot ##
windows(width = 8, height = 5)
par(mar=c(5, 5, 1, 10), xpd=TRUE)
colors = c("blue", "yellow")
b <- barplot(height = as.matrix(MLR__FS_result_lung), xlab = "Genes",ylab = "Value"                      # Grouped barplot using Base R
             ,main = "Model's Performance MLR FS Lung",beside = TRUE,col = colors,ylim=c(0,1.5),)
legend("topright", inset=c(-0.2, 0.2), 
       legend=c("R2", "R-Adj"), pch=c(19),col = c("blue","yellow"),cex=0.7)


