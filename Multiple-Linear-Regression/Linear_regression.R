library(leaps) # for feature selection

# Function to perform Multiple linear regression
new.MultipleLinearRegression <- function(Genes_Dataframe, gene,type="kid")
{
      # Extract column names from data frame
      col_names = colnames(Genes_Dataframe)
      
      # string to be the input to lm
      model_variables <-  paste0(gene," ~ ")
      
      for ( i in col_names)
      {
        if(startsWith(i, "CNV"))
          {
            model_variables <- paste0(model_variables,paste0(" ",i , " +"))
          }
      }
      
      # remove the last two characters from the string + , space; 
      # start from char 1 till string length-2
      model_variables <- substring(model_variables, 1, nchar(model_variables)-2)
      model <- lm(model_variables, data = Genes_Dataframe )
      
      # significant copy number predictors for gene.
      coefficients_df <- as.data.frame(summary(model)$coefficients)
      colnames(coefficients_df) <- c("Estimate","StdError","tvalue","pvalue")  
      significant_CNP_list <- rownames(coefficients_df[coefficients_df$pvalue <= 0.05,])
      # convert the list to a string; -1 to remove intercept; width= maximum string value
      significant_CNP <- toString(significant_CNP_list[-1], width = NULL, collapse=", ")
      
    if(type == "kid")
    {
      if (length(significant_CNP_list) > 1)
        {
           write(paste0("Gene name: ",gene,"\n",significant_CNP,"\n\n","#########","\n"),"significant_CNP_MLR.txt",append=TRUE)
        }
      else
        {
           write(paste0("Gene name: ",gene,"\n","No significant copy number predictors","\n\n","#########","\n"),"significant_CNP_MLR.txt",append=TRUE)
        }
    }
      
  else
  {
    if (length(significant_CNP_list) > 1)
      {
         write(paste0("Gene name: ",gene,"\n",significant_CNP,"\n\n","#########","\n"),"significant_CNP_MLR_lung.txt",append=TRUE)
      }
    else
      {
        write(paste0("Gene name: ",gene,"\n","No significant copy number predictors","\n\n","#########","\n"),"significant_CNP_MLR_lung.txt",append=TRUE)
      }
  }
  return(model)
}


########################################################################################

# Function to perform multiple linear regression with feature selection
new.MultipleLinearRegression_WithFeatureSelection <- function(Genes_Dataframe, gene_col, gene_name,type="kid")
{
  # Extract only the genes columns from the data frame
  x_cols <- Genes_Dataframe[,2:(length(Genes_Dataframe)-5)]
  
   
  if (type == "kid")
  {
      # regsubset to determine the most influential predictors for our model
      Best_Subset <- regsubsets(x= x_cols,
                   y= gene_col,
                   nbest = 1,      # 1 best model for each number of predictors
                   nvmax = NULL,   # NULL for no limit on number of features
                   force.in = NULL, force.out = NULL, # force features in or out
                   method = "exhaustive",)
  }
  
   
  else
  {
      # regsubset to determine the most influential predictors for our model
      Best_Subset <- regsubsets(x= x_cols,
                                y= gene_col,
                                nbest = 1,      # 1 best model for each number of predictors
                                nvmax = NULL,    # NULL for no limit on number of variables
                                force.in = NULL, force.out = NULL,
                                method = "forward",
                                )
  }
  
  # summary of the regsubset
  summary_best_subset <- summary(Best_Subset)
  
  # create data frame that contains the combinations for each subset of genes 
  result_mat <- as.data.frame(summary_best_subset$outmat)
  
  # number of predictors to use for our dataset
  best_combination <- which.max(summary_best_subset$adjr2)
  
  # Extract the row with best combinations of genes
  row_with_best_comb <- result_mat[best_combination,]
  
  # get only the genes names from the row
  best_genes <- rownames(as.data.frame(apply(row_with_best_comb,1, function(x) which(x == "*"))))
  
  
  # variable to input to the model
  model_variables <-  paste0(gene_name," ~ ")
  
  # paste the genes names to the expression
  for ( i in best_genes)
  {
      model_variables <- paste0(model_variables,paste0(" ",i , " +"))
  }
  
  # remove the last two characters from the string + , space
  model_variables <- substring(model_variables, 1, nchar(model_variables)-2)
  
  model <- lm(model_variables, data = Genes_Dataframe )
  
  # significant copy number predictors for gene.
  coefficients_df <- as.data.frame(summary(model)$coefficients)
  colnames(coefficients_df) <- c("Estimate","StdError","tvalue","pvalue")  
  significant_CNP_list <- rownames(coefficients_df[coefficients_df$pvalue <= 0.05,])
  significant_CNP <- toString(significant_CNP_list[-1], width = NULL, collapse=", ")
  
  if (type == "kid")
    {
        if (length(significant_CNP_list) > 1)
          {
            write(paste0("Gene name: ",gene_name,"\n",significant_CNP,"\n\n","#########","\n"),"significant_CNP_MLR_FS.txt",append=TRUE)
          }
        else
          {
            write(paste0("Gene name: ",gene_name,"\n","No significant copy number predictors","\n\n","#########","\n"),"significant_CNP_MLR_FS.txt",append=TRUE)
          }
   }
  
  else
  {
        if (length(significant_CNP_list) > 1)
        {
          write(paste0("Gene name: ",gene_name,"\n",significant_CNP,"\n\n","#########","\n"),"significant_CNP_MLR_FS_lung.txt",append=TRUE)
        }
       else
        {
          write(paste0("Gene name: ",gene_name,"\n","No significant copy number predictors","\n\n","#########","\n"),"significant_CNP_MLR_FS_lung.txt",append=TRUE)
        }
  }
  
  return(model)
  
}