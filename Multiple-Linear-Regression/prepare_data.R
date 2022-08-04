# Function to prepare the data for the Linear Regression algorithm #
new.Prepare_data <- function(Genes_DF , CNA_DF , Sig_Genes_DF) {
  
  # replace "." with "-" in column names
  colnames(Genes_DF) <-  gsub("\\.", "-", colnames(Genes_DF)) 
  
  # pick the 5 most significant Genes
  diff_Expressed_genes <- Sig_Genes_DF[1:5,]$Gene_name 
  
  # get only the rows matching with the 5 significant genes 
  Five_Genes_df <- Genes_DF[Genes_DF$Hugo_Symbol %in% diff_Expressed_genes,] 
  
  # intersection between the Genes DF and the CNA_DF
  intersection <- CNA_DF[CNA_DF$feature %in% colnames(Genes_DF),]
  
  # Get patients names
  patients <- intersection$feature
  
  # for loop to add 5 new columns to our intersection df representing the 5 genes
  for (gene in Five_Genes_df$Hugo_Symbol)
  {
    
    Genes_vector <- rep() # Empty vector to append the gene value for each patient
    
    for (patient in patients)
    {
      # for each patient select the patient and the genes values
      patient_Genes_df <- subset(Five_Genes_df, select=c("Hugo_Symbol", patient))
      
      # get the gene level
      value <- patient_Genes_df[patient_Genes_df$Hugo_Symbol == gene ,2]
      
      # append the gene level value in the vector
      Genes_vector <- append(Genes_vector,value)
    }
    
    # After for loops add column representing the gene and its values for each patient
    intersection[gene] <- Genes_vector
  }
  
  # drop Empty column "X"
  intersection <- subset(intersection, select = -c(X) )
  
  return(intersection)
  
}


# Function to drop columns with 50% zeros or NAN
new.Drop_zeros_NA <- function(Dataframe ,thresh = 0.50)
{ 
  # Drop columns with 0 or NA more than or equal 50% of the column length
  Cleaned <- Dataframe[, sapply(Dataframe, function(x) length(x[!(is.na(x) | x == 0)]) / length(x) >= thresh)];
  return(Cleaned)
}

