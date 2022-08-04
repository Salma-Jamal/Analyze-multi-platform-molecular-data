# prepare data for plot
new.Prepare_For_plot <- function(genes_list, models_list)
{
  
  R_Squared_vector <- rep()
  R_Adja_vector <- rep()
  
  for ( i in models_list )
  {
    model_summary <- summary(i)
    R_Squared_vector <- append(R_Squared_vector,model_summary$r.squared)
    R_Adja_vector <- append(R_Adja_vector,model_summary$adj.r.squared)
  }
  
  data <- data.frame(genes_list,R_Squared_vector, R_Adja_vector)
  data_trans <- as.data.frame(t(data))
  colnames(data_trans) <- data_trans[1,]
  data_trans <- data_trans[-1,]
  
  data_trans[, c(1:5)] <- sapply(data_trans[, c(1:5)], as.double)
  return(data_trans)
}
