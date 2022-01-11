## Mastherthesis, Melanie Fattohi
## correlation analysis 
## test if predicted cell type proportions are correlated with clinically relevant characteristics;
## also: visualization using heatmaps and box plots (boxplots for p-values)

#factor_to_numeric <-


correlation_analysis <- function(decon_output, cell_types = NULL, clinical_characteristic){
  # clinical characteristic has to be of type either numeric or factor
  
  if(!(is.factor(clinical_characteristic) || is.numeric(clinical_characteristic))){
    stop("clinical_characteristic has to be of type numeric or factor")
  } else if(is.factor(clinical_characteristic)){
    ## need to find out how to compute correlation if clinical characteristic is categorical (i.e. factor)
  }
  
  if(!is.data.frame(clinical_characteristics)){
    #clinical_characteristics <- as.data.frame(clinical_characteristics)
    stop("clinical_characteristics has to be of type data.frame")
  }
  
  ct_prop <- decon_output$decon_res$prop.est.mvw
  if(!is.null(cell_types)){
    ct_prop <- as.data.frame(ct_prop[,cell_types])
    colnames(ct_prop) <- cell_types
  }
  
  #return the correlations, the p-values and the two plots
  
}