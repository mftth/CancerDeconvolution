## Mastherthesis, Melanie Fattohi
## correlation analysis 
## test if predicted cell type proportions are correlated with clinically relevant characteristics;
## also: visualization using heatmaps and box plots (boxplots for p-values)


correlation_analysis <- function(decon_output, cell_types = NULL, clinical_characteristic){
  ct_prop <- decon_output$decon_res$prop.est.mvw
  if(!is.null(cell_types)){
    ct_prop <- as.data.frame(ct_prop[,cell_types])
    colnames(ct_prop) <- cell_types
  }
  
  # clinical characteristic has to be of type either numeric or factor
  if(!(is.factor(clinical_characteristic) || is.numeric(clinical_characteristic))){
    stop("clinical_characteristic has to be of type numeric or factor")
    
  } else if(is.factor(clinical_characteristic)){
    # for each ct proportion one anova
    aov_data <- cbind(ct_prop, clinical_characteristic)
    aov(as.formula(paste0("clinical_characteristic ~ ", paste0(colnames(ct_prop), collapse = " + "), collape = "")),
        data = aov_data)
    
  } else if(is.numeric(clinical_characteristic)){
      cor_ct_prop <- as.numeric(cor(ct_prop, clinical_characteristic))
      cor_ct_prop_pval <- sapply(1:ncol(ct_prop), function(x) cor.test(ct_prop[,x],
                                                                       clinical_characteristic)$p.value)
      names(cor_ct_prop) <- colnames(ct_prop)
      names(cor_ct_prop_pval) <- colnames(ct_prop)
      return(list("corr" = cor_ct_prop,
                  "p-values" = cor_ct_prop_pval))
  }
  
  #return the correlations, the p-values and the two plots
}