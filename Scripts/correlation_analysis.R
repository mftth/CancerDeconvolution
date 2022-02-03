## Mastherthesis, Melanie Fattohi
## correlation analysis 
## test if predicted cell type proportions are correlated with clinically relevant characteristics;

correlation_analysis <- function(decon_output, cell_types = NULL, clinical_characteristic){
  ct_prop <- decon_output$decon_res$prop.est.mvw
  if(!is.null(cell_types)){
    ct_prop <- as.data.frame(ct_prop[,cell_types])
    colnames(ct_prop) <- cell_types
  }
  
  # clinical characteristic has to be of type either numeric or character
  if(!(is.character(clinical_characteristic) || is.numeric(clinical_characteristic))){
    stop("clinical_characteristic has to be of type numeric or character")
    
  } else if(is.character(clinical_characteristic)){
    aov_data <- as.data.frame(cbind(ct_prop, clinical_characteristic))
    aov_call <- lapply(1:(ncol(aov_data)-1), 
                       function(x) aov(as.formula(paste0(colnames(aov_data)[x], 
                                                         " ~ clinical_characteristic", collapse = "")),
                                       data = aov_data))
    names(aov_call) <- colnames(aov_data)[1:(ncol(aov_data)-1)]
    aov_summary <- lapply(aov_call, function(x) summary(x))
    aov_pvalue <- sapply(aov_summary, function(x) x[[1]]$`Pr(>F)`[1])
    return(list("aov_pvalue" = aov_pvalue))
      
  } else if(is.numeric(clinical_characteristic)){
    cor_ct_prop <- as.numeric(cor(ct_prop, clinical_characteristic))
    cor_ct_prop_pval <- sapply(1:ncol(ct_prop), function(x) cor.test(ct_prop[,x],
                                                                     clinical_characteristic)$p.value)
    names(cor_ct_prop) <- colnames(ct_prop)
    names(cor_ct_prop_pval) <- colnames(ct_prop)
    return(list("corr" = cor_ct_prop,
                "corr_pvalue" = cor_ct_prop_pval))
  }
  
}
