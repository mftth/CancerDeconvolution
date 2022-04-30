## Mastherthesis, Melanie Fattohi
## correlation analysis 
## test if predicted cell type proportions are correlated with clinically relevant characteristics;

correlation_analysis <- function(decon_output, cell_types = NULL, clinical_characteristic){
  ct_prop <- data.frame(decon_output$decon_res$prop.est.mvw)
  
  if(any(is.na(clinical_characteristic))){
    ct_prop <- ct_prop[- which(is.na(clinical_characteristic)),]
    clinical_characteristic <- clinical_characteristic[-which(is.na(clinical_characteristic))]
  }
  if(!is.null(cell_types)){
    ct_prop <- as.data.frame(ct_prop[,cell_types])
    colnames(ct_prop) <- cell_types
  }
  
  # clinical characteristic has to be of type either numeric or character
  if(!(is.character(clinical_characteristic) || is.numeric(clinical_characteristic))){
    stop("clinical_characteristic has to be of type numeric or character")
    
  } else if(is.character(clinical_characteristic)){
    aov_data <- ct_prop
    aov_data$clinical_characteristic <- clinical_characteristic
    aov_call <- lapply(1:(ncol(aov_data)-1), 
                       function(x) aov(as.formula(paste0(colnames(aov_data)[x], 
                                                         " ~ clinical_characteristic", collapse = "")),
                                       data = aov_data))
    names(aov_call) <- colnames(aov_data)[1:(ncol(aov_data)-1)]
    aov_summary <- lapply(aov_call, function(x) summary(x))
    aov_pvalue <- sapply(aov_summary, function(x) x[[1]]$`Pr(>F)`[1])
    sign_pval <- (which(aov_pvalue < 0.05))
    comparison_list <- rje::powerSet(unique(clinical_characteristic))
    comparison_list <- comparison_list[which(sapply(comparison_list, function(x) length(x)) == 2)]
    
    aov_plots <- lapply(sign_pval, function(celltype){ 
      ggplot(aov_data, aes(x = clinical_characteristic, y = aov_data[,celltype])) +
      geom_boxplot() +
      scale_x_discrete() + xlab("clinical_characteristic") +
      ylab(paste0(colnames(aov_data)[celltype], " proportions", collapse = "")) +
      geom_signif(comparisons = comparison_list, map_signif_level=TRUE) +
      theme_bw() + theme(text = element_text(size = 11))  
    })
    
    violin_plots <- lapply(sign_pval, function(celltype){ 
      ggplot(aov_data, aes(x = clinical_characteristic, y = aov_data[,celltype])) +
        geom_violin() +
        #scale_x_discrete() + 
        xlab("clinical_characteristic") +
        ylab(paste0(colnames(aov_data)[celltype], " proportions", collapse = "")) +
        geom_signif(comparisons = comparison_list, map_signif_level=TRUE) +
        theme_bw() + theme(text = element_text(size = 11))  +
        stat_summary(fun=median, geom="point", size=2, color="red") +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)
    })
    
    return(list("aov_pvalue" = aov_pvalue, "aov_plots" = aov_plots, "violin_plots" = violin_plots))
      
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
