## Mastherthesis, Melanie Fattohi
## survival analysis 
## test if relevant cell type proportions are predictive in regards to disease-related survival;
## relevant cell types are determined through Machine Learning (feature importance) or otherwise
## just testing cell types of interest
## also: visualization using Kaplan-Meier plots

library(survival)
library(survminer)
library(ggplot2)


continuous_to_discrete <- function(con_vec, col_name){
  dis_vec <- rep(NA, length(con_vec))
  dis_vec[which(con_vec <= median(con_vec))] <- paste0(col_name, "_low", collapse = "")
  dis_vec[which(con_vec > median(con_vec))] <- paste0(col_name, "_high", collapse = "")
  return(dis_vec)
}

survival_analysis <- function(decon_output, cell_types = NULL, OS, censor, clinical_characteristics, ...){
  ## cell_types is vector of cell types of interest (i.e. those investigated in survival analysis)
  ## OS is numeric, censor integer vector
  ## visualize only those with p-value < 0.05
  
  if(!is.data.frame(clinical_characteristics)){
    #clinical_characteristics <- as.data.frame(clinical_characteristics)
    stop("clinical_characteristics has to be of type data.frame")
  }
  
  ct_prop <- decon_output$decon_res$prop.est.mvw
  if(!is.null(cell_types)){
    ct_prop <- as.data.frame(ct_prop[,cell_types])
    colnames(ct_prop) <- cell_types
  }
  ## transform numerical cell type proportion into categorical variable
  ct_prop_categories <- sapply(1:ncol(ct_prop), 
                               function(x) continuous_to_discrete(ct_prop[,x], colnames(ct_prop)[x]))
  rownames(ct_prop_categories) <- rownames(ct_prop)
  colnames(ct_prop_categories) <- sapply(colnames(ct_prop), 
                                         function(x) paste0(x, "_high_low", collapse = ""))
  
  ## transform numerical characteristics into categorical variable
  num_characteristics <- which(sapply(1:ncol(clinical_characteristics), 
                                      function(x) class(clinical_characteristics[,x])) == "numeric")
  
  
  if(length(num_characteristics)==0){
    clinical_characteristics_cat <- clinical_characteristics
  } else {
    clinical_characteristics_cat <- sapply(num_characteristics, 
                                           function(x) continuous_to_discrete(clinical_characteristics[,x],
                                                                              colnames(clinical_characteristics)[x]))
    rownames(clinical_characteristics_cat) <- rownames(clinical_characteristics)
    colnames(clinical_characteristics_cat) <- sapply(colnames(clinical_characteristics)[num_characteristics], 
                                                     function(x) paste0(x, "_high_low", collapse = ""))
    clinical_characteristics_cat <- cbind(clinical_characteristics[,-num_characteristics],
                                          clinical_characteristics_cat)
  }
  
  ## fit survfit object, create survival curve and manipulate formula
  survival_meta <- as.data.frame(cbind(ct_prop_categories, clinical_characteristics_cat))
  survival_fit <- lapply(1:ncol(survival_meta), 
                         function(x) survfit(as.formula(paste0("Surv(OS, censor) ~ ", colnames(survival_meta)[x], 
                                                               collapse = "")), 
                                             data = survival_meta))
  names(survival_fit) <- colnames(survival_meta)
  for (idx in 1:length(survival_fit)) {
    survival_fit[[idx]]$call$formula <- as.formula(paste0("Surv(OS, censor) ~ ", colnames(survival_meta)[idx], 
                                                          collapse = ""))
  }
  
  ## calculate pvalue from survfit object
  survival_pval <- sapply(survival_fit, 
                          function(x) survminer::surv_pvalue(x, data = survival_meta)$pval)
  
  ## create Kaplan-Meier plots
  kp_plot <- ggsurvplot(survival_fit, data = survival_meta, combine = TRUE, ...)
  single_kp_plots <- lapply(survival_fit, 
                            function(x) ggsurvplot(x, data = survival_meta, pval = TRUE,
                                                   legend.title = ""))
  
  return(list("survfit_objects" = survival_fit,
              "survival_pvals" = survival_pval,
              "combined_kp" = kp_plot,
              "single_kp" = single_kp_plots))
}
