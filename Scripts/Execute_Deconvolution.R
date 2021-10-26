## Mastherthesis, Melanie Fattohi
## execute deconvolution with SCDC
## input is same input as for Quality_Control(), Deconvolve_SCDC() and/or Calculate_pvalue()
## need decision about if QC should take place or no and if p-value should be calculated or not 
## ouput is new SCDC "model", deconvolution result and possibly p-value of deconvolution

Execute_Deconvolution <- function(qc, calc_pval, ...){
  if(qc){
    sc_qc <- Quality_Control(...)
    if(calc_pval){
      decon_res <- Calculate_pvalue(...) # 7 features + cell types
    } else {
      decon_res <- Deconvolve_SCDC(...) # 3 features + cell types
    }
  } else {
    sc_qc <- NULL
    if(calc_pval){
      decon_res <- Calculate_pvalue(...) # 7 features + cell types
    } else {
      decon_res <- Deconvolve_SCDC(...) # 3 features + cell types
    }
  }
  
  decon_output <- list("sc_qc" = sc_qc, "decon_res" = decon_res)
  return(decon_output)
}

prepare_decon_res <- function(decon_res, clinical_char, p_value){
  ## check if pval calc took place
  if (p_value){
    decon_res_prepped <- as.data.frame(cbind(decon_res$decon_res$prop.est.mvw, t(Repset_scdc_baron$decon_res$yeval$spearmany.sample.table), 
                                             t(Repset_scdc_baron$decon_res$yeval$mADy.sample.table), t(Repset_scdc_baron$decon_res$yeval$RMSDy.sample.table),
                                             Repset_scdc_baron$p_value_per_sample, repset_meta$Grading))
  } else {
    
  }
  ## form a dataframe out of the deco output --> maybe do this already in Execute_Deconvolution
  decon_res <- as.data.frame(cbind(Repset_scdc_baron$decon_res$prop.est.mvw, t(Repset_scdc_baron$decon_res$yeval$spearmany.sample.table), 
                                   t(Repset_scdc_baron$decon_res$yeval$mADy.sample.table), t(Repset_scdc_baron$decon_res$yeval$RMSDy.sample.table),
                                   Repset_scdc_baron$p_value_per_sample, repset_meta$Grading))
  colnames(decon_res) <- c(colnames(Repset_scdc_baron$decon_res$prop.est.mvw), "spearman", "mad", "rmsd",
                           "pearson_pval", "spearman_pval", "mad_pval", "rmsd_pval", "response") ## p-value
  decon_res[,1:(ncol(Repset_scdc_baron$decon_res$prop.est.mvw)+7)] <- sapply(1:(ncol(Repset_scdc_baron$decon_res$prop.est.mvw)+7), 
                                                                             function(x) as.numeric(decon_res[,x]))
}