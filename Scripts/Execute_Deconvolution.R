## Mastherthesis, Melanie Fattohi
## execute deconvolution with SCDC
## input is same input as for Quality_Control(), Deconvolve_SCDC() and/or Calculate_pvalue()
## need decision about if QC should take place or no and if p-value should be calculated or not 
## ouput is new SCDC "model", deconvolution result and possibly p-value of deconvolution

source("~/SCDC/SCDC/R/Basic_Functions.R")
source("~/SCDC/SCDC/R/Deconvolution.R")


prepare_decon_res <- function(decon_res, clinical_char, p_value){
  if (p_value){
    decon_res_prepped <- as.data.frame(cbind(decon_res$decon_res$prop.est.mvw, t(decon_res$decon_res$yeval$spearmany.sample.table), 
                                             t(decon_res$decon_res$yeval$mADy.sample.table), t(decon_res$decon_res$yeval$RMSDy.sample.table),
                                             decon_res$p_value_per_sample, clinical_char))
    colnames(decon_res_prepped) <- c(colnames(decon_res$decon_res$prop.est.mvw), "spearman", "mad", "rmsd", 
                                     "pearson_pval", "spearman_pval", "mad_pval", "rmsd_pval", "response")
    decon_res_prepped[,1:(ncol(decon_res$decon_res$prop.est.mvw)+7)] <- sapply(1:(ncol(decon_res$decon_res$prop.est.mvw)+7), 
                                                                               function(x) as.numeric(decon_res_prepped[,x]))
    decon_res_prepped$response <- as.factor(decon_res_prepped$response)
  } else {
    decon_res_prepped <- as.data.frame(cbind(decon_res$prop.est.mvw, t(decon_res$yeval$spearmany.sample.table), 
                                             t(decon_res$yeval$mADy.sample.table), t(decon_res$yeval$RMSDy.sample.table),
                                             clinical_char))
    colnames(decon_res_prepped) <- c(colnames(decon_res$prop.est.mvw), "spearman", "mad", "rmsd", "response")
    decon_res_prepped[,1:(ncol(decon_res$prop.est.mvw)+3)] <- sapply(1:(ncol(decon_res$prop.est.mvw)+3), 
                                                                     function(x) as.numeric(decon_res_prepped[,x]))
    decon_res_prepped$response <- as.factor(decon_res_prepped$response)
  }
  
  return(decon_res_prepped)
}


# Execute_Deconvolution <- function(qc, calc_pval, sc_data = NULL, sc_meta, sc_path, multiple_donors, cell_types, clinical_char, ...){
#   if(qc){
#     sc_qc <- Quality_Control(sc_data = sc_data, sc_meta = sc_meta, sc_path = sc_path, 
#                              multiple_donors = multiple_donors, cell_types = cell_types)
#     if(calc_pval){
#       decon_res <- Calculate_pvalue(sc_data = sc_qc$sc.eset.qc, ...)
#       decon_res_prepped <- prepare_decon_res(decon_res = decon_res, clinical_char = clinical_char, 
#                                              p_value = TRUE)
#     } else {
#       sc_basis <- Create_basis(sc_data = sc_qc$sc.eset.qc)
#       decon_res <- Deconvolve_SCDC(sc_data = sc_data, multiple_donors = multiple_donors, cell_types = cell_types, ...)
#       decon_res_prepped <- prepare_decon_res(decon_res = decon_res, clinical_char = clinical_char, 
#                                              p_value = FALSE)
#     }
#   } else {
#     sc_qc <- NULL
#     if(calc_pval){
#       decon_res <- Calculate_pvalue(...)
#       decon_res_prepped <- prepare_decon_res(decon_res = decon_res, clinical_char = clinical_char, 
#                                              p_value = TRUE)
#     } else {
#       decon_res <- Deconvolve_SCDC(...)
#       decon_res_prepped <- prepare_decon_res(decon_res = decon_res, clinical_char = clinical_char, 
#                                              p_value = FALSE)
#     }
#   }
# 
#   decon_output <- list("sc_qc" = sc_qc, "decon_res" = decon_res)
#   return(decon_output)
# }