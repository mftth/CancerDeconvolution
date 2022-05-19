## Mastherthesis, Melanie Fattohi

source("~/Masterthesis/CancerDeconvolution/Scripts/SCDC_Basic_Functions.R")
#source("~/SCDC/SCDC/R/Deconvolution.R")


prepare_decon_res <- function(p_value, decon_res, clinical_char, observed_statistics = NULL){
  ## p_value : boolean. True, if p-value was estimated; False, if p-value was not estimated
  ## decon_res: output of Calculate_pvalue
  ## clinical_char: numeric or character vector of one clinical variable
  ## observed_statistics : output of get_test_statistics_vec. only necessary, if p-value = FALSE 
  
  if (p_value){
    decon_res_prepped <- as.data.frame(cbind(decon_res$decon_res$prop.est.mvw, decon_res$statistics_observed$pearson_vec, 
                                             decon_res$statistics_observed$spearman_vec, decon_res$statistics_observed$mad_vec,
                                             decon_res$statistics_observed$rmsd_vec, decon_res$p_value_per_sample, clinical_char))
    colnames(decon_res_prepped) <- c(colnames(decon_res$decon_res$prop.est.mvw), "pearson", "spearman", "mad", "rmsd", 
                                     "pearson_pval", "spearman_pval", "mad_pval", "rmsd_pval", "response")
    decon_res_prepped[,1:(ncol(decon_res$decon_res$prop.est.mvw)+8)] <- sapply(1:(ncol(decon_res$decon_res$prop.est.mvw)+8), 
                                                                               function(x) as.numeric(decon_res_prepped[,x]))
    decon_res_prepped$response <- as.factor(decon_res_prepped$response)
  } else {
    decon_res_prepped <- as.data.frame(cbind(decon_res$prop.est.mvw, observed_statistics$pearson_vec, 
                                             observed_statistics$spearman_vec, observed_statistics$mad_vec,
                                             observed_statistics$rmsd_vec, clinical_char))
    colnames(decon_res_prepped) <- c(colnames(decon_res$prop.est.mvw), "pearson", "spearman", "mad", "rmsd", "response")
    decon_res_prepped[,1:(ncol(decon_res$prop.est.mvw)+4)] <- sapply(1:(ncol(decon_res$prop.est.mvw)+4), 
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