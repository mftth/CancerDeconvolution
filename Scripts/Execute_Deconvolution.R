## Mastherthesis, Melanie Fattohi
## execute deconvolution with SCDC
## input is same input as for Quality_Control(), Deconvolve_SCDC() and/or Calculate_pvalue()
## need decision about if QC should take place or no and if p-value should be calculated or not 
## ouput is new SCDC "model", deconvolution result and possibly p-value of deconvolution

Execute_Deconvolution <- function(qc, calc_pval, ...){
  if(qc){
    sc_qc <- Quality_Control(...)
    if(calc_pval){
      decon_res <- Calculate_pvalue(...)
    } else {
      decon_res <- Deconvolve_SCDC(...)
    }
  } else {
    sc_qc <- NULL
    if(calc_pval){
      decon_res <- Calculate_pvalue(...)
    } else {
      decon_res <- Deconvolve_SCDC(...)
    }
  }
  
  decon_output <- list("sc_qc" = sc_qc, "decon_res" = decon_res)
  return(decon_output)
}