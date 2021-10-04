## Mastherthesis, Melanie Fattohi
## calculate p-value of deconvolution with SCDC using Monte Carlo simulation, inspired by CIBERSORT
## input for deconvolution 
### a) bulk RNA-seq dataset (numeric matrix, samples in columns, genes in rows)
### b) metadata for bulk (matrix, metadata in columns, samples in rows)
### c) path to Expression Set object(s) of qc'ed scRNA-seq dataset(s) (vector)
### d) cell types (vector)
### e) Ensemble true/false
### f) multiple_donors true/false (maybe get the info with c))
## input for Monte Carlo simulation
### g (?)
## one function that samples a) and executes Deconvolve_SCDC()
## one function that executes Deconvolve_SCDC() of whole a) and then performs Monte Carlo simulation
## for calculation of p-value
## output:
### SCDC object with predicted cell type proportions and reconstruction error and p-value

library(SCDC)
#library(MonteCarlo)
library(parallel)

## write utility function that calcs the ensemble prediction with highest mean weight
get_ensemble_res <- function(ensemble_output){
  weights <- w_table[1:5, 1:(ncol(ensemble_output$w_table)-4)]
  ensemble_res <- ensemble_output$prop.list[[which.max(colMeans(weights))]]
  return(ensemble_res)
}

get_spearman_vec <- function(decon_res){
  if(is.null(decon_res$ensemble_res)){
    decon_res <- decon_res$scdc_res
  } else if(is.null(decon_res$scdc_res)){
    decon_res <- get_ensemble_res(decon_res$ensemble_res)
  }
  spearman_vec <- as.vector(decon_res$yeval$spearmany.sample.table)
  return(spearman_vec)
}

get_pearson_vec <- function(bulk_data, decon_res){
  if(is.null(decon_res$ensemble_res)){
    decon_res <- decon_res$scdc_res
  } else if(is.null(decon_res$scdc_res)){
    decon_res <- get_ensemble_res(decon_res$ensemble_res)
  }
  bulk_data_est <- decon_res$basis.mvw  %*% t(decon_res$prop.est.mvw)
  pearson_vec <- sapply(1:ncol(bulk_data), function(x) cor(bulk_data[,x], bulk_data_est[,x]))
  return(pearson_vec)
}

## g = #genes that are present in both bulk and scrna
get_statistic_decon_sampled_bulk <- function(g, bulk_data, ...){
  ## generate sampled bulk_data
  sampled_genes <- sample(rownames(bulk_data), g)
  sampled_bulk <- bulk_data[match(sampled_genes, rownames(bulk_data)),]
  
  ## perform decon of sampled bulk_data with Deconvolve_SCDC
  message(paste("Number of genes in sampled bulk dataset: ", g, sep = ""))
  decon_res <- Deconvolve_SCDC(bulk_data = sampled_bulk, ...)
  #decon_res <- decon_res[[which(!sapply(decon_res, is.null))]]
  
  spearman_vec <- get_spearman_vec(decon_res = decon_res)
  pearson_vec <- get_pearson_vec(bulk_data = sampled_bulk, decon_res = decon_res)
  
  statistics <- list(spearman_vec, pearson_vec)
  ## output: Spearman correlation
  return(statistics)
}


Calculate_pvalue <- function(g = c(5000, 10000, 15000), nrep = 500, ncores = 5, silent = TRUE, bulk_data, ...) { 
  
  ## Deconvolution of whole bulk RNA-seq dataset
  message("Executing deconvolution with the whole bulk RNA-seq dataset ..")
  decon_res <- Deconvolve_SCDC(bulk_data, ...)
  spearman_vec_whole <- get_spearman_vec(decon_res = decon_res)
  pearson_vec_whole <- get_pearson_vec(bulk_data = bulk_data, decon_res = decon_res)
  
  ## Deconvolution of sampled bulk RNA-seq dataset
  message("Calculation of p-value. This step takes some time ..")
  spearman_matrix_sampled <- 
    do.call(cbind, mclapply(1:nrep, function(x) sapply(g, 
                                                       function(y) get_statistic_decon_sampled_bulk(y, 
                                                                                                   bulk_data = bulk_data, ...)), 
                            mc.cores = ncores, mc.silent = silent))
  spearman_matrix_sampled <- apply(spearman_matrix_sampled, MARGIN= 2, FUN = sort)
  #param_list <- list("g" = g, "bulk_data" = bulk_data)
  #MC_result <- MonteCarlo(func=sample_bulk, nrep = 500, param_list = param_list, ncpus = 5) 

  ## calculating p-value
  p_value <- apply(spearman_matrix_sampled, MARGIN = 2, function(x) {
    1 - (which.min(abs(x - spearman_vec_whole)) / length(spearman_matrix_sampled))})
  message("Done.")

  decon_res_pval <- list("decon_res" = decon_res, "p_value" = p_value)
  return(decon_res_pval)
  
} 
