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
library(MonteCarlo)

## write utility function that calcs the ensemble prediction with highest mean weight
get_ensemble_res <- function(ensemble_output){
  weights <- w_table[1:5, 1:(ncol(ensemble_output$w_table)-4)]
  ensemble_res <- ensemble_output$prop.list[[which.max(colMeans(weights))]]
  return(ensemble_res)
}


## g = #genes that are present in both bulk and scrna
sample_bulk <- function(g, bulk_data, ...){#g = NULL, ...){
  ## generate sampled bulk_data
  if(is.null(g)){
    message("g cannot be NULL!")
    #sc_data <- readRDS(file = sc_path)
    #g <- unname(nrow(sc_data$sc.eset.qc@featureData))
  }

  sampled_genes <- sample(rownames(bulk_data), g)
  sampled_bulk <- bulk_data[match(sampled_genes, rownames(bulk_data)),]
  
  ## perform decon of sampled bulk_data with Deconvolve_SCDC
  deco_res <- Deconvolve_SCDC(bulk_data = sampled_bulk, ...)
  
  if(is.null(deco_res$ensemble_res)){
     
    spearman_vec <- as.vector(deco_res$scdc_res$yeval$spearmany.sample.table)
    
  } else if(is.null(deco_res$scdc_res)){
    
    deco_res_ensemble <- get_ensemble_res(deco_res$ensemble_res)
    spearman_vec <- as.vector(deco_res_ensemble$yeval$spearmany.sample.table)
  }
  ## output: Spearman correlation
  return(list("spearman_vec" = spearman_vec))
}



Calculate_pvalue <- function(g = NULL, bulk_data, ...) { 
  ## g is a vector of integers
  
  message("Executing deconvolution with the whole bulk RNA-seq dataset ..")
  deco_res <- Deconvolve_SCDC(bulk_data, bulk_meta, sc_path, cell_types, ensemble, multiple_donors)
  
  if(is.null(g)){
    message("g cannot be NULL!")
    param_list = list("g" = c(5000, 10000, 15000))#, bulk_data = bulk_data) 
  } else {
    param_list = list("g" = g, bulk_data = bulk_data) 
  }
  
  message("Create Monte Carlo simulation for calculation of p-value. This step takes some time ..")
  MC_result <- MonteCarlo(func=sample_bulk, nrep = 10, param_list = param_list, ncpus = 5) 
  
  ## get the R's from MC_result (nulldist) vector/matrix
  ## get the R's from deco_res (mix_r) vector
  ## pval = pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
  
  ## return pval and deco_res
  
} 