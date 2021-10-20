## Mastherthesis, Melanie Fattohi
## calculate p-value of deconvolution with SCDC
## input for deconvolution 
### a) bulk RNA-seq dataset (numeric matrix, samples in columns, genes in rows)
### b) metadata for bulk (matrix, metadata in columns, samples in rows)
### c) path to Expression Set object(s) of qc'ed scRNA-seq dataset(s) (vector)
### d) cell types (vector)
### e) Ensemble true/false
### f) multiple_donors true/false (maybe get the info with c))
## input for p-value calculation 
### (ns_genes = # genes to be sampled from bulk_data)
### nrep = # permutations
### ncores = # cores used for calculation of p-value
### silent = should output of calculation be shown
## output:
### SCDC object with predicted cell type proportions and reconstruction error and p-value (based on Pearson and Spearman)

library(SCDC)
#library(MonteCarlo)
library(parallel)
library(robustbase)

## write utility function that calcs the ensemble prediction with highest mean weight
# get_ensemble_res <- function(ensemble_output){
#   weights <- ensemble_output$w_table[1:5, 1:(ncol(ensemble_output$w_table)-4)]
#   ensemble_res <- ensemble_output$prop.list[[which.max(colMeans(weights))]]
#   return(ensemble_res)
# }

## function to calculate correlation of bulk matrix M and M'= basis * cell type proportions
get_test_statistics_vec <- function(bulk_data, decon_res){
  #if(is.null(decon_res$ensemble_res)){
  #  decon_res <- decon_res$scdc_res
  #} else if(is.null(decon_res$scdc_res)){
  #  decon_res <- get_ensemble_res(decon_res$ensemble_res)
  #}
  
  ## matching genes of bulk_data and scRNA-seq reference
  common_genes <- intersect(rownames(bulk_data), rownames(decon_res$basis.mvw))
  bulk_data_obs <- bulk_data[common_genes,]
  #bulk_data_obs <- as.matrix(bulk_data_obs, ncol = ncol(bulk_data_obs))
  bulk_data_obs <- getCPM0(bulk_data_obs)
  bulk_data_est <- decon_res$basis.mvw[common_genes,]  %*% t(decon_res$prop.est.mvw)
  bulk_data_est <- getCPM0(bulk_data_est)
  
  ## calculate pearson and spearman correlation of bulk_data and scRNA-seq reference %*% cell type proportions
  pearson_vec <- sapply(1:ncol(bulk_data_obs), function(x) cor(bulk_data_obs[,x], bulk_data_est[,x], method = "pearson"))
  spearman_vec <- sapply(1:ncol(bulk_data_obs), function(x) cor(bulk_data_obs[,x], bulk_data_est[,x], method = "spearman"))
  mad_vec <- colMedians(abs((bulk_data_est - bulk_data_obs) - median(bulk_data_est - bulk_data_obs))) # median absolute deviation
  rmsd_vec <- sqrt(colMeans((bulk_data_est - bulk_data_obs)^2))
  
  
  test_statistics_vec <- list("pearson_vec" = pearson_vec, "spearman_vec" = spearman_vec,
                          "mad_vec" = mad_vec, "rmsd_vec" = rmsd_vec)
  return(test_statistics_vec)
}


get_statistics_decon_sampled_bulk <- function(ns_genes, bulk_data, ...){
  ## generate sampled bulk_data by generating ns_genes random rows
  sampled_idx <- sample(1:nrow(bulk_data), ns_genes)
  sampled_bulk <- bulk_data[sampled_idx,]
  ## shuffle gene labels
  #sampled_genes <- sample(rownames(sampled_bulk))
  #rownames(sampled_bulk) <- sampled_genes
  #sampled_bulk <- bulk_data[match(sampled_genes, rownames(bulk_data)),]
  
  ## perform decon of sampled bulk_data with Deconvolve_SCDC
  message(paste("Number of genes in sampled bulk dataset: ", ns_genes, sep = ""))
  decon_sampled_res <- Deconvolve_SCDC(bulk_data = sampled_bulk, ...)
  #decon_res <- decon_res[[which(!sapply(decon_res, is.null))]]

  statistics_vec <- get_test_statistics_vec(bulk_data = sampled_bulk, decon_res = decon_sampled_res)
  #pearson_vec <- get_corr_vec(bulk_data = sampled_bulk, decon_res = decon_sampled_res)$pearson_vec  
  #spearman_vec <- get_corr_vec(bulk_data = sampled_bulk, decon_res = decon_sampled_res)$spearman_vec

  #statistics <- list("pearson_vec" = pearson_vec, "spearman_vec" = spearman_vec)
  return(statistics_vec)
}


get_number_marker_genes <- function(decon_res){
  quantiles <- sapply(1:ncol(decon_res$basis.mvw), function(x) quantile(decon_res$basis.mvw[,x], seq(0,1,0.01))[100])
  marker_genes_ct <- lapply(1:length(quantiles), function(x) {
    which(decon_res$basis.mvw[,x]>quantiles[x])
  })
  n_marker_genes <- length(Reduce(union, marker_genes_ct))
  #n_marker_genes <- length(marker_genes)
  
  #list_markergenes <- list("marker_genes" = marker_genes, "n_marker_genes" = n_marker_genes)
  return(n_marker_genes)
}


Calculate_pvalue <- function(nrep = 500, ncores = 5, silent = TRUE, bulk_data, ...) { 
  
  ## Deconvolution of whole bulk RNA-seq dataset
  message("Executing deconvolution of the whole bulk RNA-seq dataset ..")
  decon_res <- Deconvolve_SCDC(bulk_data = bulk_data, ...)
  pearson_vec_whole <- get_test_statistics_vec(bulk_data = bulk_data, decon_res = decon_res)$pearson_vec
  spearman_vec_whole <- get_test_statistics_vec(bulk_data = bulk_data, decon_res = decon_res)$spearman_vec
  mad_vec_whole <- get_test_statistics_vec(bulk_data = bulk_data, decon_res = decon_res)$mad_vec
  rmsd_vec_whole <- get_test_statistics_vec(bulk_data = bulk_data, decon_res = decon_res)$rmsd_vec
  
  ## Deconvolution of sampled bulk RNA-seq dataset
  message("Calculation of p-value. This step takes some time ..")
  ## Calculate number of marker genes 
  #marker_genes <- get_number_marker_genes(decon_res = decon_res)
  ns_genes <- get_number_marker_genes(decon_res = decon_res)
  #bulk_no_markers <- bulk_data[- marker_genes$marker_genes,]
  #ns_genes <- marker_genes$n_marker_genes
  
  ## Approximate distribution of test statistic
  statistics_sampled <- mclapply(1:nrep, function(x) get_statistics_decon_sampled_bulk(ns_genes = ns_genes, bulk_data = bulk_data, ...),
                                 mc.cores = ncores, mc.silent = silent)
  pearson_matrix_sampled <- sapply(statistics_sampled, function(x) x$pearson_vec)
  spearman_matrix_sampled <- sapply(statistics_sampled, function(x) x$spearman_vec)
  mad_matrix_sampled <- sapply(statistics_sampled, function(x) x$mad_vec)
  rmsd_matrix_sampled <- sapply(statistics_sampled, function(x) x$rmsd_vec)

  ## Calculate p-val
  p_value_wy_pearson <- (sum(colMedians(abs(pearson_matrix_sampled)) >= median(abs(pearson_vec_whole)))+1)/(ncol(pearson_matrix_sampled)+1)
  p_value_wy_spearman <- (sum(colMedians(abs(spearman_matrix_sampled)) >= median(abs(spearman_vec_whole)))+1)/(ncol(spearman_matrix_sampled)+1)
  p_value_wy_mad <- (sum(colMedians(abs(mad_matrix_sampled)) <= median(abs(mad_vec_whole)))+1)/(ncol(mad_matrix_sampled)+1)
  p_value_wy_rmsd <- (sum(colMedians(abs(rmsd_matrix_sampled)) <= median(abs(rmsd_vec_whole)))+1)/(ncol(rmsd_matrix_sampled)+1)
  
  # pearson_matrix_sampled_sort <- apply(pearson_matrix_sampled, MARGIN= 2, FUN = sort)
  # spearman_matrix_sampled_sort <- apply(spearman_matrix_sampled, MARGIN= 2, FUN = sort)
  # p_value_cibersort_pearson <- apply(pearson_matrix_sampled_sort, MARGIN = 2, function(x) {
  #   1 - (which.min(abs(x - pearson_vec_whole)) / nrep)})
  # p_value_cibersort_spearman <- apply(spearman_matrix_sampled_sort, MARGIN = 2, function(x) {
  #   1 - (which.min(abs(x - spearman_vec_whole)) / nrep)})
  
  message("Done.")

  decon_res_pval <- list("decon_res" = decon_res, 
                         "statistics_sampled" = statistics_sampled,
                         "p_value_wy_pearson" = p_value_wy_pearson,
                         "p_value_wy_spearman" = p_value_wy_spearman,
                         "p_value_wy_mad" = p_value_wy_mad,
                         "p_value_wy_rmsd" = p_value_wy_rmsd)
                         #"p_value_cibersort_pearson" = p_value_cibersort_pearson,
                         #"p_value_cibersort_spearman" = p_value_cibersort_spearman)
  return(decon_res_pval)
} 
