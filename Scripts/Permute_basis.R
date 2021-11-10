## Mastherthesis, Melanie Fattohi
## permute basis for calculation of p-value
## 1) scdc_qc
## 2) scdc_basis
## 3) scdc_prop
## 4) calc coeff
## 5) basis permutieren, save it
## 6) scdc_prop 
## 7) calc coeff
## 8) repeat steps 5-7 500x
## 9) calc p-val

source("~/SCDC/SCDC/R/Basic_Functions.R")
source("~/SCDC/SCDC/R/Deconvolution.R")
source("~/SCDC/SCDC/R/ENSEMBLE.R")

library(parallel)
library(robustbase)

calc_ens_res <- function(ensemble_output){
  ensemble_res <- list()
  weights <- ensemble_output$w_table[1, 1:(ncol(ensemble_output$w_table)-4)] # take measure inverse_SSE as recommended by Dong et al --> hard coded
  prop_list <- lapply(ensemble_output$prop.list, function(x) x$prop.est.mvw)
  ensemble_res$prop.est.mvw <- wt_prop(weights, prop_list)
  
  basis_list <- lapply(ensemble_output$prop.list, function(x) x$basis.mvw)
  basistmp <- lapply(1:length(basis_list), function(x) basis_list[[x]]*weights[x])
  common_genes <- Reduce(intersect, lapply(basistmp, function(x) rownames(x)))
  basistmp <- lapply(basistmp, function(x) x[common_genes,])
  ensemble_res$basis.mvw <- Reduce("+", basistmp)
  
  return(ensemble_res)
}


Deconvolve_SCDC <- function(bulk_data, bulk_meta, sc_data, sc_basis, cell_types, ensemble, multiple_donors, ...) { 

  message("Creating ExpressionSet object of the bulk RNA-seq dataset ..")
  ## matching samples of bulk_data with samples of bulk_meta
  
  if(ncol(bulk_data) >= nrow(bulk_meta)){
    
    bulk_data <- bulk_data[, match(rownames(bulk_meta), colnames(bulk_data))]
    
  } else {
    
    bulk_meta <- bulk_meta[match(colnames(bulk_data), rownames(bulk_meta)), ]
    
  }
  
  bulk_eset <- getESET(exprs = bulk_data, fdata = rownames(bulk_data), pdata = bulk_meta)
  
  ## performing deconvolution
  if(ensemble){
    
    if(length(sc_data) < 2){
      
      stop("Deconvolution with ENSEMBLE requires at least two scRNA-seq references!")  ## error
      
    } else {
      
      message("Performing deconvolution with ENSEMBLE ..")
      #sc_list <- lapply(sc_data, function(x) x$sc.eset.qc)
      decon_res <- SCDC_ENSEMBLE(bulk.eset = bulk_eset, sc.eset.list = sc_data, sc.basis.list = sc_basis,
                                 ct.varname = "cluster", sample = "sample", ## die können auch zu ... werden. dann ist man flexibel, was die spaltennamen angeht -> make it soft coded
                                 ct.sub =  cell_types, ...) 
      decon_res <- calc_ens_res(ensemble_output =  decon_res)
      message("Done.")
      
    }
    
  } else {
    
    if(length(sc_data) > 1){
      
      stop("Deconvolution with SCDC requires only one scRNA-seq reference!") ## error
      
    } else {
      
      #sc_data <- sc_data[[1]]
      
      if(multiple_donors){
        
        message("Performing deconvolution of multiple donors with SCDC ..")  
        
        decon_res <- SCDC_prop(bulk.eset = bulk_eset, sc.eset = sc_data, sc.basis = sc_basis,
                               ct.varname = "cluster", sample = "sample",  ## die können auch zu ... werden. dann ist man flexibel, was die spaltennamen angeht -> make it soft coded
                               ct.sub =  cell_types, ...) 
        message("Done.")
        
      } else {
        
        message("Performing deconvolution of one donor with SCDC ..")  
        
        decon_res <- SCDC_prop_ONE(bulk.eset = bulk_eset, sc.eset = sc_data, sc.basis = sc_basis,
                                   ct.varname = "cluster", sample = "sample", ## die können auch zu ... werden. dann ist man flexibel, was die spaltennamen angeht -> make it soft coded
                                   ct.sub =  cell_types, ...) 
        message("Done.")
        
      }
      
    }
    
  }
  
  return(decon_res)
}


Create_basis <- function(sc_data, cell_types, ensemble, multiple_donors){
  if(ensemble){
    
    sc_basis <- lapply(1:length(sc_data), function(idx){
      if(multiple_donors[[idx]]){
        SCDC_basis(x = sc_data[[idx]], ct.sub = cell_types, ct.varname = "cluster", sample = "sample") 
      } else {
        SCDC_basis_ONE(x = sc_data[[idx]], ct.sub = cell_types,  ct.varname = "cluster", sample = "sample")
      }
    })
    
  } else {
    
    if(multiple_donors){
      sc_basis <- SCDC_basis(x = sc_data, ct.sub = cell_types, ct.varname = "cluster", sample = "sample")  
    } else {
      sc_basis <- SCDC_basis_ONE(x = sc_data, ct.sub = cell_types,  ct.varname = "cluster", sample = "sample")  
    }
    
  }
  return(sc_basis)
}


get_test_statistics_vec <- function(bulk_data, decon_res){

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

get_marker_genes <- function(decon_res){
  quantiles <- sapply(1:ncol(decon_res$basis.mvw), function(x) quantile(decon_res$basis.mvw[,x], seq(0,1,0.01))[100])
  marker_genes_ct <- lapply(1:length(quantiles), function(x) {
    which(decon_res$basis.mvw[,x]>quantiles[x])
  })
  marker_genes <- Reduce(union, lapply(marker_genes_ct, function(x) names(x)))
  marker_genes_idx <- Reduce(union, marker_genes_ct)
  n_marker_genes <- length(marker_genes)
  
  list_marker_genes <- list("marker_genes" = marker_genes, "marker_genes_idx" = marker_genes_idx, "n_marker_genes" = n_marker_genes)
  return(list_marker_genes)
}


get_permuted_basis_statistics <- function(marker_genes, bulk_data, bulk_meta, sc_data, sc_basis, cell_types, ensemble, multiple_donors, ...){
  
  if(!ensemble){
    sc_basis <- list(sc_basis)
  }
  
  ## generate permuted basis by shuffling gene labels of marker genes
  ## basis of only marker genes
  sc_basis_marker<- sc_basis
  for (i in 1:length(sc_basis)) {
    sc_basis_marker[[i]]$basis <- sc_basis[[i]]$basis[marker_genes$marker_genes,]
    sc_basis_marker[[i]]$basis.mvw <- sc_basis[[i]]$basis.mvw[marker_genes$marker_genes,]
  }
  sc_basis_marker_rn <- lapply(sc_basis_marker, function(x) rownames(x$basis.mvw))
  # lapply(1:length(sc_basis), function(idx){
  #   sc_basis_marker[[idx]]$basis <- sc_basis[[idx]]$basis[marker_genes$marker_genes,]
  #   sc_basis_marker[[idx]]$basis.mvw <- sc_basis[[idx]]$basis.mvw[marker_genes$marker_genes,]
  # })
  
  ## basis of only non-marker genes
  sc_basis_nmarker <- sc_basis
  nmarker_genes <- lapply(1:length(sc_basis), function(x) setdiff(rownames(sc_basis[[x]]$basis), 
                                                                  rownames(sc_basis_marker[[x]]$basis)))
  for (i in 1:length(sc_basis)) {
    sc_basis_nmarker[[i]]$basis <- sc_basis[[i]]$basis[nmarker_genes[[i]],]
    sc_basis_nmarker[[i]]$basis.mvw <- sc_basis[[i]]$basis.mvw[nmarker_genes[[i]],]
  }
  
  ## shuffle marker-gene-many genes from non-marker genes and set as labels for marker genes
  shuffled_genes_idx <- lapply(1:length(sc_basis_nmarker), function(x) sample(1:nrow(sc_basis_nmarker[[x]]$basis.mvw), 
                                                                              length(sc_basis_marker_rn[[x]])))
  for (i in 1:length(sc_basis)) {
    rownames(sc_basis_marker[[i]]$basis.mvw) <- rownames(sc_basis_nmarker[[i]]$basis.mvw)[shuffled_genes_idx[[i]]]
    rownames(sc_basis_marker[[i]]$basis) <- rownames(sc_basis_nmarker[[i]]$basis.mvw)[shuffled_genes_idx[[i]]]
  }
  
  ## set marker genes as labels for shuffled non-marker genes, i.e. switch labels between
  shuffled_genes <- lapply(sc_basis_marker_rn, function(x) sample(x))
  for (i in 1:length(sc_basis)) {
    rownames(sc_basis_nmarker[[i]]$basis.mvw)[shuffled_genes_idx[[i]]] <- shuffled_genes[[i]]
    rownames(sc_basis_nmarker[[i]]$basis)[shuffled_genes_idx[[i]]] <- shuffled_genes[[i]]
  }
  
  ## put both matrices together and shuffle rows
  shuffled_basis <- lapply(1:length(sc_basis), function(x) rbind(sc_basis_nmarker[[x]]$basis, sc_basis_marker[[x]]$basis))
  shuffled_basis.mvw <- lapply(1:length(sc_basis), function(x) rbind(sc_basis_nmarker[[x]]$basis.mvw, sc_basis_marker[[x]]$basis.mvw))
  shuffled_idx <- lapply(1:length(sc_basis), function(x) sample(1:nrow(shuffled_basis.mvw[[x]])))
  shuffled_basis <- lapply(1:length(sc_basis), function(x) shuffled_basis[[x]][shuffled_idx[[x]],])
  shuffled_basis.mvw <- lapply(1:length(sc_basis), function(x) shuffled_basis.mvw[[x]][shuffled_idx[[x]],])
  sc_basis_shuffled <- sc_basis
  for (i in 1:length(sc_basis)) {
    sc_basis_shuffled[[i]]$basis <- shuffled_basis[[i]]
    sc_basis_shuffled[[i]]$basis.mvw <- shuffled_basis.mvw[[i]]
  }
  
  if(!ensemble){
    sc_basis_shuffled <- sc_basis_shuffled[[1]]
  }
  
  #shuffled_genes <- sample(rownames(sc_basis_shuffled$basis.mvw))
  #rownames(sc_basis_shuffled$basis.mvw) <- shuffled_genes
  #rownames(sc_basis_shuffled$basis) <- shuffled_genes
  
  ## perform deconvolution of shuffled basis with Deconvolve_SCDC
  decon_shuffled_res <- Deconvolve_SCDC(bulk_data = bulk_data, bulk_meta = bulk_meta, sc_data = sc_data, sc_basis = sc_basis_shuffled,  
                                        cell_types = cell_types, ensemble = ensemble, multiple_donors = multiple_donors, ...)
  
  statistics_vec <- get_test_statistics_vec(bulk_data = bulk_data, decon_res = decon_shuffled_res)
  
  return(statistics_vec)
}


Calculate_pvalue <- function(nrep = 500, ncores = 5, silent = TRUE, bulk_data, bulk_meta, sc_data, cell_types, ensemble, multiple_donors, ...){
  
  ## Generate original basis matrix
  sc_basis <- Create_basis(sc_data = sc_data,  cell_types = cell_types, ensemble = ensemble, multiple_donors = multiple_donors)
  
  ## Deconvolution of whole bulk RNA-seq dataset with original basis
  message("Executing deconvolution of the whole bulk RNA-seq dataset ..")
  decon_res <- Deconvolve_SCDC(bulk_data = bulk_data, bulk_meta = bulk_meta, sc_data = sc_data, sc_basis = sc_basis,  
                               cell_types = cell_types, ensemble = ensemble, multiple_donors = multiple_donors, ...)
  statistics_obs <- get_test_statistics_vec(bulk_data = bulk_data, decon_res = decon_res)
  
  ## Deconvolution of shuffled basis
  message("Calculating p-value. This step takes some time ..")
  marker_genes <- get_marker_genes(decon_res = decon_res)
  
  ## Approximate distribution of test statistic
  statistics_sampled <- mclapply(1:nrep, function(x) get_permuted_basis_statistics(marker_genes = marker_genes,
                                                                                   bulk_data = bulk_data,
                                                                                   bulk_meta = bulk_meta,
                                                                                   sc_data = sc_data,
                                                                                   sc_basis = sc_basis, 
                                                                                   cell_types = cell_types, 
                                                                                   ensemble = ensemble,
                                                                                   multiple_donors = multiple_donors, ...),
                                 mc.cores = ncores, mc.silent = silent)
  pearson_matrix_sampled <- sapply(statistics_sampled, function(x) x$pearson_vec)
  spearman_matrix_sampled <- sapply(statistics_sampled, function(x) x$spearman_vec)
  mad_matrix_sampled <- sapply(statistics_sampled, function(x) x$mad_vec)
  rmsd_matrix_sampled <- sapply(statistics_sampled, function(x) x$rmsd_vec)
  
  ## Calculate p-value
  p_value_wy_pearson <- (sum(colMedians(abs(pearson_matrix_sampled)) >= median(abs(statistics_obs$pearson_vec)))+1)/(ncol(pearson_matrix_sampled)+1)
  p_value_wy_spearman <- (sum(colMedians(abs(spearman_matrix_sampled)) >= median(abs(statistics_obs$spearman_vec)))+1)/(ncol(spearman_matrix_sampled)+1)
  p_value_wy_mad <- (sum(colMedians(abs(mad_matrix_sampled)) <= median(abs(statistics_obs$mad_vec)))+1)/(ncol(mad_matrix_sampled)+1)
  p_value_wy_rmsd <- (sum(colMedians(abs(rmsd_matrix_sampled)) <= median(abs(statistics_obs$rmsd_vec)))+1)/(ncol(rmsd_matrix_sampled)+1)
  
  p_value_wy_pearson_per_sample <- sapply(1:nrow(pearson_matrix_sampled), 
                                          function(x) (sum(abs(pearson_matrix_sampled[x,]) >= abs(statistics_obs$pearson_vec[x]))+1)/(ncol(pearson_matrix_sampled)+1))
  p_value_wy_spearman_per_sample <- sapply(1:nrow(spearman_matrix_sampled), 
                                           function(x) (sum(abs(spearman_matrix_sampled[x,]) >= abs(statistics_obs$spearman_vec[x]))+1)/(ncol(spearman_matrix_sampled)+1))
  p_value_wy_mad_per_sample <- sapply(1:nrow(mad_matrix_sampled), 
                                      function(x) (sum(abs(mad_matrix_sampled[x,]) <= abs(statistics_obs$mad_vec[x]))+1)/(ncol(mad_matrix_sampled)+1))
  p_value_wy_rmsd_per_sample <- sapply(1:nrow(rmsd_matrix_sampled), 
                                       function(x) (sum(abs(rmsd_matrix_sampled[x,]) <= abs(statistics_obs$rmsd_vec[x]))+1)/(ncol(rmsd_matrix_sampled)+1))
  
  p_value_per_sample <- data.frame(Pearson = p_value_wy_pearson_per_sample,
                                   Spearman = p_value_wy_spearman_per_sample,
                                   mAD = p_value_wy_mad_per_sample,
                                   RMSD = p_value_wy_rmsd_per_sample,
                                   row.names = rownames(decon_res$prop.est.mvw))
  
  message("Done.")
  
  decon_res_pval <- list("decon_res" = decon_res, 
                         "statistics_observed" = statistics_obs,
                         "statistics_sampled" = statistics_sampled,
                         "p_value_per_sample" = p_value_per_sample,
                         "p_value_wy_pearson" = p_value_wy_pearson,
                         "p_value_wy_spearman" = p_value_wy_spearman,
                         "p_value_wy_mad" = p_value_wy_mad,
                         "p_value_wy_rmsd" = p_value_wy_rmsd)
  return(decon_res_pval)
}

