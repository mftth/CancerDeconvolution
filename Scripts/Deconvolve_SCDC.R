## Mastherthesis, Melanie Fattohi
## deconvolve bulk RNA-seq data with SCDC based on scRNA-seq
## input: 
### a) bulk RNA-seq dataset (numeric matrix, samples in columns, genes in rows)
### b) metadata for bulk (matrix, metadata in columns, samples in rows)
### c) path to Expression Set object(s) of qc'ed scRNA-seq dataset(s) (vector)
### d) cell types (vector)
### e) Ensemble true/false
### f) multiple_donors true/false (maybe get the info with c))
## create Expression Set object of a) and b) using getESET()
## multiple options for deconvolution
### 1) ENSEMBLE of a) using c) if e) == TRUE
### 2) SCDC of a) using c) if e) == FALSE
### 2.1) SCDC_prop of a) using c) if f) == TRUE
### 2.2.) SCDC_prop_ONE of a) using c) if f) == FALSE
## output:
### SCDC object with predicted cell type proportions and reconstruction error

library(SCDC)


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
  #ensemble_res$basis.mvw <- basis_list[[which.max(weights)]]
  
  return(ensemble_res)
}


Deconvolve_SCDC <- function(bulk_data, bulk_meta, sc_path, cell_types, ensemble, multiple_donors, ...) {
  
  if(length(sc_path) == 1){
    
    message("Importing the ExpressionSet object of the scRNA-seq dataset ..")
    sc_data <- list(readRDS(file = sc_path))
    
  } else {
    
    message("Importing the ExpressionSet objects of the scRNA-seq datasets ..")
    sc_data <- lapply(sc_path, function(x) readRDS(file = x))
    
  }
  
  message("Creating ExpressionSet object of the bulk RNA-seq dataset ..")
  #if(!all(colnames(bulk_data) == rownames(bulk_meta))){
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
      #bulk_ensemble <- NULL
      #bulk_scdc <- NULL
      
    } else {
      
      message("Performing deconvolution with ENSEMBLE ..")
      sc_list <- lapply(sc_data, function(x) x$sc.eset.qc)
      decon_res <- SCDC_ENSEMBLE(bulk.eset = bulk_eset, sc.eset.list = sc_list, 
                                 ct.varname = "cluster", sample = "sample", ## die können auch zu ... werden. dann ist man flexibel, was die spaltennamen angeht -> make it soft coded
                                 ct.sub =  cell_types, ...) ## hier auch
      decon_res <- calc_ens_res(ensemble_output =  decon_res)
      message("Done.")
      
    }
    
  } else {
    
    if(length(sc_data) > 1){
      
      stop("Deconvolution with SCDC requires only one scRNA-seq reference!") ## error
      #bulk_scdc <- NULL
      #bulk_ensemble <- NULL
      
    } else {
      
      sc_data <- sc_data[[1]]
      
      if(multiple_donors){
        
        message("Performing deconvolution of multiple donors with SCDC ..")  
        
        decon_res <- SCDC_prop(bulk.eset = bulk_eset, sc.eset = sc_data$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", ## die können auch zu ... werden. dann ist man flexibel, was die spaltennamen angeht -> make it soft coded
                               ct.sub = cell_types, ...)  ## hier auch
        message("Done.")
        
      } else {
        
        message("Performing deconvolution of one donor with SCDC ..")  
        
        decon_res <- SCDC_prop_ONE(bulk.eset = bulk_eset, 
                                   sc.eset = sc_data$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", ## die können auch zu ... werden. dann ist man flexibel, was die spaltennamen angeht -> make it soft coded
                                   ct.sub = cell_types, ...)  ## hier auch
        message("Done.")
        
      }
      
    }
    
  }
  
  #decon_res <- list("ensemble_res" = bulk_ensemble, "scdc_res" = bulk_scdc)
  return(decon_res)
}
