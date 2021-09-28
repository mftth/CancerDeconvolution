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
      
      message("Deconvolution with ENSEMBLE requires at least two scRNA-seq references!")  ## error
      bulk_ensemble <- NULL
      bulk_scdc <- NULL
      
    } else {
      
      message("Performing deconvolution with ENSEMBLE ..")
      sc_list <- lapply(sc_data, function(x) x$sc.eset.qc)
      
      bulk_ensemble <- SCDC_ENSEMBLE(bulk.eset = bulk_eset, sc.eset.list = sc_list, 
                                     ct.varname = "cluster", sample = "sample", 
                                     ct.sub =  cell_types, ...)
      message("Done.")
      bulk_scdc <- NULL
      
    }
    
  } else {
    
    if(length(sc_data) > 1){
      
      message("Deconvolution with SCDC requires only one scRNA-seq reference!") ## error
      bulk_scdc <- NULL
      bulk_ensemble <- NULL
      
    } else {
      
      sc_data <- sc_data[[1]]
      
      if(multiple_donors){
        
        message("Performing deconvolution of multiple donors with SCDC ..")  
        
        bulk_scdc <- SCDC_prop(bulk.eset = bulk_eset, sc.eset = sc_data$sc.eset.qc, 
                               ct.varname = "cluster", sample = "sample", 
                               ct.sub = cell_types, ...) 
        message("Done.")
        
      } else {
        
        message("Performing deconvolution of one donor with SCDC ..")  
        
        bulk_scdc <- SCDC_prop_ONE(bulk.eset = bulk_eset, 
                                   sc.eset = sc_data$sc.eset.qc, 
                                   ct.varname = "cluster", sample = "sample", 
                                   ct.sub = cell_types, ...) 
        message("Done.")
        
      }
      
      bulk_ensemble <- NULL
      
    }
    
  }
  
  decon_res <- list("ensemble_res" = bulk_ensemble, "scdc_res" = bulk_scdc)
  return(decon_res)
  
}
