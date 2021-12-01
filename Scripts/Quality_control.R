## Mastherthesis, Melanie Fattohi
## perform QC of scRNA-seq dataset using SCDC_qc() and save qc'ed scRNA-seq dataset as ExpressionSet object
## input: 
### a) scRNA-seq dataset (numeric matrix, cells in columns, genes in rows)
### b) metadata for scRNA (matrix, metadata in columns, cells in rows)
### c) path+name to save Expression Set object(s) of qc'ed scRNA-seq dataset(s) (vector)
### d) cell types (vector)
### f) multiple_donors true/false (maybe get the info with c))
## output:
### qc'ed scRNA-seq dataset saved in path+name

source("~/SCDC/SCDC/R/Basic_Functions.R")
source("~/SCDC/SCDC/R/Deconvolution.R")
library(pheatmap)

Quality_control <- function(sc_data, sc_meta, sc_path, cell_types, multiple_donors, ...){
  message("Creating ExpressionSet object of the scRNA-seq dataset ..")
  
  if(ncol(sc_data) >= nrow(sc_meta)){
    
    sc_data <- sc_data[, match(rownames(sc_meta), colnames(sc_data))]
    
  } else {
    
    sc_meta <- sc_meta[match(colnames(sc_data), rownames(sc_meta)), ]
    
  }
  
  ## filter genes
  message("Filtering genes ..")
  genes_var <- apply(sc_data, MARGIN = 1, FUN = var)
  sc_data <- sc_data[genes_var != 0, ]
  sc_data <- sc_data[rowSums(sc_data) >= 1, ]
  sc_eset <- getESET(exprs = sc_data, fdata = rownames(sc_data), pdata = sc_meta)
  
  if(multiple_donors){
    
    message("Performing QC using SCDC_qc() of scRNA-seq dataset with multiple donors. This step may take some time ..")  
    
    sc_qc <- SCDC_qc(sc.eset = sc_eset, ct.varname = "cluster", sample = "sample", ct.sub = cell_types, ...) 
    message("Done.")
    
  } else {
    
    message("Performing QC using SCDC_qc() of scRNA-seq dataset with one donor. This step may take some time ..")  
    
    sc_qc <- SCDC_qc_ONE(sc.eset = sc_eset, ct.varname = "cluster", sample = "sample", ct.sub = cell_types, ...) 
    message("Done.")
    
  }
  
  saveRDS(sc_qc, file = sc_path)
  
  return(sc_qc)
}
