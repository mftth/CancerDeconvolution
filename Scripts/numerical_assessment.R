## Mastherthesis, Melanie Fattohi
## numerical assessment of SCDC
## pseudo data
## add noise
## see in scdc how they did it
## pwert uber noise auflegen

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(RMThreshold)
set.seed(4)
## real bulk RNA-seq data and meta data
repset <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")

## qc'ed scRNA-seq data
qc_baron <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_baron_exo.RDS")
qc_segerstolpe <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/Segerstolpe_qc_exo.RDS")
qc_lawlor <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/Lawlor_qc_exo.RDS")

cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")

## generate pseudo bulk RNA-seq data from scRNA-seq data, of which we know underlying cell type proportions
## use generateBulk_allcells from SCDC
pseudo_bulk <- generateBulk_allcells(qc_segerstolpe$sc.eset.qc, ct.varname = "cluster", sample = "sample", ct.sub = cts)
pseudo_bulk_truep <- pseudo_bulk$truep
pseudo_bulk_data <- pseudo_bulk$pseudo_eset@assayData$exprs

pseudo_bulk_data_noise <- lapply(1:10, function(x) sapply(1:ncol(pseudo_bulk_data), 
                                                          function(y) jitter(pseudo_bulk_data[,y], factor = y)))
# add_noise <- function(matr, times){
#   if(times == 1){
#     return(matr)
#   } else {
#     add_noise(add.Gaussian.noise(matr, symm = FALSE), times-1)
#     return(add_noise(add.Gaussian.noise(matr, symm = FALSE), times-1))
#   }
# }


#pseudo_bulk_data_noise <- replicate(5, add.Gaussian.noise(pseudo_bulk_data, symm = FALSE), simplify = FALSE)
#pseudo_bulk_data_noise_red <- Reduce("+", pseudo_bulk_data_noise, accumulate = TRUE)
