## Mastherthesis, Melanie Fattohi
## numerical assessment of SCDC
## pseudo data
## add noise
## see in scdc how they did it
## pwert uber noise auflegen

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(bseqsc)
#library(RMThreshold)

set.seed(4)

add_noise <- function(matr, times, res_list = NULL){
  if(is.null(res_list)){
    res_list <- list()
    #vecs <- 1:times
  }
  if(times == 0){
    #res_list[[times]] <- matr
    return(rev(res_list))
  } else {
    #res_list[[times]] <- add.Gaussian.noise(matr, symm = FALSE, stddev = times*10, mean = times*10+1)
    #res_list[[times]] <- jitter(matr, factor = times*10, amount = times*10+1)
    #res_list[[length(vecs)-times+1]] <- matr+1
    #res_list[[times]] <- matr+1
    res_list[[times]] <- matr + matrix(runif(nrow(matr)*ncol(matr), min = 5, max = 500), # rnorm doesnt make enough noise
                                       nrow=nrow(matr))
    return(add_noise(res_list[[times]], times-1, res_list))
  }
}

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

## add noise to pseudo bulk RNA-seq
#pseudo_bulk_data_noise <- lapply(seq(0, 100, 10), function(x) sapply(1:ncol(pseudo_bulk_data), 
#                                                          function(y) jitter(pseudo_bulk_data[,y], factor = y)))
#for (i in 1:length(pseudo_bulk_data_noise)) {
#  colnames(pseudo_bulk_data_noise[[i]]) <- colnames(pseudo_bulk_data)
#}
pseudo_bulk_data_noise <- add_noise(matr = pseudo_bulk_data, times = 10)

## for each matrix (i.e. original and noised ones) perform decon with framework
pseudo_decon <- Calculate_pvalue(ncores = 10, silent = FALSE, bulk_data = pseudo_bulk_data,
                                 bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                 sc_data = qc_baron$sc.eset.qc, cell_types = cts,
                                 ensemble = FALSE, multiple_donors = TRUE)


pseudo_decon_noise <- lapply(pseudo_bulk_data_noise, 
                             function(x) Calculate_pvalue(ncores = 15, silent = FALSE, bulk_data = x,
                                                          bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                                          sc_data = qc_baron$sc.eset.qc, cell_types = cts,
                                                          ensemble = FALSE, multiple_donors = TRUE))
pseudo_decon_all <- list(pseudo_decon)
pseudo_decon_all <- c(pseudo_decon_all, pseudo_decon_noise)

ct_props <- lapply(pseudo_decon_all, function(x) x$decon_res$prop.est.mvw)
pearson_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_pearson)
spearman_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_spearman)
mad_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_mad)
rmsd_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_rmsd)
scdc_metrics <- lapply(ct_props, 
                       function(x) SCDC_peval(ptrue = pseudo_bulk_truep, pest = x, 
                                              pest.names = "pseudo_bulk")$evals.table)

# plot(do.call(c, pearson_pval))
# plot(do.call(c, spearman_pval))
# plot(do.call(c, mad_pval))
# plot(do.call(c, rmsd_pval))
# plot(sapply(scdc_metrics, function(x) x[1]))
# plot(sapply(scdc_metrics, function(x) x[2]))
# plot(sapply(scdc_metrics, function(x) x[3]))

## for each matrix (i.e. original and noised ones) perform decon with cibersort/bseqsc
setwd("~/artdeco/artdeco")
devtools::load_all()
bseqsc_config('~/artdeco/CIBERSORT.R')

pseudo_decon_cibersort <- Deconvolve_transcriptome(transcriptome_data = pseudo_bulk_data,
                                                   deconvolution_algorithm = "bseqsc",
                                                   models = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron")

pseudo_decon_noise_cibersort <- mclapply(pseudo_bulk_data_noise, 
                                        function(x) Deconvolve_transcriptome(transcriptome_data = x,
                                                                            deconvolution_algorithm = "bseqsc",
                                                                            models = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron"),
                                        mc.cores = 3, mc.silent = FALSE)

save.image("~/Masterthesis/Workspaces/numerical_assessment.R")