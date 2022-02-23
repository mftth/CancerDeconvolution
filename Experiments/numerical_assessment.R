## Mastherthesis, Melanie Fattohi
## numerical assessment of SCDC

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
#library(RMThreshold)

#set.seed(4)

# add_noise <- function(matr, times, res_list = NULL){
#   if(is.null(res_list)){
#     res_list <- list()
#     #vecs <- 1:times
#   }
#   if(times == 0){
#     #res_list[[times]] <- matr
#     return(res_list)
#   } else {
#     #res_list[[times]] <- add.Gaussian.noise(matr, symm = FALSE, stddev = times*10, mean = times*10+1)
#     #res_list[[times]] <- jitter(matr, factor = times*10, amount = times*10+1)
#     #res_list[[length(vecs)-times+1]] <- matr+1
#     #res_list[[times]] <- matr+1
#     res_list[[times]] <- matr + matrix(runif(nrow(matr)*ncol(matr), min = 5, max = 500), # rnorm doesnt make enough noise
#                                        nrow=nrow(matr))
#     return(add_noise(res_list[[times]], times-1, res_list))
#   }
# }

## add noise differently via for loop
add_noise_iteratively <- function(matr, times){
  uniform_noise <- lapply(1:times, function(x) matrix(runif(nrow(matr)*ncol(matr), 
                                                            min = 5, max = 500),nrow=nrow(matr)))
  #uniform_noise <- lapply(1:times, function(x) matrix(rnorm(nrow(matr)*ncol(matr), 
  #                                                         mean = 4, sd = 8),nrow=nrow(matr)))
  
  added_noise <- list()
  for (i in 1:times) {
    if(i == 1){
      added_noise[[i]] <- matr + uniform_noise[[i]]
    } else {
      added_noise[[i]] <- added_noise[[i-1]] + uniform_noise[[i]] 
    }
  }
  return(added_noise)
}

## real bulk RNA-seq data and meta data
repset <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")

## qc'ed scRNA-seq data
qc_baron <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_baron_exo.RDS")
qc_segerstolpe <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/Segerstolpe_qc_exo.RDS")
#qc_lawlor <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/Lawlor_qc_exo.RDS")

cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")

## generate pseudo bulk RNA-seq data from scRNA-seq data, of which we know underlying cell type proportions
## use generateBulk_allcells from SCDC
pseudo_bulk <- generateBulk_allcells(qc_segerstolpe$sc.eset.qc, ct.varname = "cluster", sample = "sample", ct.sub = cts)
pseudo_bulk_truep <- pseudo_bulk$truep
pseudo_bulk_data <- pseudo_bulk$pseudo_eset@assayData$exprs


#pseudo_bulk_data_noise <- add_noise(matr = pseudo_bulk_data, times = 10)
pseudo_bulk_data_noise <- add_noise_iteratively(matr = pseudo_bulk_data, times = 10)


## for each matrix (i.e. original and noised ones) perform decon with framework
pseudo_decon <- Calculate_pvalue(nrep = 1000, ncores = 15, silent = FALSE, bulk_data = pseudo_bulk_data,
                                 bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                 sc_data = qc_baron$sc.eset.qc, cell_types = cts,
                                 ensemble = FALSE, multiple_donors = TRUE)


pseudo_decon_noise <- lapply(pseudo_bulk_data_noise, 
                             function(x) Calculate_pvalue(nrep = 1000, ncores = 15, silent = FALSE, bulk_data = x,
                                                          bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                                          sc_data = qc_baron$sc.eset.qc, cell_types = cts,
                                                          ensemble = FALSE, multiple_donors = TRUE))
pseudo_decon_all <- list("original" = pseudo_decon)
pseudo_decon_all <- c(pseudo_decon_all, pseudo_decon_noise)

#ct_props <- lapply(pseudo_decon_all, function(x) x$decon_res$prop.est.mvw)
#pearson_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_pearson)
#spearman_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_spearman)
#mad_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_mad)
#rmsd_pval <- lapply(pseudo_decon_all, function(x) x$p_value_wy_rmsd)
#scdc_metrics <- lapply(ct_props, 
#                       function(x) SCDC_peval(ptrue = pseudo_bulk_truep, pest = x, 
#                                              pest.names = "pseudo_bulk")$evals.table)
pearson_scdc <- sapply(pseudo_decon_all, function(x) x$statistics_observed$pearson_vec)
rmsd_scdc <- sapply(pseudo_decon_all, function(x) x$statistics_observed$rmsd_vec)


# plot(do.call(c, pearson_pval))
# plot(do.call(c, spearman_pval))
# plot(do.call(c, mad_pval))
# plot(do.call(c, rmsd_pval))
# plot(sapply(scdc_metrics, function(x) x[1]))
# plot(sapply(scdc_metrics, function(x) x[2]))
# plot(sapply(scdc_metrics, function(x) x[3]))

## df: v1 = iteration, v2 = sample, v3 = pvalue
spearman_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$Spearman))
spearman_pval_all$sample <- colnames(pseudo_bulk_data)
spearman_pval_all <- melt(spearman_pval_all)
spearman_pval_all$value <- rev(spearman_pval_all$value)
spearman_pval_all_plot <- ggplot(spearman_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC Spearman p-value") +  xlab("noise iteration")

pearson_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$Pearson))
pearson_pval_all$sample <- colnames(pseudo_bulk_data)
pearson_pval_all <- melt(pearson_pval_all)
#pearson_pval_all$value <- rev(pearson_pval_all$value)
pearson_pval_all_plot <- ggplot(pearson_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC Pearson p-value") +  xlab("noise iteration")

mad_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$mAD))
mad_pval_all$sample <- colnames(pseudo_bulk_data)
mad_pval_all <- melt(mad_pval_all)
mad_pval_all$value <- rev(mad_pval_all$value)
mad_pval_all_plot <- ggplot(mad_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC mAD p-value") +  xlab("noise iteration")

rmsd_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$RMSD))
rmsd_pval_all$sample <- colnames(pseudo_bulk_data)
rmsd_pval_all <- melt(rmsd_pval_all)
#rmsd_pval_all$value <- rev(rmsd_pval_all$value)
rmsd_pval_all_plot <- ggplot(rmsd_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC RMSD p-value") +  xlab("noise iteration")

pearson_scdc <- as.data.frame(pearson_scdc)
pearson_scdc$sample <- colnames(pseudo_bulk_data)
pearson_scdc <- melt(pearson_scdc)
pearson_scdc$value <- rev(pearson_scdc$value)
pearson_scdc_plot <- ggplot(pearson_scdc, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC Pearson correlation") +  xlab("noise iteration")

rmsd_scdc <- as.data.frame(rmsd_scdc)
rmsd_scdc$sample <- colnames(pseudo_bulk_data)
rmsd_scdc <- melt(rmsd_scdc)
rmsd_scdc$value <- rev(rmsd_scdc$value)
rmsd_scdc_plot <- ggplot(rmsd_scdc, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC RMSD") +  xlab("noise iteration")


#save.image("~/Masterthesis/Workspaces/numerical_assessment.RData")
#save.image("~/Masterthesis/test2.RData")

## for each matrix (i.e. original and noised ones) perform decon with cibersort/bseqsc
library(bseqsc)
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
pseudo_decon_all_cibersort <- list("original" = pseudo_decon_cibersort)
pseudo_decon_all_cibersort <- c(pseudo_decon_all_cibersort, pseudo_decon_noise_cibersort)

error_cibersort <- which(sapply(pseudo_decon_all_cibersort, 
                                function(x) class(x))!="try-error")
#ct_props_cibersort <- pseudo_decon_all_cibersort
#ct_props_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
#                                              function(x) x[,1:(ncol(x)-3)])
pearson_pval_cibersort <- pseudo_decon_all_cibersort
pearson_pval_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
                                          function(x) x$P_value)
rmse_cibersort <- pseudo_decon_all_cibersort
rmse_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
                                          function(x) x$RMSE)
correlation_cibersort <- pseudo_decon_all_cibersort
correlation_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
                                                 function(x) x$Correlation)
#scdc_metrics_cibersort <- pseudo_decon_all_cibersort
#ct_props_cibersort[error_cibersort] <- lapply(ct_props_cibersort[error_cibersort], 
#                                              function(x) as.matrix(x))
#scdc_metrics_cibersort[error_cibersort] <- lapply(ct_props_cibersort[error_cibersort], 
#                                                  function(x) SCDC_peval(ptrue = pseudo_bulk_truep, pest = x, 
#                                              pest.names = "pseudo_bulk")$evals.table)


# plot(sapply(pearson_pval_cibersort[error_cibersort], function(x) mean(x)))
# plot(sapply(rmse_cibersort[error_cibersort], function(x) mean(x)))
# plot(sapply(correlation_cibersort[error_cibersort], function(x) mean(x)))
# plot(sapply(scdc_metrics_cibersort[error_cibersort], function(x) x[1]))
# plot(sapply(scdc_metrics_cibersort[error_cibersort], function(x) x[2]))
# plot(sapply(scdc_metrics_cibersort[error_cibersort], function(x) x[3]))

cibersort_pearson_pval_all <- as.data.frame(do.call(cbind, pearson_pval_cibersort))
cibersort_pearson_pval_all$sample <- colnames(pseudo_bulk_data)
cibersort_pearson_pval_all <- melt(cibersort_pearson_pval_all)
cibersort_pearson_pval_all_plot <- ggplot(cibersort_pearson_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("CIBERSORT Pearson p-value") +  xlab("noise iteration")

cibersort_rmse_all <- as.data.frame(do.call(cbind, rmse_cibersort))
cibersort_rmse_all$sample <- colnames(pseudo_bulk_data)
cibersort_rmse_all <- melt(cibersort_rmse_all)
cibersort_rmse_all_plot <- ggplot(cibersort_rmse_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("CIBERSORT RMSE") +  xlab("noise iteration")

cibersort_correlation_all <- as.data.frame(do.call(cbind, correlation_cibersort))
cibersort_correlation_all$sample <- colnames(pseudo_bulk_data)
cibersort_correlation_all <- melt(cibersort_correlation_all)
cibersort_correlation_all_plot <- ggplot(cibersort_correlation_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("CIBERSORT Pearson correlation") +  xlab("noise iteration")


#save.image("~/Masterthesis/test2.RData")

#######################################
library(ggplot2)
library(ggpubr)

# compare pearson pvalue of canceRdeconvolution and p-value of artdeco
# pearson_pval_comparison <- data.frame(method = rep(c("SCDC", "bseqsc"), 
#                                                    each = length(error_cibersort)), 
#                                       noise_iteration = rep(error_cibersort, 2), 
#                                       pearson_pvalue = NA)
# pearson_pval_comparison$pearson_pvalue[1:length(error_cibersort)] <- 
#   sapply(pseudo_decon_all[error_cibersort], function(x) mean(x$p_value_per_sample$Pearson))
# pearson_pval_comparison$pearson_pvalue[(length(error_cibersort)+1):nrow(pearson_pval_comparison)] <- 
#   sapply(pearson_pval_cibersort[error_cibersort], function(x) mean(x))
# 
# pearson_pval_comparison_plot <- ggplot(pearson_pval_comparison,
#                                        aes(x = noise_iteration, y = pearson_pvalue, color = method)) + 
#                                 geom_point()
#pearson_pval_comparison_plot
# 
# # compare rmsd of canceRdeconvolution and rmse of artdeco
# rmsd_comparison <- data.frame(method = rep(c("SCDC", "bseqsc"), 
#                                            each = length(error_cibersort)), 
#                               noise_iteration = rep(error_cibersort, 2), 
#                               rmsd = NA)
# rmsd_comparison$rmsd[1:length(error_cibersort)] <- 
#   sapply(pseudo_decon_all[error_cibersort], function(x) mean(x$statistics_observed$rmsd_vec))
# rmsd_comparison$rmsd[(length(error_cibersort)+1):nrow(rmsd_comparison)] <- 
#   sapply(rmse_cibersort[error_cibersort], function(x) mean(x))
# 
# rmsd_comparison_plot <- ggplot(rmsd_comparison,
#                                aes(x = noise_iteration, y = rmsd, color = method)) + 
#   geom_point() # differences in each method are not visible anymore, plot them individually
#rmsd_comparison_plot
# 
# # compare pearson correlation of canceRdeconvolution and correlation of artdeco
# corr_comparison <- data.frame(method = rep(c("SCDC", "bseqsc"), 
#                                            each = length(error_cibersort)), 
#                               noise_iteration = rep(error_cibersort, 2), 
#                               corr = NA)
# corr_comparison$corr[1:length(error_cibersort)] <- 
#   sapply(pseudo_decon_all[error_cibersort], function(x) mean(x$statistics_observed$pearson_vec))
# corr_comparison$corr[(length(error_cibersort)+1):nrow(corr_comparison)] <- 
#   sapply(correlation_cibersort[error_cibersort], function(x) mean(x))
# 
# corr_comparison_plot <- ggplot(corr_comparison,
#                                aes(x = noise_iteration, y = corr, color = method)) + 
#   geom_point()
#corr_comparison_plot

# compare scdc metrics of canceRdeconvolution and scdc metrics of artdeco
# metrics_comparison <- data.frame(method = rep(c("SCDC", "bseqsc"), 
#                                            each = length(error_cibersort)), 
#                                  noise_iteration = rep(error_cibersort, 2), 
#                                  rmsd = NA, mad = NA, R = NA)
# metrics_comparison$rmsd[1:length(error_cibersort)] <- 
#   sapply(scdc_metrics[error_cibersort], function(x) x[1])
# metrics_comparison$rmsd[(length(error_cibersort)+1):nrow(metrics_comparison)] <- 
#   sapply(scdc_metrics_cibersort[error_cibersort], function(x) x[1])
# metrics_comparison$mad[1:length(error_cibersort)] <- 
#   sapply(scdc_metrics[error_cibersort], function(x) x[2])
# metrics_comparison$mad[(length(error_cibersort)+1):nrow(metrics_comparison)] <- 
#   sapply(scdc_metrics_cibersort[error_cibersort], function(x) x[2])
# metrics_comparison$R[1:length(error_cibersort)] <- 
#   sapply(scdc_metrics[error_cibersort], function(x) x[3])
# metrics_comparison$R[(length(error_cibersort)+1):nrow(metrics_comparison)] <- 
#   sapply(scdc_metrics_cibersort[error_cibersort], function(x) x[3])
# 
# metrics_rmsd_comparison_plot <- ggplot(metrics_comparison,
#                                        aes(x = noise_iteration, y = rmsd, color = method)) + 
#   geom_point()
# metrics_rmsd_comparison_plot
# 
# metrics_mad_comparison_plot <- ggplot(metrics_comparison,
#                                       aes(x = noise_iteration, y = mad, color = method)) + 
#   geom_point()
# metrics_mad_comparison_plot
# 
# metrics_R_comparison_plot <- ggplot(metrics_comparison,
#                                     aes(x = noise_iteration, y = R, color = method)) + 
#   geom_point()
# metrics_R_comparison_plot


ggarrange(pearson_pval_all_plot, spearman_pval_all_plot, 
          #mad_pval_all_plot, 
          #rmsd_pval_all_plot,
          cibersort_pearson_pval_all_plot, 
          ncol = 1)
ggarrange(rmsd_scdc_plot, cibersort_rmse_all_plot, ncol = 1)
ggarrange(pearson_scdc_plot, cibersort_correlation_all_plot, ncol = 1)

save.image("~/Masterthesis/Workspaces/numerical_assessment.RData")
