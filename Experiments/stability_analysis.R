## Mastherthesis, Melanie Fattohi
## stability analysis of SCDC

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/visualization.R")

set.seed(25)

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


## qc'ed scRNA-seq data
qc_baron <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_baron_exo.RDS")
qc_segerstolpe <- readRDS("~/Masterthesis/CancerDeconvolution/Data/SingleCell/Segerstolpe_qc_exo.RDS")

cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")

## generate pseudo bulk RNA-seq data from scRNA-seq data, of which we know underlying cell type proportions
## use generateBulk_allcells from SCDC
pseudo_bulk <- generateBulk_allcells(qc_segerstolpe$sc.eset.qc, 
                                     ct.varname = "cluster", sample = "sample", ct.sub = cts)
pseudo_bulk_truep <- pseudo_bulk$truep
pseudo_bulk_data <- pseudo_bulk$pseudo_eset@assayData$exprs


#pseudo_bulk_data_noise <- add_noise(matr = pseudo_bulk_data, times = 10)
pseudo_bulk_data_noise <- add_noise_iteratively(matr = pseudo_bulk_data, times = 100)
names(pseudo_bulk_data_noise) <- c(paste("00", 1:9, sep = ""), paste("0", 10:99, sep = ""), "100")
pseudo_bulk_data_noise_red <- pseudo_bulk_data_noise[c(1,seq(5, 100, 5))]

#######################################
## for each matrix (i.e. original and noised ones) perform decon with framework
pseudo_decon <- Calculate_pvalue(nrep = 500, ncores = 15, silent = FALSE, bulk_data = pseudo_bulk_data,
                                 bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                 sc_data = qc_baron$sc.eset.qc, cell_types = cts,
                                 ensemble = FALSE, multiple_donors = TRUE)


pseudo_decon_noise <- lapply(pseudo_bulk_data_noise_red, 
                             function(x) Calculate_pvalue(nrep = 500, ncores = 15, silent = FALSE, bulk_data = x,
                                                          bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                                          sc_data = qc_baron$sc.eset.qc, cell_types = cts,
                                                          ensemble = FALSE, multiple_donors = TRUE))
pseudo_decon_all <- list("original" = pseudo_decon)
pseudo_decon_all <- c(pseudo_decon_all, pseudo_decon_noise)
#names(pseudo_decon_all) <- c("original", as.character(1:50))

## df: v1 = iteration, v2 = sample, v3 = pvalue
spearman_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$Spearman))
spearman_pval_all$sample <- colnames(pseudo_bulk_data)
spearman_pval_all <- melt(spearman_pval_all)
spearman_pval_all$value <- -log10(spearman_pval_all$value)
#spearman_pval_all$value <- rev(spearman_pval_all$value)
spearman_pval_all_plot <- ggplot(spearman_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC -log10(Spearman p-value)") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

pearson_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$Pearson))
pearson_pval_all$sample <- colnames(pseudo_bulk_data)
pearson_pval_all <- melt(pearson_pval_all)
pearson_pval_all$value <- -log10(pearson_pval_all$value)
#pearson_pval_all$value <- rev(pearson_pval_all$value)
pearson_pval_all_plot <- ggplot(pearson_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC -log10(Pearson p-value)") +  xlab("noise iteration")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

mad_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$mAD))
mad_pval_all$sample <- colnames(pseudo_bulk_data)
mad_pval_all <- melt(mad_pval_all)
mad_pval_all$value <- -log10(mad_pval_all$value)
#mad_pval_all$value <- rev(mad_pval_all$value)
mad_pval_all_plot <- ggplot(mad_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC -log10(mAD p-value)") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

rmsd_pval_all <- as.data.frame(sapply(pseudo_decon_all, function(x) x$p_value_per_sample$RMSD))
rmsd_pval_all$sample <- colnames(pseudo_bulk_data)
rmsd_pval_all <- melt(rmsd_pval_all)
rmsd_pval_all$value <- -log10(rmsd_pval_all$value)
#rmsd_pval_all$value <- rev(rmsd_pval_all$value)
rmsd_pval_all_plot <- ggplot(rmsd_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC -log10(RMSD p-value)") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")


pearson_scdc <- sapply(pseudo_decon_all, function(x) x$statistics_observed$pearson_vec)
rmsd_scdc <- sapply(pseudo_decon_all, function(x) x$statistics_observed$rmsd_vec)
spearman_scdc <- sapply(pseudo_decon_all, function(x) x$statistics_observed$spearman_vec)
mad_scdc <- sapply(pseudo_decon_all, function(x) x$statistics_observed$mad_vec)

pearson_scdc <- as.data.frame(pearson_scdc)
pearson_scdc$sample <- colnames(pseudo_bulk_data)
pearson_scdc <- melt(pearson_scdc)
#pearson_scdc$value <- rev(pearson_scdc$value)
pearson_scdc_plot <- ggplot(pearson_scdc, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC Pearson correlation") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rmsd_scdc <- as.data.frame(rmsd_scdc)
rmsd_scdc$sample <- colnames(pseudo_bulk_data)
rmsd_scdc <- melt(rmsd_scdc)
#rmsd_scdc$value <- rev(rmsd_scdc$value)
rmsd_scdc_plot <- ggplot(rmsd_scdc, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC RMSD") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

spearman_scdc <- as.data.frame(spearman_scdc)
spearman_scdc$sample <- colnames(pseudo_bulk_data)
spearman_scdc <- melt(spearman_scdc)
spearman_scdc_plot <- ggplot(spearman_scdc, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC Spearman") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mad_scdc <- as.data.frame(mad_scdc)
mad_scdc$sample <- colnames(pseudo_bulk_data)
mad_scdc <- melt(mad_scdc)
mad_scdc_plot <- ggplot(mad_scdc, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("SCDC mAD") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

save.image("~/Masterthesis/Workspaces/stability_analysis.RData")

#######################################
## for each matrix (i.e. original and noised ones) perform decon with cibersort/bseqsc
library(bseqsc)
setwd("~/artdeco/artdeco")
devtools::load_all()
bseqsc_config('~/artdeco/CIBERSORT.R')

pseudo_decon_cibersort <- Deconvolve_transcriptome(transcriptome_data = pseudo_bulk_data,
                                                   deconvolution_algorithm = "bseqsc",
                                                   models = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
                                                   nr_permutations = 500)

pseudo_decon_noise_cibersort <- mclapply(pseudo_bulk_data_noise_red, 
                                        function(x) Deconvolve_transcriptome(transcriptome_data = x,
                                                                            deconvolution_algorithm = "bseqsc",
                                                                            models = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
                                                                            nr_permutations = 500),
                                        mc.cores = 15, mc.silent = FALSE)
pseudo_decon_all_cibersort <- list("original" = pseudo_decon_cibersort)
pseudo_decon_all_cibersort <- c(pseudo_decon_all_cibersort, 
                                pseudo_decon_noise_cibersort)

error_cibersort <- which(sapply(pseudo_decon_all_cibersort, 
                                function(x) class(x))!="try-error")

pearson_pval_cibersort <- pseudo_decon_all_cibersort
pearson_pval_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
                                          function(x) x$P_value)
rmse_cibersort <- pseudo_decon_all_cibersort
rmse_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
                                          function(x) x$RMSE)
correlation_cibersort <- pseudo_decon_all_cibersort
correlation_cibersort[error_cibersort] <- lapply(pseudo_decon_all_cibersort[error_cibersort], 
                                                 function(x) x$Correlation)

cibersort_pearson_pval_all <- as.data.frame(do.call(cbind, pearson_pval_cibersort))
cibersort_pearson_pval_all$sample <- colnames(pseudo_bulk_data)
cibersort_pearson_pval_all <- melt(cibersort_pearson_pval_all)
cibersort_pearson_pval_all$value <- -log10(cibersort_pearson_pval_all$value)
cibersort_pearson_pval_all_plot <- ggplot(cibersort_pearson_pval_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("CIBERSORT -log10(Pearson p-value)") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

cibersort_rmse_all <- as.data.frame(do.call(cbind, rmse_cibersort))
cibersort_rmse_all$sample <- colnames(pseudo_bulk_data)
cibersort_rmse_all <- melt(cibersort_rmse_all)
cibersort_rmse_all_plot <- ggplot(cibersort_rmse_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("CIBERSORT RMSE") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cibersort_correlation_all <- as.data.frame(do.call(cbind, correlation_cibersort))
cibersort_correlation_all$sample <- colnames(pseudo_bulk_data)
cibersort_correlation_all <- melt(cibersort_correlation_all)
cibersort_correlation_all_plot <- ggplot(cibersort_correlation_all, aes(x=variable, y=value)) + 
  geom_boxplot() + ylab("CIBERSORT Pearson correlation") +  xlab("noise iteration") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

save.image("~/Masterthesis/Workspaces/stability_analysis.RData")

#######################################
library(ggplot2)
library(ggpubr)

ggarrange(rmsd_pval_all_plot,
          pearson_pval_all_plot, 
          #spearman_pval_all_plot, 
          #mad_pval_all_plot, 
          cibersort_pearson_pval_all_plot, 
          ncol = 1)
ggarrange(rmsd_scdc_plot, cibersort_rmse_all_plot, ncol = 1)
ggarrange(pearson_scdc_plot, cibersort_correlation_all_plot, ncol = 1)

save.image("~/Masterthesis/Workspaces/stability_analysis.RData")
