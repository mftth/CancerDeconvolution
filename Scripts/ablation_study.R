## Mastherthesis, Melanie Fattohi
## compute power set of a set of cell types
## perform deconvolution for each subset of cell types
## observe differences in p-value and performance

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")

library(parallel)
library(robustbase)
library(rje)
library(ggplot2)
library(pheatmap)

ablation_study <- function(ct_combi_list = NULL, ct_set = NULL, sub_ct_set = NULL, res_path, ...) { 
  ## remaining input of Calculate_pvalue() except for cell_types
  
  if (is.null(ct_combi_list) && is.null(ct_set) && is.null(sub_ct_set)){
    stop("Either provide a set of cell types to be investigated or the combinations of cell types of interest!")
  }
  
  if (!is.null(ct_combi_list)){
    cts <- ct_combi_list
  } else {
    ## compute power set of the cell types which should be experimented with
    message("Create all possible combinations of the cell types ..")
    power_ct_set <- rje::powerSet(sub_ct_set)
    
    if (is.null(ct_set)){
      cts <- power_ct_set[-1]
    } else {
      cts <- lapply(power_ct_set, function(x) c(ct_set, x))
    }
  }
  
  names_cts <- lapply(cts, function(x) paste(unname(sapply(x, function(y) substr(y, 1, 3))), collapse = "_"))
  names(cts) <- unlist(names_cts)
  
  ## perform deconvolution for each combination
  message("Perform deconvolution for each cell type combination. This step takes time ..")
  decon_res_per_subset <- lapply(1:length(cts), function(x){
    decon_res_x <- Calculate_pvalue(cell_types = unname(cts[x])[[1]], ...)
    saveRDS(decon_res_x, file = paste(res_path, "/", names(cts)[x], ".RDS", sep = ""))
    })
  
  decon_res_per_subset <- lapply(1:length(cts), function(x) readRDS(file = paste(res_path, "/", names(cts)[x], ".RDS", sep = "")))
  names(decon_res_per_subset) <- names(cts)
  
  message("Ablation study done.")
  message("Visualize results ..")
  ## box plots
  p_value_per_subset <- lapply(decon_res_per_subset, function(x) x$p_value_per_sample)
  
  pval_pearson <- matrix(ncol = 2)
  colnames(pval_pearson) <- c("decon_res", "value")
  pval_spearman <- matrix(ncol = 2)
  colnames(pval_spearman) <- c("decon_res", "value")q
  pval_mad <- matrix(ncol = 2)
  colnames(pval_mad) <- c("decon_res", "value")
  pval_rmsd <- matrix(ncol = 2)
  colnames(pval_rmsd) <- c("decon_res", "value")
  
  for (idx in 1:length(decon_res_per_subset)) {
    tmp_matrix_pearson <- matrix(c(rep(names(decon_res_per_subset)[idx], nrow(p_value_per_subset[[idx]])), 
                                   p_value_per_subset[[idx]]$Pearson),
                                 ncol = 2, nrow = nrow(p_value_per_subset[[idx]])) 
    tmp_matrix_spearman <- matrix(c(rep(names(decon_res_per_subset)[idx], nrow(p_value_per_subset[[idx]])), 
                                   p_value_per_subset[[idx]]$Spearman),
                                 ncol = 2, nrow = nrow(p_value_per_subset[[idx]]))
    tmp_matrix_mad <- matrix(c(rep(names(decon_res_per_subset)[idx], nrow(p_value_per_subset[[idx]])), 
                                   p_value_per_subset[[idx]]$mAD),
                                 ncol = 2, nrow = nrow(p_value_per_subset[[idx]]))
    tmp_matrix_rmsd <- matrix(c(rep(names(decon_res_per_subset)[idx], nrow(p_value_per_subset[[idx]])), 
                                   p_value_per_subset[[idx]]$RMSD),
                                 ncol = 2, nrow = nrow(p_value_per_subset[[idx]]))
    
    pval_pearson <- rbind(pval_pearson, tmp_matrix_pearson)
    pval_spearman <- rbind(pval_spearman, tmp_matrix_spearman)
    pval_mad <- rbind(pval_mad, tmp_matrix_mad)
    pval_rmsd <- rbind(pval_rmsd, tmp_matrix_rmsd)
  }
  
  pval_metrics <- list("pval_pearson" = pval_pearson, "pval_spearman" = pval_spearman,
                       "pval_mad" = pval_mad, "pval_rmsd" = pval_rmsd)
  pval_metrics <- lapply(pval_metrics, function(x) as.data.frame(x[-1,]))
  
  pval_metrics$pval_pearson$value <- as.numeric(pval_metrics$pval_pearson$value)
  pval_metrics$pval_spearman$value <- as.numeric(pval_metrics$pval_spearman$value)
  pval_metrics$pval_mad$value <- as.numeric(pval_metrics$pval_mad$value)
  pval_metrics$pval_rmsd$value <- as.numeric(pval_metrics$pval_rmsd$value)
  
  pval_metrics$pval_pearson$decon_res <- factor(pval_metrics$pval_pearson$decon_res, 
                                                levels = names(p_value_per_subset))
  pval_metrics$pval_spearman$decon_res <- factor(pval_metrics$pval_spearman$decon_res, 
                                                 levels = names(p_value_per_subset))
  pval_metrics$pval_mad$decon_res <- factor(pval_metrics$pval_mad$decon_res, 
                                            levels = names(p_value_per_subset))
  pval_metrics$pval_rmsd$decon_res <- factor(pval_metrics$pval_rmsd$decon_res, 
                                             levels = names(p_value_per_subset))
  
  box_pearson <- ggplot(pval_metrics$pval_pearson, aes(x = decon_res, y = value)) +
    geom_boxplot() +
    labs(title = "P-value of Pearson correlation of deconvolution result") + 
    ylab("p-value") +
    xlab("cell type set") +
    coord_flip()
  
  box_spearman <- ggplot(pval_metrics$pval_spearman, aes(x = decon_res, y = value)) +
    geom_boxplot() +
    labs(title = "P-value of Spearman correlation of deconvolution result") + 
    ylab("p-value") +
    xlab("cell type set") +
    coord_flip()
  
  box_mad <- ggplot(pval_metrics$pval_mad, aes(x = decon_res, y = value)) +
    geom_boxplot() +
    labs(title = "P-value of mAD of deconvolution result") +
    ylab("p-value") +
    xlab("cell type set") +
    coord_flip()
  
  box_rmsd <- ggplot(pval_metrics$pval_rmsd, aes(x = decon_res, y = value)) +
    geom_boxplot() +
    labs(title = "P-value of RMSD of deconvolution result") + 
    ylab("p-value") +
    xlab("cell type set") +
    coord_flip()
  
  ## heatmaps with pheatmap
  pearson_pheat <- Reduce(rbind, lapply(p_value_per_subset, function(x) t(x$Pearson)))
  colnames(pearson_pheat) <- rownames(p_value_per_subset[[1]])
  rownames(pearson_pheat) <- names(p_value_per_subset)
  pheatmap_pearson <- pheatmap(pearson_pheat, 
                               main = "P-value of Pearson correlation of deconvolution results") 
  spearman_pheat <- Reduce(rbind, lapply(p_value_per_subset, function(x) t(x$Spearman)))
  colnames(spearman_pheat) <- rownames(p_value_per_subset[[1]])
  rownames(spearman_pheat) <- names(p_value_per_subset)
  pheatmap_spearman <- pheatmap(spearman_pheat, 
                               main = "P-value of Spearman correlation of deconvolution results") 
  mad_pheat <- Reduce(rbind, lapply(p_value_per_subset, function(x) t(x$mAD)))
  colnames(mad_pheat) <- rownames(p_value_per_subset[[1]])
  rownames(mad_pheat) <- names(p_value_per_subset)
  pheatmap_mad <- pheatmap(mad_pheat,  
                               main = "P-value of mAD correlation of deconvolution results") 
  rmsd_pheat <- Reduce(rbind, lapply(p_value_per_subset, function(x) t(x$RMSD)))
  colnames(rmsd_pheat) <- rownames(p_value_per_subset[[1]])
  rownames(rmsd_pheat) <- names(p_value_per_subset)
  pheatmap_rmsd <- pheatmap(rmsd_pheat, 
                               main = "P-value of RMSD correlation of deconvolution results") 
  
  # give deconres with lowest median p value of each kind
  pearson_min_median_pval <- sapply(decon_res_per_subset, 
                                              function(x) median(x$p_value_per_sample$Pearson))
  pearson_min_median_pval <- pearson_min_median_pval[which(pearson_min_median_pval == 
                                                             min(pearson_min_median_pval))]
  
  spearman_min_median_pval <- sapply(decon_res_per_subset, 
                                              function(x) median(x$p_value_per_sample$Spearman))
  spearman_min_median_pval <- spearman_min_median_pval[which(spearman_min_median_pval == 
                                                               min(spearman_min_median_pval))]
  
  mad_min_median_pval <- sapply(decon_res_per_subset, 
                                              function(x) median(x$p_value_per_sample$mAD))
  mad_min_median_pval <- mad_min_median_pval[which(mad_min_median_pval ==
                                                     min(mad_min_median_pval))]
  
  rmsd_min_median_pval <- sapply(decon_res_per_subset, 
                                              function(x) median(x$p_value_per_sample$RMSD))
  rmsd_min_median_pval <- rmsd_min_median_pval[which(rmsd_min_median_pval ==
                                                       min(rmsd_min_median_pval))]
    
  # give deconres with best median metric of each kind 
  pearson_max_median <- sapply(decon_res_per_subset, 
                                         function(x) median(x$statistics_observed$pearson_vec))
  pearson_max_median <- pearson_max_median[which(pearson_max_median ==
                                                   max(pearson_max_median))]
  
  spearman_max_median <- sapply(decon_res_per_subset, 
                                          function(x) median(x$statistics_observed$spearman_vec))
  spearman_max_median <- spearman_max_median[which(spearman_max_median ==
                                               max(spearman_max_median))]
  
  mad_min_median <- sapply(decon_res_per_subset, 
                                     function(x) median(x$statistics_observed$mad_vec))
  mad_min_median <- mad_min_median[which(mad_min_median ==
                                           min(mad_min_median))]
  
  rmsd_min_median <- sapply(decon_res_per_subset, 
                                      function(x) median(x$statistics_observed$rmsd_vec))
  rmsd_min_median <- rmsd_min_median[which(rmsd_min_median ==
                                             min(rmsd_min_median))]
  
  ablation_result <- list("decon_res_per_subset" = decon_res_per_subset,
                          "boxplot_pvals" = 
                            list("boxplot_pval_pearson" = box_pearson,
                                 "boxplot_pval_spearman" = box_spearman,
                                 "boxplot_pval_mad" = box_mad,
                                 "boxplot_pval_rmsd" = box_rmsd),
                          "heatmap_pvals" =
                            list("pheatmap_pval_pearson" = pheatmap_pearson,
                                 "pheatmap_pval_spearman" = pheatmap_spearman,
                                 "pheatmap_pval_mad" = pheatmap_mad,
                                 "pheatmap_pval_rmsd" = pheatmap_rmsd),
                          "best_pval_models" =
                            list("min_median_pval_pearson" = pearson_min_median_pval,
                                 "min_median_pval_spearman" = spearman_min_median_pval,
                                 "min_median_pval_mad" = mad_min_median_pval,
                                 "min_median_pval_rmsd" = rmsd_min_median_pval),
                          "best_models" = 
                            list("max_median_pearson" = pearson_max_median,
                                 "max_median_spearman" = spearman_max_median,
                                 "min_median_mad" = mad_min_median,
                                 "min_medianl_rmsd" = rmsd_min_median)
                          )
  return(ablation_result)
}
