## Mastherthesis, Melanie Fattohi
## compute power set of a set of cell types
## perform deconvolution for each subset of cell types
## observe differences in p-value and performance

# - man kann grundmenge an ct angeben (default NULL)
# - man gibt an welche ct ausprobiert werden soll
# vec <- c("as", "ai", "ar")
# rje::powerSet(vec)
# - jede kombi erhält einen  namen
# - für jede ausgewählte kombi wird calc_pval ausgeführt (verteilt über #threads)
#         - es werden abbildungen erstellt über die performance jeder kombi erstellt (pwert der vier metriken)
#         - gibt es signifikante unterschiede zwischen den kombis?

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")

library(parallel)
library(robustbase)
library(rje)


repset <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")
qc_baron <- readRDS(file = "~/Praktikum/Data/Baron/qc_baron_exo.RDS")
ct_set <- c("alpha", "beta", "gamma", "delta")
sub_ct_set <- c("acinar", "ductal")

ablation_study <- function(ct_set = NULL, sub_ct_set, ...) { ## remaining input of Calculate_pvalue() except for cell_types
  
  ## compute power set of the cell types which should be experimented with
  power_ct_set <- rje::powerSet(sub_ct_set)
  
  if (is.null(ct_set)){
    cts <- power_ct_set
  } else {
    cts <- lapply(power_ct_set, function(x) c(ct_set, x))
  }
  
  names_cts <- lapply(cts, function(x) paste(unname(sapply(x, function(y) substr(y, 1, 3))), collapse = "_"))
  names(cts) <- unlist(names_cts)
  
  ## perform deconvolution for each combination
  decon_res_per_subset <- lapply(cts, function(x) Calculate_pvalue(cell_types = x, ...))
  p_value_per_subset <- lapply(decon_res_per_subset, function(x) x$p_value_per_sample)
  
  pval_pearson <- matrix(ncol = 2)
  colnames(pval_pearson) <- c("decon_res", "value")
  pval_spearman <- matrix(ncol = 2)
  colnames(pval_spearman) <- c("decon_res", "value")
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
  
  pval_metrics$pval_pearson$decon_res <- factor(pval_metrics$pval_pearson$decon_res, levels = names(p_value_per_subset))
  pval_metrics$pval_spearman$decon_res <- factor(pval_metrics$pval_spearman$decon_res, levels = names(p_value_per_subset))
  pval_metrics$pval_mad$decon_res <- factor(pval_metrics$pval_mad$decon_res, levels = names(p_value_per_subset))
  pval_metrics$pval_rmsd$decon_res <- factor(pval_metrics$pval_rmsd$decon_res, levels = names(p_value_per_subset))
  
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
  
  # give deconres with lowest median p value of each kind
  # give deconres with best median metric of each kind 
  
  
  ablation_result <- list("decon_res_per_subset" = decon_res_per_subset,
                          "boxplot_pvals" = 
                            list("boxplot_pval_pearson" = box_pearson,
                                 "boxplot_pval_spearman" = box_spearman,
                                 "boxplot_pval_mad" = box_mad,
                                 "boxplot_pval_rmsd" = box_rmsd),
                          "best_pval_models" =
                            list(),
                          "best_models" = 
                            list()
                          )
  
}





