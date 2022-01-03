## Mastherthesis, Melanie Fattohi
## generate different plots for visualization of results of analyses
## 1) survival plots
## 2) correlations plots (boxplots, heatmaps)
## 3) cell type proportions plots (as heatmap or bar plots)
## 4) heatmap correlation plots of marker genes annotated with cell type proportions and clinical characteristics
## 5) ROC curve with AUC of ML analysis
## 6) accuracy, sensitivity, specificity

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pROC)

## 3)
heatmap_proportions <- function(decon_output, bulk_annotation, ...){
  heatmap_proportions <- pheatmap(decon_output$decon_res$prop.est.mvw,
                                  annotation_row = bulk_annotation,
                                  show_rownames = FALSE, ...)
  
  return(heatmap_proportions)
}

barplot_proportions <- function(decon_output, bulk_annotation_vec){
  # bulk_annotation_vec is a character vector of length = nrow(decon_output$decon_res$prop.est.mvw)
  decon_res <- decon_output$decon_res$prop.est.mvw
  decon_res_molten <- reshape2::melt(decon_res)
  decon_res_molten <- cbind(decon_res_molten, rep(bulk_annotation_vec, ncol(decon_res)))
  colnames(decon_res_molten) <- c("sample", "cell_type", "value", "clinical_characteristic")
  
  barplot_proportions <- ggplot(decon_res_molten, aes(fill = cell_type, y = value, 
                                                      x = clinical_characteristic)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(barplot_proportions)
}


## 4)
heatmap_corr_genes <- function(decon_output = NULL, bulk_data, bulk_annotation, 
                               marker_genes = NULL, colnames = FALSE, ...){
  # decon_output has to be given if there are no marker genes. Otherwise it is not required
  # bulk_annotation has to be a data frame, having the colnames of bulk_data as rownames
  
  if(is.null(marker_genes) && is.null(decon_output)){
    stop("Heatmap can only be generated if marker genes are given or computed from a deconvolution result!")
  }
  
  if(is.null(marker_genes)){
    marker_genes <- get_marker_genes(decon_output$decon_res)$marker_genes
  }
  
  if(!is.null(decon_output)){
    annotation_df <- cbind(bulk_annotation, decon_output$decon_res$prop.est.mvw)
  } else {
    annotation_df <- bulk_annotation
  }
  
  bulk_marker <- bulk_data[which(rownames(bulk_data) %in% marker_genes),]
  bulk_corr <- cor(bulk_marker)
  
  heatmap_corr_genes <- pheatmap(mat = bulk_corr,
                                 show_rownames = FALSE, show_colnames = colnames,
                                 annotation_col = annotation_df, ...)
  
  return(heatmap_corr_genes)
}


## 5)
roc_curve <- function(labels, predictions, levels, ...){
  roc_obj <- roc(labels, ordered(predictions, levels = levels))
  roc_curve_plot <- plot.roc(roc_obj, print.auc = TRUE, ...)
  
  return(roc_curve_plot)
}
