## Mastherthesis, Melanie Fattohi
## generate different plots for visualization of results of analyses
## 1) survival plots in survival_analysis.R
## 2) correlations plots (boxplots, heatmaps) in correlation_analysis.R
## 3) cell type proportions plots (as heatmap or bar plots)
## 4) heatmap correlation plots of marker genes annotated with cell type proportions and clinical characteristics
## 5) ROC curve with AUC of ML analysis
## 6) accuracy, sensitivity, specificity

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pROC)
library(ggpubr)

## 3) cell type proportions plots (as heatmap or bar plots)
heatmap_proportions <- function(decon_output, clinical_characteristics = NA, ...){
  # clinical char is dataframe (one or more) with rownames as meta data
  # deconoutput ist output von framework
  heatmap_proportions <- pheatmap(decon_output$decon_res$prop.est.mvw,
                                  annotation_row = clinical_characteristics,
                                  show_rownames = FALSE, ...)
  
  return(heatmap_proportions)
}

barplot_proportions <- function(decon_output, clinical_characteristics_vec){
  # clinical_characteristics_vec is a character vector of length = nrow(decon_output$decon_res$prop.est.mvw)
  decon_res <- decon_output$decon_res$prop.est.mvw
  decon_res_molten <- reshape2::melt(decon_res)
  decon_res_molten <- cbind(decon_res_molten, rep(clinical_characteristics_vec, ncol(decon_res)))
  colnames(decon_res_molten) <- c("sample", "cell_type", "value", "clinical_characteristic")
  
  barplot_proportions <- ggplot(decon_res_molten, aes(fill = cell_type, y = value, 
                                                      x = clinical_characteristic)) +
    geom_bar(position = "fill", stat = "identity")
  
  return(barplot_proportions)
}


## p-value plot
boxplot_pvalue <- function(decon_output_list, pvalue_type = "spearman"){
  ## decon_output_list = named list of multiple decon outputs, preferably of the same sc rna-seq dataset
  ## pvalue_type = c("pearson", "spearman", "mad", "rmsd)
  
  if(is.null(names(decon_output_list))){
    stop("The list of deconvolution outputs has to be named!")
  }
  
  if(pvalue_type == "pearson"){
    idx <- 1
  } else if (pvalue_type == "spearman"){
    idx <- 2
  } else if (pvalue_type == "mad"){
    idx <- 3
  } else if (pvalue_type == "rmsd"){
    idx <- 4
  }
  
decon_pval <- lapply(decon_output_list, function(x) x$p_value_per_sample[[idx]])
sample_names <- lapply(decon_output_list, function(x) rownames(x$decon_res$prop.est.mvw))
bulk_names <- lapply(1:length(decon_pval), 
                     function(x) rep(names(decon_pval)[x], 
                                     times = length(decon_pval[[x]])))
boxplot_df <- data.frame("neg_log10_pvalue" = -log10(Reduce(c, decon_pval)), 
                         "sample" = Reduce(c, sample_names),
                         "bulk_dataset" = Reduce(c, bulk_names)) 
pval_all_plot <- ggplot(boxplot_df, aes(x=bulk_dataset, y=neg_log10_pvalue)) + 
  geom_boxplot() + ylab("-log10(pvalue)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

return(pval_all_plot)
}

## 4) heatmap correlation plots of marker genes annotated with cell type proportions and clinical characteristics
heatmap_corr_genes <- function(decon_output = NULL, bulk_data, clinical_characteristics, 
                               marker_genes = NULL, colnames = FALSE, ...){
  # decon_output has to be given if there are no marker genes. Otherwise it is not required
  # clinical_characteristics has to be a data frame, having the colnames of bulk_data as rownames
  
  if(is.null(marker_genes) && is.null(decon_output)){
    stop("Heatmap can only be generated if marker genes are given or computed from a deconvolution result!")
  }
  
  if(is.null(marker_genes)){
    marker_genes <- get_marker_genes(decon_output$decon_res)$marker_genes
  }
  
  if(!is.null(decon_output)){
    annotation_df <- cbind(clinical_characteristics, decon_output$decon_res$prop.est.mvw)
  } else {
    annotation_df <- clinical_characteristics
  }
  
  bulk_marker <- bulk_data[which(rownames(bulk_data) %in% marker_genes),]
  bulk_corr <- cor(bulk_marker)
  
  heatmap_corr_genes <- pheatmap(mat = bulk_corr,
                                 show_rownames = FALSE, show_colnames = colnames,
                                 annotation_col = annotation_df, ...)
  
  return(heatmap_corr_genes)
}


## 5) ROC curve with AUC of ML analysis
## only for classification
roc_curve <- function(labels, predictions, levels, ...){
  roc_obj <- roc(labels, ordered(predictions, levels = levels))
  roc_curve_plot <- plot.roc(roc_obj, print.auc = TRUE, ...)
  
  return(roc_curve_plot)
}


## 6) accuracy, sensitivity, specificity
## only for classification
## visualize these three metrics for each model
## models are given in a named list
barplot_ML_evaluation <- function(model_evaluation_list){
  if(is.null(names(model_evaluation_list))){
    names(model_evaluation_list) <- sapply(1:length(model_evaluation_list), 
                                           function(x) paste0("model", x, collapse = ""))
  }
  
  eval_table_list <- lapply(model_evaluation_list, function(x) x$byClass)
  eval_classes <- lapply(model_evaluation_list, function(x) rownames(x$table))
  # create data frame ready for ggplot2
  for (idx in 1:length(eval_table_list)) {
    if(class(eval_table_list[[idx]])[1] == "numeric"){
      eval_table_list[[idx]] <- as.data.frame(t(eval_table_list[[idx]]))
      eval_table_list[[idx]]$Model <- names(eval_table_list)[idx]
      eval_table_list[[idx]]$Class <- paste0(c(eval_classes[[idx]]), collapse = "_")
    } else {
      eval_table_list[[idx]] <- as.data.frame(eval_table_list[[idx]])
      eval_table_list[[idx]]$Model <- names(eval_table_list)[idx]
      eval_table_list[[idx]]$Class <- eval_classes[[idx]]
      #sapply(strsplit(rownames(eval_table_list[[idx]]), split = "Class: "), 
                                    #       function(x) x[2])
    }
    rownames(eval_table_list[[idx]]) <- NULL
  }
  eval_table <- as.data.frame(Reduce(rbind, eval_table_list))
  colnames(eval_table) <- gsub(" ", "_", colnames(eval_table))
  
  barplot_acc <- ggplot(data = eval_table, aes(x=Class, y=Balanced_Accuracy, fill=Model)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    #ggtitle("Accuracy") + 
    ylab("Accuracy") + 
    ylim(0,1) +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  barplot_sens <- ggplot(data = eval_table, aes(x=Class, y=Sensitivity, fill=Model)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    #ggtitle("Sensitivity") + 
    ylab("Sensitivity") + 
    ylim(0,1) +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  barplot_spec <- ggplot(data = eval_table, aes(x=Class, y=Specificity, fill=Model)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    #ggtitle("Specificity") + 
    ylab("Specificity") + 
    ylim(0,1) +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  barplot_ML_eval <- ggarrange(barplot_acc, barplot_sens, barplot_spec, nrow = 1, common.legend = T)
  
  return(barplot_ML_eval)
}

