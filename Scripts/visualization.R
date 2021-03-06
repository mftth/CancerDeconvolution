## Mastherthesis, Melanie Fattohi
## generate different plots for visualization of results of analyses

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pROC)
library(ggpubr)
library(umap)

## cell type proportions plots (as heatmap or bar plots)
heatmap_proportions <- function(decon_output, clinical_characteristics = NA, ...){
  ## decon_output: output of Calculate_pvalue
  ## clinical_characteristics: data.frame with clinical variables with samples in rows as decon_output
  ## ...: further input of pheatmap
  
  heatmap_proportions <- pheatmap(decon_output$decon_res$prop.est.mvw,
                                  annotation_row = clinical_characteristics,
                                  show_rownames = FALSE, ...)
  
  return(heatmap_proportions)
}

barplot_proportions <- function(decon_output, clinical_characteristics_vec){
  ## decon_output: output of Calculate_pvalue
  ## clinical_characteristics_vec: numeric or character vector of one clinical variable

  decon_res <- decon_output$decon_res$prop.est.mvw
  if(any(is.na(clinical_characteristics_vec))){
    decon_res <- decon_res[- which(is.na(clinical_characteristics_vec)),]
    clinical_characteristics_vec <- clinical_characteristics_vec[-which(is.na(clinical_characteristics_vec))]
  }
  decon_res_molten <- reshape2::melt(decon_res)
  decon_res_molten <- cbind(decon_res_molten, rep(clinical_characteristics_vec, ncol(decon_res)))
  colnames(decon_res_molten) <- c("sample", "cell_type", "value", "clinical_characteristic")
  
  barplot_proportions <- ggplot(decon_res_molten, aes(fill = cell_type, y = value, 
                                                      x = clinical_characteristic)) + theme_bw() +
    geom_bar(position = "fill", stat = "identity")
  
  return(barplot_proportions)
}

boxplot_proportions <- function(decon_output, clinical_characteristics_vec, cell_types = NULL){
  ## decon_output: output of Calculate_pvalue
  ## clinical_characteristics_vec: numeric or character vector of one clinical variable
  ## cell_types: character vector of cell types of interest
  
  if(!is.null(cell_types)){
    decon_res <- decon_output$decon_res$prop.est.mvw[,cell_types]
  } else {
    decon_res <- decon_output$decon_res$prop.est.mvw
  }
  
  decon_res_molten <- reshape2::melt(decon_res)
  decon_res_molten <- cbind(decon_res_molten, rep(clinical_characteristics_vec, ncol(decon_res)))
  colnames(decon_res_molten) <- c("sample", "cell_type", "value", "clinical_characteristic")
  
  boxplot_proportions <- ggplot(decon_res_molten, aes(x = clinical_characteristic,
                                                      y = value, 
                                                      fill = cell_type)) +
    geom_boxplot() + ylab("cell type proportion") + theme_bw()
  
  return(boxplot_proportions)
}


## p-value plot
boxplot_pvalue <- function(decon_output_list, pvalue_type = "Spearman", technology = NULL){
  ## decon_output_list: named list of multiple decon outputs, preferably of the same scrna-seq dataset
  ## pvalue_type: c("pearson", "spearman", "mad", "rmsd)
  ## technology: vector of technologies used to create each bulk dataset
  
  if(is.null(names(decon_output_list))){
    stop("The list of deconvolution outputs has to be named!")
  }
  
  if(pvalue_type == "Pearson"){
    idx <- 1
  } else if (pvalue_type == "Spearman"){
    idx <- 2
  } else if (pvalue_type == "mAD"){
    idx <- 3
  } else if (pvalue_type == "RMSD"){
    idx <- 4
  }
  
decon_pval <- lapply(decon_output_list, function(x) x$p_value_per_sample[[idx]])
sample_names <- lapply(decon_output_list, function(x) rownames(x$decon_res$prop.est.mvw))
bulk_names <- lapply(1:length(decon_pval), 
                     function(x) rep(names(decon_pval)[x], 
                                     times = length(decon_pval[[x]])))
if(!is.null(technology)){
  technology <- lapply(1:length(decon_pval), 
                       function(x) rep(technology[x], 
                                       times = length(decon_pval[[x]])))
  boxplot_df <- data.frame("neg_log10_pvalue" = -log10(Reduce(c, decon_pval)), 
                           "sample" = Reduce(c, sample_names),
                           "bulk_dataset" = Reduce(c, bulk_names),
                           "technology" = Reduce(c, technology)) 
} else {
  boxplot_df <- data.frame("neg_log10_pvalue" = -log10(Reduce(c, decon_pval)), 
                           "sample" = Reduce(c, sample_names),
                           "bulk_dataset" = Reduce(c, bulk_names)) 
}
pval_all_plot <- ggplot(boxplot_df, aes(x=bulk_dataset, y=neg_log10_pvalue, fill = technology)) + 
  geom_boxplot() + ylab( paste0("-log10(", pvalue_type, " pvalue)", collapse = "")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

return(pval_all_plot)
}

## heatmap correlation plots of marker genes annotated with cell type proportions and clinical characteristics
heatmap_corr_genes <- function(decon_output = NULL, bulk_data, clinical_characteristics, cell_types = NULL,
                               marker_genes = NULL, colnames = FALSE, ...){
  ## decon_output: output of Calculate_pvalue
  ## bulk_data: data.frame used for deconvolution
  ## clinical_characteristics: data.frame with clinical variables with samples in rows as decon_output
  ## cell_types: character vector of cell types of interest
  ## marker_genes: vector of marker gene names
  ## colnames: should sample names be printed in the heatmap?
  ## ...: further input of pheatmap
  
  ## decon_output has to be given if there are no marker genes. Otherwise it is not required
  
  if(is.null(marker_genes) && is.null(decon_output)){
    stop("Heatmap can only be generated if marker genes are given or computed from a deconvolution result!")
  }
  
  if(is.null(marker_genes)){
    marker_genes <- get_marker_genes(decon_output$decon_res)$marker_genes
  }
  
  if(!is.null(decon_output)){
    if(!is.null(cell_types)){
      annotation_df <- cbind(clinical_characteristics, decon_output$decon_res$prop.est.mvw[,cell_types])
    } else {
      annotation_df <- cbind(clinical_characteristics, decon_output$decon_res$prop.est.mvw)
    }
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


## UMAP plot
## for clustering of bulk samples based on predicted cell type proportions
## do the cell type proportions separate the clinical characteristics into distinct clusters?
umap_plot <- function(decon_output, clinical_characteristic_vec, cell_types = NULL){
  ## decon_output: output of Calculate_pvalue
  ## clinical_characteristics_vec: numeric or character vector of one clinical variable
  ## cell_types: character vector of cell types of interest
  
  if(!is.null(cell_types)){
    decon_res <- decon_output$decon_res$prop.est.mvw[,cell_types]
  } else {
    decon_res <- decon_output$decon_res$prop.est.mvw
  }
  
  if(any(is.na(clinical_characteristic_vec))){
    decon_res <- decon_res[- which(is.na(clinical_characteristic_vec)),]
    clinical_characteristic_vec <- clinical_characteristic_vec[-which(is.na(clinical_characteristic_vec))]
  }
  
  decon_umap <- umap(decon_res)
  decon_umap_df <- data.frame(decon_umap$layout, 
                              "clinical_characteristic" = clinical_characteristic_vec,
                              row.names = rownames(decon_res))
  
  umap.plot <- ggplot(decon_umap_df, aes(x = X1, y = X2, color = clinical_characteristic)) +
    geom_point() + xlab("UMAP 1") + ylab("UMAP 2") + theme_bw() + 
    theme(legend.position="top") + labs(color = "Clinical characteristic") 
    #stat_ellipse()
  return(umap.plot)
}


## ROC curve with AUC of ML analysis
## only for classification
roc_curve <- function(labels, predictions, levels, ...){
  ## labels: observed labels of clinical variable
  ## predictions: predicted labels of variable
  ## levels: levels of clinical variable
  ## ...: further input of plot.roc
  
  roc_obj <- roc(labels, ordered(predictions, levels = levels))
  roc_curve_plot <- plot.roc(roc_obj, print.auc = TRUE, ...)
  
  return(roc_curve_plot)
}


## accuracy, sensitivity, specificity
## only for classification
## visualize these three metrics for each model
## models are given in a named list
barplot_ML_evaluation <- function(model_evaluation_list){
  ## model_evaluation_list: list of evaluations of models
  
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
    ylim(0,1) + theme_bw() +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  barplot_sens <- ggplot(data = eval_table, aes(x=Class, y=Sensitivity, fill=Model)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    #ggtitle("Sensitivity") + 
    ylab("Sensitivity") + 
    ylim(0,1) + theme_bw() +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  barplot_spec <- ggplot(data = eval_table, aes(x=Class, y=Specificity, fill=Model)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    #ggtitle("Specificity") + 
    ylab("Specificity") + 
    ylim(0,1) + theme_bw() +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  barplot_ML_eval <- ggarrange(barplot_acc, barplot_sens, barplot_spec, nrow = 1, common.legend = T)
  
  return(barplot_ML_eval)
}


boxplot_ML_sd <- function(ml_model_list, folds = 5, repeats = 10, levels = NULL){
  ## ml_model_list is a list of models given by train_ML_model$rf_model_* function
  ## folds: number of folds used for cv
  ## repeatds: number of repetitions performed for cv
  
  if(is.null(names(ml_model_list))){
    names(ml_model_list) <- sapply(1:length(ml_model_list), 
                                           function(x) paste0("model", x, collapse = ""))
  }
  
  cv_bestTune <- lapply(ml_model_list, function(x) as.numeric(x$bestTune))
  cv_list <- lapply(1:length(ml_model_list), 
                    function(x) ml_model_list[[x]]$pred[ml_model_list[[x]]$pred$mtry == cv_bestTune[[x]],])
  names(cv_list) <- names(ml_model_list)
  
  for (i in 1:length(cv_list)) {
    cv_list[[i]]$fold <- as.character(Reduce(rbind, strsplit(cv_list[[i]]$Resample, "\\."))[,1])
    cv_list[[i]]$rep <- as.character(Reduce(rbind, strsplit(cv_list[[i]]$Resample, "\\."))[,2])
    
    cv_list_repeat_acc <- c()
    cv_list_repeat_sens <- c()
    cv_list_repeat_spec <- c()
    for (j in 1:repeats) {
      if(j < 10){
        repeat_str <- paste0("Rep0", j, sep= "")
      } else {
        repeat_str <- paste0("Rep", j, sep= "")
      }
      repeat_table <- cv_list[[i]][cv_list[[i]]$rep == repeat_str,]
      repeat_confusion_matrix <- caret::confusionMatrix(repeat_table$pred, 
                                                        repeat_table$obs, 
                                                        mode = "everything")$byClass
      if(is.matrix(repeat_confusion_matrix)){
        repeat_accuracy <- unname(repeat_confusion_matrix[,'Balanced Accuracy'])
        cv_list_repeat_acc <- c(cv_list_repeat_acc, repeat_accuracy)
        repeat_sensitivity <- unname(repeat_confusion_matrix[,'Sensitivity'])
        cv_list_repeat_sens <- c(cv_list_repeat_sens, repeat_sensitivity)
        repeat_specificty <- unname(repeat_confusion_matrix[,'Specificity']) 
        cv_list_repeat_spec <- c(cv_list_repeat_spec, repeat_specificty)
      } else {
        repeat_accuracy <- unname(repeat_confusion_matrix['Balanced Accuracy'])
        cv_list_repeat_acc <- c(cv_list_repeat_acc, repeat_accuracy)
        repeat_sensitivity <- unname(repeat_confusion_matrix['Sensitivity'])
        cv_list_repeat_sens <- c(cv_list_repeat_sens, repeat_sensitivity)
        repeat_specificty <- unname(repeat_confusion_matrix['Specificity']) 
        cv_list_repeat_spec <- c(cv_list_repeat_spec, repeat_specificty)
      }
      
    }
    
    cv_list_repeat_perf <- data.frame("Accuracy" = cv_list_repeat_acc,
                                      "Sensitivity" = cv_list_repeat_sens,
                                      "Specificity" = cv_list_repeat_spec)
    
    cv_list[[i]] <- list("table_complete" = cv_list[[i]],
                         "performance_cv" = cv_list_repeat_perf)
  }
  
  perf_all_models <- Reduce(rbind,lapply(cv_list, function(x) x$performance_cv))
  perf_all_models$Model <- Reduce(c,lapply(1:length(cv_list), 
                                           function(x) rep(names(cv_list)[x], 
                                                           nrow(cv_list[[x]]$performance_cv))))
  perf_all_models_molten <- melt(perf_all_models)
  
  boxplot_ML_eval <- ggplot(perf_all_models_molten, aes(x=variable, y=value, fill = Model)) + 
    geom_boxplot() + theme_bw() + #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(0,1) + ylab("Performance value") + xlab("Performance characteristic")
  
  if(!is.null(levels)){
    roc_curves <- lapply(cv_list, function(x) roc_curve(labels = x$table_complete$obs, 
                                                        predictions = x$table_complete$pred, 
                                                        levels = levels))
    return(list("boxplots" = boxplot_ML_eval,
                "ROCcurves" = roc_curves))
  } else {
    return(list("boxplots" = boxplot_ML_eval))
  }
  
}
