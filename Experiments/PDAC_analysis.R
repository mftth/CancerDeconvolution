source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/visualization.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/survival_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/correlation_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_Deconvolution.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_MachineLearning.R")

## 62 samples; survival; tumor subtype (basal, classical, hybrid); RNA-seq
Guo_bulk <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_bulk.RDS")
Guo_meta <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_metadata.RDS")
Guo_baron_decon <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Baron/Guo_decon.RDS")
Guo_tosti_decon <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti/Guo_decon.RDS")

## 357 samples; survival; tumor subtype; microarray
## 0 NA, 1 Classical, 2 Basal
Moffitt_array_bulk <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_bulk.RDS")
Moffitt_array_meta <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_metadata.RDS")
Moffitt_array_baron_decon <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Baron/Moffitt_array_decon.RDS")
Moffitt_array_tosti_decon <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti/Moffitt_array_decon.RDS")
Moffitt_array_meta$tumor_subtype <- rep(NA, nrow(Moffitt_array_meta))
Moffitt_array_meta$tumor_subtype[Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2` == "1"] <- "Classical" 
Moffitt_array_meta$tumor_subtype[Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2` == "2"] <- "Basal" 
## exlude metastasis samples
Moffitt_array_metastasis <- which(Moffitt_array_meta$`tissue type:ch2` == "Metastasis")
Moffitt_array_bulk <- Moffitt_array_bulk[,-Moffitt_array_metastasis]
Moffitt_array_meta <- Moffitt_array_meta[-Moffitt_array_metastasis,]
Moffitt_array_baron_decon$p_value_per_sample <- Moffitt_array_baron_decon$p_value_per_sample[-Moffitt_array_metastasis,]
Moffitt_array_baron_decon$decon_res$prop.est.mvw <- Moffitt_array_baron_decon$decon_res$prop.est.mvw[-Moffitt_array_metastasis,]
Moffitt_array_baron_decon$statistics_observed$pearson_vec <- Moffitt_array_baron_decon$statistics_observed$pearson_vec[-Moffitt_array_metastasis]
Moffitt_array_baron_decon$statistics_observed$spearman_vec <- Moffitt_array_baron_decon$statistics_observed$spearman_vec[-Moffitt_array_metastasis]
Moffitt_array_baron_decon$statistics_observed$mad_vec <- Moffitt_array_baron_decon$statistics_observed$mad_vec[-Moffitt_array_metastasis]
Moffitt_array_baron_decon$statistics_observed$rmsd_vec <- Moffitt_array_baron_decon$statistics_observed$rmsd_vec[-Moffitt_array_metastasis]
Moffitt_array_tosti_decon$p_value_per_sample <- Moffitt_array_tosti_decon$p_value_per_sample[-Moffitt_array_metastasis,]
Moffitt_array_tosti_decon$decon_res$prop.est.mvw <- Moffitt_array_tosti_decon$decon_res$prop.est.mvw[-Moffitt_array_metastasis,]
Moffitt_array_tosti_decon$statistics_observed$pearson_vec <- Moffitt_array_tosti_decon$statistics_observed$pearson_vec[-Moffitt_array_metastasis]
Moffitt_array_tosti_decon$statistics_observed$spearman_vec <- Moffitt_array_tosti_decon$statistics_observed$spearman_vec[-Moffitt_array_metastasis]
Moffitt_array_tosti_decon$statistics_observed$mad_vec <- Moffitt_array_tosti_decon$statistics_observed$mad_vec[-Moffitt_array_metastasis]
Moffitt_array_tosti_decon$statistics_observed$rmsd_vec <- Moffitt_array_tosti_decon$statistics_observed$rmsd_vec[-Moffitt_array_metastasis]


## visualization of results
#baron_guo_heatmap_corr <- heatmap_corr_genes(decon_output = Guo_baron_decon, bulk_data = Guo_bulk,
#                                            clinical_characteristics = data.frame(Guo_meta$description,
#                                                                             row.names = rownames(Guo_meta)))
tosti_guo_heatmap_corr <- heatmap_corr_genes(decon_output = Guo_tosti_decon, bulk_data = Guo_bulk,
                                             clinical_characteristics = data.frame(Guo_meta$description,
                                                                                   row.names = rownames(Guo_meta)),
                                             clustering_method = "median")

#baron_moffitt_heatmap_corr <- heatmap_corr_genes(decon_output = Moffitt_array_baron_decon, bulk_data = Moffitt_array_bulk,
#                                             clinical_characteristics = data.frame(Moffitt_array_meta$tumor_subtype,
#                                                                                   row.names = rownames(Moffitt_array_meta)))
tosti_moffitt_heatmap_corr <- heatmap_corr_genes(decon_output = Moffitt_array_tosti_decon, bulk_data = Moffitt_array_bulk,
                                             clinical_characteristics = data.frame(Moffitt_array_meta$tumor_subtype,
                                                                                   row.names = rownames(Moffitt_array_meta)))


## cell type prop
#baron_Guo_prop_heatmap_prop <- heatmap_proportions(decon_output = Guo_baron_decon,
#                                              clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description, row.names = rownames(Guo_meta)))
tosti_Guo_prop_heatmap_prop <- heatmap_proportions(decon_output = Guo_tosti_decon,
                                                   clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description, 
                                                                                         row.names = rownames(Guo_meta)),
                                                   clustering_method = "complete")

#baron_Guo_prop_bar <- barplot_proportions(decon_output = decon_baron$Guo,
#                                          clinical_characteristics = Guo_meta$description)
#tosti_Guo_prop_bar <- barplot_proportions(decon_output = decon_tosti$Guo,
#                                          clinical_characteristics = Guo_meta$description)

#baron_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Moffitt_array,
#                                                        clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`, row.names = rownames(Moffitt_array_meta)))
tosti_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = Moffitt_array_tosti_decon,
                                                        clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$tumor_subtype, 
                                                                                              row.names = rownames(Moffitt_array_meta)),
                                                        clustering_method = "complete")

#baron_Moffitt_array_prop_bar <- barplot_proportions(decon_output = decon_baron$Moffitt_array,
#                                                    clinical_characteristics = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`)
#tosti_Moffitt_array_prop_bar <- barplot_proportions(decon_output = decon_tosti$Moffitt_array,
#                                                    clinical_characteristics = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`)


## survival analysis
Guo_OS <- Guo_meta$Days
Guo_Zensur <- Guo_meta$status.1
Guo_Zensur[Guo_Zensur == 1] <- 0
Guo_Zensur[Guo_Zensur == 2] <- 1
#guo_hybrid <- which(Guo_meta$description == "Hybrid")
#guo_baron_decon_surv <- Guo_baron_decon
#guo_baron_decon_surv$decon_res$prop.est.mvw <- guo_baron_decon_surv$decon_res$prop.est.mvw[-guo_hybrid,]
#guo_tosti_decon_surv <- Guo_tosti_decon
#guo_tosti_decon_surv$decon_res$prop.est.mvw <- guo_tosti_decon_surv$decon_res$prop.est.mvw[-guo_hybrid,]
#Guo_OS <- Guo_OS[-guo_hybrid]
#Guo_Zensur <- Guo_Zensur[-guo_hybrid]

baron_guo_survival <- survival_analysis(decon_output = guo_baron_decon_surv, OS = Guo_OS, censor = Guo_Zensur, 
                                        clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description[-guo_hybrid], 
                                                                              row.names = rownames(Guo_meta)[-guo_hybrid]))
tosti_guo_survival <- survival_analysis(decon_output = guo_tosti_decon_surv, OS = Guo_OS, censor = Guo_Zensur, 
                                        clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description[-guo_hybrid], 
                                                                              row.names = rownames(Guo_meta)[-guo_hybrid]))


## correlation analysis
baron_guo_correlation <- correlation_analysis(decon_output = Guo_baron_decon, 
                                              clinical_characteristic = Guo_meta$description)
tosti_guo_correlation <- correlation_analysis(decon_output = Guo_tosti_decon, 
                                              clinical_characteristic = Guo_meta$description)

baron_moffitt_correlation <- correlation_analysis(decon_output = Moffitt_array_baron_decon, 
                                                  clinical_characteristic = Moffitt_array_meta$tumor_subtype)
tosti_moffitt_correlation <- correlation_analysis(decon_output = Moffitt_array_tosti_decon, 
                                                  clinical_characteristic = Moffitt_array_meta$tumor_subtype)


## ML analysis
baron_guo_prepped <- prepare_decon_res(p_value = TRUE, decon_res = Guo_baron_decon, 
                                       clinical_char = Guo_meta$description)
baron_guo_trainRowNumbers <- createDataPartition(baron_guo_prepped$response, p = 0.8, list = FALSE)
baron_guo_train <- baron_guo_prepped[baron_guo_trainRowNumbers,]
baron_guo_test <- baron_guo_prepped[-baron_guo_trainRowNumbers,]
baron_guo_ml_model <- train_ML_model(trainData = baron_guo_train)
baron_guo_ml_pred <- test_ML_model(train_output = baron_guo_ml_model, 
                                   testData = baron_guo_test[,-ncol(baron_guo_test)], 
                                   truth_vec = baron_guo_test$response)
baron_guo_roc <- roc_curve(labels = baron_guo_test$response, 
                           predictions = baron_guo_ml_pred$predicted_reduced, levels = c("Basal", "Hybrid", "Classical"))
baron_guo_prepped2 <- baron_guo_prepped[-guo_hybrid,]
baron_guo_prepped2$response <- factor(baron_guo_prepped2$response, levels = c("Basal", "Classical"))
baron_guo_ml_model <- train_ML_model(trainData = baron_guo_prepped2)


tosti_guo_prepped <- prepare_decon_res(p_value = TRUE, decon_res = Guo_tosti_decon, 
                                       clinical_char = Guo_meta$description)
tosti_guo_trainRowNumbers <- createDataPartition(tosti_guo_prepped$response, p = 0.8, list = FALSE)
tosti_guo_train <- tosti_guo_prepped[tosti_guo_trainRowNumbers,]
tosti_guo_test <- tosti_guo_prepped[-tosti_guo_trainRowNumbers,]
tosti_guo_ml_model <- train_ML_model(trainData = tosti_guo_train)
tosti_guo_ml_pred <- test_ML_model(train_output = tosti_guo_ml_model, 
                                   testData = tosti_guo_test[,-ncol(tosti_guo_test)], 
                                   truth_vec = tosti_guo_test$response)
tosti_guo_roc <- roc_curve(labels = tosti_guo_test$response, 
                           predictions = tosti_guo_ml_pred$predicted_reduced, levels = c("Basal", "Hybrid", "Classical"))
tosti_guo_prepped2 <- tosti_guo_prepped[-guo_hybrid,]
tosti_guo_prepped2$response <- factor(tosti_guo_prepped2$response, levels = c("Basal", "Classical"))
tosti_guo_ml_model <- train_ML_model(trainData = tosti_guo_prepped2)

barplot_ML_evaluation(list(#"baron_guo_whole" = baron_guo_ml_pred$evaluation_whole,
                           "baron_guo_reduced" = baron_guo_ml_pred$evaluation_reduced,
                           #"tosti_guo_whole" = tosti_guo_ml_pred$evaluation_whole,
                           "tosti_guo_reduced" = tosti_guo_ml_pred$evaluation_reduced))

boxplot_ML_sd(list("baron_guo" = baron_guo_ml_model$rf_model_whole,
                   #"baron_guo_reduced" = baron_guo_ml_model$rf_model_reduced,
                   "tosti_guo" = tosti_guo_ml_model$rf_model_whole))#,
                   #"tosti_guo_reduced" = tosti_guo_ml_model$rf_model_reduced))


baron_moffitt_prepped <- prepare_decon_res(p_value = TRUE, decon_res = Moffitt_array_baron_decon, 
                                          clinical_char = Moffitt_array_meta$tumor_subtype)
baron_moffitt_prepped <- baron_moffitt_prepped[!is.na(baron_moffitt_prepped$response),]
#baron_moffitt_trainRowNumbers <- createDataPartition(baron_moffitt_prepped$response, p = 0.8, list = FALSE)
#baron_moffitt_train <- baron_moffitt_prepped[baron_moffitt_trainRowNumbers,]
#baron_moffitt_test <- baron_moffitt_prepped[-baron_moffitt_trainRowNumbers,]
#baron_moffitt_ml_model <- train_ML_model(trainData = baron_moffitt_train)
#baron_moffitt_ml_pred <- test_ML_model(train_output = baron_moffitt_ml_model, 
#                                       testData = baron_moffitt_test[,-ncol(baron_moffitt_test)], 
#                                       truth_vec = baron_moffitt_test$response)
baron_moffitt_ml_model <- train_ML_model(trainData = baron_moffitt_prepped)
#baron_moffitt_roc <- roc_curve(labels = baron_moffitt_test$response, 
#                           predictions = baron_moffitt_ml_pred$predicted_reduced, levels = c("Basal", "Hybrid", "Classical"))

tosti_moffitt_prepped <- prepare_decon_res(p_value = TRUE, decon_res = Moffitt_array_tosti_decon, 
                                           clinical_char = Moffitt_array_meta$tumor_subtype)
tosti_moffitt_prepped <- tosti_moffitt_prepped[!is.na(tosti_moffitt_prepped$response),]
#tosti_moffitt_trainRowNumbers <- createDataPartition(tosti_moffitt_prepped$response, p = 0.8, list = FALSE)
#tosti_moffitt_train <- tosti_moffitt_prepped[tosti_moffitt_trainRowNumbers,]
#tosti_moffitt_test <- tosti_moffitt_prepped[-tosti_moffitt_trainRowNumbers,]
#tosti_moffitt_ml_model <- train_ML_model(trainData = tosti_moffitt_train)
#tosti_moffitt_ml_pred <- test_ML_model(train_output = tosti_moffitt_ml_model, 
#                                   testData = tosti_moffitt_test[,-ncol(tosti_moffitt_test)], 
#                                   truth_vec = tosti_moffitt_test$response)
tosti_moffitt_ml_model <- train_ML_model(trainData = tosti_moffitt_prepped)
#tosti_moffitt_roc <- roc_curve(labels = tosti_moffitt_test$response, 
#                           predictions = tosti_moffitt_ml_pred$predicted_reduced, levels = c("Basal", "Hybrid", "Classical"))

boxplot_ML_sd(list(#"baron_moffitt_whole" = baron_moffitt_ml_model$rf_model_whole,
                   "baron_moffitt_reduced" = baron_moffitt_ml_model$rf_model_reduced,
                   #"tosti_moffitt_whole" = tosti_moffitt_ml_model$rf_model_whole,
                   "tosti_moffitt_reduced" = tosti_moffitt_ml_model$rf_model_reduced))
