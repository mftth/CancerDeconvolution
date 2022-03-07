## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
## bulk RNA-seq data by Guo and Moffitt
## analysis of association between ct props and tumor subtypes in PDAC (Basal-like, Classical)

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
PAAD_tosti_decon <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti/PAAD_decon.RDS")
Yang_tosti_decon <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti/Yang_decon.RDS")
PAAD_bulk <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_raw_bulk.RDS")
PAAD_meta <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_metadata.RDS")
Yang_bulk <- readRDS("~/Masterthesis/Data/Bulk/Yang/Yang_bulk.RDS")
Yang_meta <- readRDS("~/Masterthesis/Data/Bulk/Yang/Yang_metadata.RDS")

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

tosti_guo_signature_heatmap <- heatmap_corr_genes(decon_output = Guo_tosti_decon, 
                                                  bulk_data = Guo_bulk,
                                                  clinical_characteristics = data.frame(Guo_meta$description,
                                                                                        row.names = rownames(Guo_meta)),
                                                  marker_genes = basal_classical_signature$Classical)
Guo_tosti_decon2 <-  Guo_tosti_decon
Guo_tosti_decon2$decon_res$prop.est.mvw <-  Guo_tosti_decon2$decon_res$prop.est.mvw[,c("gamma", "mductal", "racinar", "sacinar")]
Guo_tosti_prop_discreate <- sapply(1:ncol(Guo_tosti_decon2$decon_res$prop.est.mvw), 
                                   function(x) continuous_to_discrete(Guo_tosti_decon2$decon_res$prop.est.mvw[,x], 
                                                                      col_name = colnames(Guo_tosti_decon2$decon_res$prop.est.mvw)[x]))
colnames(Guo_tosti_prop_discreate) <- colnames(Guo_tosti_decon2$decon_res$prop.est.mvw)
tosti_guo_signature_heatmap2 <- heatmap_corr_genes(decon_output = Guo_tosti_decon2, 
                                                   bulk_data = Guo_bulk,
                                                   clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                                        #Guo_tosti_prop_discreate,
                                                                                        row.names = rownames(Guo_meta)),
                                                  marker_genes = guo_signature$V1)

PAAD_signature_heatmap <- heatmap_corr_genes(bulk_data = PAAD_bulk, 
                                             clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                   row.names = rownames(PAAD_meta)),
                                             marker_genes = guo_signature$V1)
Yang_signature_heatmap <- heatmap_corr_genes(bulk_data = Yang_bulk, 
                                             clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`,
                                                                                   row.names = rownames(Yang_meta)),
                                             marker_genes = c(basal_classical_signature$Classical, basal_classical_signature$Basal_like))

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
Moffitt_survival <- data.frame("OS" = as.numeric(Moffitt_array_meta$`survival_months:ch2`), 
                               "censor" = as.numeric(Moffitt_array_meta$`death_event_1death_0censor:ch2`), 
                               row.names = rownames(Moffitt_array_meta))
Moffitt_survival_NA <- is.na(Moffitt_survival$OS)
Moffitt_baron_decon_surv <- Moffitt_array_baron_decon
Moffitt_baron_decon_surv$decon_res$prop.est.mvw <- Moffitt_baron_decon_surv$decon_res$prop.est.mvw[!Moffitt_survival_NA,]
Moffitt_tosti_decon_surv <- Moffitt_array_tosti_decon
Moffitt_tosti_decon_surv$decon_res$prop.est.mvw <- Moffitt_tosti_decon_surv$decon_res$prop.est.mvw[!Moffitt_survival_NA,]
Moffitt_survival <- Moffitt_survival[!Moffitt_survival_NA,]

baron_guo_survival <- survival_analysis(decon_output = guo_baron_decon_surv, OS = Guo_OS, censor = Guo_Zensur, 
                                        clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description[-guo_hybrid], 
                                                                              row.names = rownames(Guo_meta)[-guo_hybrid]))
tosti_guo_survival <- survival_analysis(decon_output = guo_tosti_decon_surv, OS = Guo_OS, censor = Guo_Zensur, 
                                        clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description[-guo_hybrid], 
                                                                              row.names = rownames(Guo_meta)[-guo_hybrid]))

baron_moffitt_survival <- survival_analysis(decon_output = Moffitt_baron_decon_surv, OS = Moffitt_survival$OS, 
                                            censor = Moffitt_survival$censor, 
                                            clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$tumor_subtype[!Moffitt_survival_NA], 
                                                                                  row.names = rownames(Moffitt_survival)))
tosti_moffitt_survival <- survival_analysis(decon_output =  Moffitt_tosti_decon_surv, OS = Moffitt_survival$OS, 
                                            censor = Moffitt_survival$censor, 
                                            clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$tumor_subtype[!Moffitt_survival_NA], 
                                                                                  row.names = rownames(Moffitt_survival)))


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
## import basal_classical signatures to create comparative baseline model
basal_classical_signature <- data.frame(t(read.table("~/Masterthesis/Data/Bulk/Moffitt/basal_classical_signature.txt",
                                                     header = FALSE, row.names = 1, sep = "\t")))
guo_signature <- data.frame(read.table("~/Masterthesis/Data/Bulk/Guo/signature_genes.txt",
                                       header = FALSE, sep = "\t"))


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


Guo_signature_genes <- Guo_bulk[c(basal_classical_signature$Basal_like, basal_classical_signature$Classical),]
Guo_signature_genes <- t(na.omit(Guo_signature_genes))
all(rownames(Guo_signature_genes) == rownames(Guo_meta))
Guo_signature_genes <- data.frame(Guo_signature_genes, "response" = Guo_meta$description)
Guo_signature_genes$response <- factor(Guo_signature_genes$response, levels = c("Basal", "Classical", "Hybrid"))
Guo_signature_trainRowNumbers <- createDataPartition(Guo_signature_genes$response, p = 0.8, list = FALSE)
Guo_signature_train <- Guo_signature_genes[Guo_signature_trainRowNumbers,]
Guo_signature_test <- Guo_signature_genes[-Guo_signature_trainRowNumbers,]
Guo_signature_test$response <- factor(Guo_signature_test$response, levels = c("Basal", "Classical", "Hybrid"))
Guo_baseline_model <- train_ML_model(trainData = Guo_signature_train)
Guo_baseline_pred <- test_ML_model(train_output = Guo_baseline_model, 
                                   testData = Guo_signature_test[,-ncol(Guo_signature_test)], 
                                   truth_vec = Guo_signature_test$response)
Guo_baseline_roc <- roc_curve(labels = Guo_signature_test$response,
                              predictions =  Guo_baseline_pred$predicted_reduced, levels = c("Basal", "Classical", "Hybrid"))
Guo_signature_genes2 <- Guo_signature_genes[-guo_hybrid,]
Guo_signature_genes2$response <- factor(Guo_signature_genes2$response, levels = c("Basal", "Classical"))
Guo_baseline_model2 <- train_ML_model(trainData = Guo_signature_genes2)
Guo_signature_genes3 <- data.frame("KRAS" = as.numeric(Guo_bulk["KRAS",]), "GATA6" = as.numeric(Guo_bulk["GATA6",]), 
                                   "response" = Guo_meta$description, row.names = rownames(Guo_meta))
Guo_signature_genes3 <- Guo_signature_genes3[-guo_hybrid,]
Guo_signature_genes3$response <- factor(Guo_signature_genes3$response, levels = c("Basal", "Classical"))
Guo_baseline_model3 <- train_ML_model(trainData = Guo_signature_genes3, feature_selection = FALSE)
Guo_signature_genes4 <- t(Guo_bulk[guo_signature$V1,])
all(rownames(Guo_signature_genes4) == rownames(Guo_meta))
Guo_signature_genes4 <- data.frame(Guo_signature_genes4, "response" = Guo_meta$description)
Guo_signature_genes4 <- Guo_signature_genes4[-guo_hybrid,]
Guo_signature_genes4$response <- factor(Guo_signature_genes4$response, levels = c("Basal", "Classical"))
Guo_baseline_model4 <- train_ML_model(trainData = Guo_signature_genes4)


barplot_ML_evaluation(list(#"baron_guo_whole" = baron_guo_ml_pred$evaluation_whole,
                           "baron_guo_reduced" = baron_guo_ml_pred$evaluation_reduced,
                           #"tosti_guo_whole" = tosti_guo_ml_pred$evaluation_whole,
                           "tosti_guo_reduced" = tosti_guo_ml_pred$evaluation_reduced,
                           "baseline_reduced" = Guo_baseline_pred$evaluation_reduced))

guo_baseline_comparison <- boxplot_ML_sd(list("baron_guo" = baron_guo_ml_model$rf_model_whole,
                                             #"baron_guo_reduced" = baron_guo_ml_model$rf_model_reduced,
                                              "tosti_guo" = tosti_guo_ml_model$rf_model_whole,
                                              "baseline_guo" = Guo_baseline_model3$rf_model_whole))#,
                                             #"tosti_guo_reduced" = tosti_guo_ml_model$rf_model_reduced))
# guo_baseline_comparison_data <- guo_baseline_comparison$data
# guo_baseline_comparison_data_acc <- guo_baseline_comparison_data[guo_baseline_comparison_data$variable == "Accuracy",]
# guo_baseline_comparison_data_acc <- guo_baseline_comparison_data_acc[guo_baseline_comparison_data_acc$Model != "baron_guo",]
# guo_baseline_comparison_data_acc$Model_code <- integer(nrow(guo_baseline_comparison_data_acc))
# guo_baseline_comparison_data_acc$Model_code[guo_baseline_comparison_data_acc$Model == "tosti_guo"] <- 1
# t.test(guo_baseline_comparison_data_acc$value~guo_baseline_comparison_data_acc$Model_code, 
#        var.equal = FALSE, alternative = "two.sided")
# guo_baseline_comparison_data_sens <- guo_baseline_comparison_data[guo_baseline_comparison_data$variable == "Sensitivity",]
# guo_baseline_comparison_data_sens <- guo_baseline_comparison_data_sens[guo_baseline_comparison_data_sens$Model != "baron_guo",]
# guo_baseline_comparison_data_sens$Model_code <- integer(nrow(guo_baseline_comparison_data_sens))
# guo_baseline_comparison_data_sens$Model_code[guo_baseline_comparison_data_sens$Model == "tosti_guo"] <- 1
# t.test(guo_baseline_comparison_data_sens$value~guo_baseline_comparison_data_sens$Model_code, 
#        var.equal = FALSE, alternative = "two.sided")
# guo_baseline_comparison_data_spec <- guo_baseline_comparison_data[guo_baseline_comparison_data$variable == "Specificity",]
# guo_baseline_comparison_data_spec <- guo_baseline_comparison_data_spec[guo_baseline_comparison_data_spec$Model != "baron_guo",]
# guo_baseline_comparison_data_spec$Model_code <- integer(nrow(guo_baseline_comparison_data_spec))
# guo_baseline_comparison_data_spec$Model_code[guo_baseline_comparison_data_spec$Model == "tosti_guo"] <- 1
# t.test(guo_baseline_comparison_data_spec$value~guo_baseline_comparison_data_spec$Model_code, 
#        var.equal = FALSE, alternative = "two.sided")



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


Moffitt_signature_genes <- Moffitt_array_bulk[c(basal_classical_signature$Basal_like, basal_classical_signature$Classical),]
Moffitt_signature_genes <- t(na.omit(Moffitt_signature_genes))
all(rownames(Moffitt_signature_genes) == rownames(Moffitt_array_meta))
Moffitt_signature_genes <- data.frame(Moffitt_signature_genes, "response" = Moffitt_array_meta$tumor_subtype)
Moffitt_signature_genes <- Moffitt_signature_genes[!is.na(Moffitt_signature_genes$response),]
Moffitt_signature_genes$response <- factor(Moffitt_signature_genes$response, levels = c("Basal", "Classical"))
Moffitt_baseline_model <- train_ML_model(trainData = Moffitt_signature_genes)

Moffitt_signature_genes2 <- data.frame("KRAS" = as.numeric(Moffitt_array_bulk["KRAS",]), "GATA6" = as.numeric(Moffitt_array_bulk["GATA6",]), 
                                   "response" = Moffitt_array_meta$tumor_subtype, row.names = rownames(Moffitt_array_meta))
Moffitt_signature_genes2 <- Moffitt_signature_genes2[!is.na(Moffitt_signature_genes2$response),]
Moffitt_signature_genes2$response <- factor(Moffitt_signature_genes2$response, levels = c("Basal", "Classical"))
Moffitt_baseline_model2 <- train_ML_model(trainData = Moffitt_signature_genes2, feature_selection = FALSE)
#Moffitt_signature_trainRowNumbers <- createDataPartition(Moffitt_signature_genes$response, p = 0.8, list = FALSE)
#Moffitt_signature_train <- Moffitt_signature_genes[Moffitt_signature_trainRowNumbers,]
#Moffitt_signature_test <- Moffitt_signature_genes[-Moffitt_signature_trainRowNumbers,]
#Moffitt_signature_test$response <- factor(Moffitt_signature_test$response, levels = c("Basal", "Classical", "Hybrid"))
#Moffitt_baseline_model <- train_ML_model(trainData = Moffitt_signature_train)
#Moffitt_baseline_pred <- test_ML_model(train_output = Moffitt_baseline_model, 
#                                   testData = Moffitt_signature_test[,-ncol(Moffitt_signature_test)], 
#                                   truth_vec = Moffitt_signature_test$response)
#Moffitt_baseline_roc <- roc_curve(labels = Moffitt_signature_test$response,
#                              predictions =  Moffitt_baseline_pred$predicted_reduced, levels = c("Basal", "Classical", "Hybrid"))



boxplot_ML_sd(list(#"baron_moffitt_whole" = baron_moffitt_ml_model$rf_model_whole,
                   "baron_moffitt_reduced" = baron_moffitt_ml_model$rf_model_reduced,
                   #"tosti_moffitt_whole" = tosti_moffitt_ml_model$rf_model_whole,
                   "tosti_moffitt_reduced" = tosti_moffitt_ml_model$rf_model_reduced,
                   "baseline_moffitt_whole" = Moffitt_baseline_model2$rf_model_whole))
                   #"baseline_moffitt_reduced" = Moffitt_baseline_model$rf_model_reduced))

boxplot_ML_sd(list("baron_guo" = baron_guo_ml_model$rf_model_whole,
                   #"baron_guo_reduced" = baron_guo_ml_model$rf_model_reduced,
                   "tosti_guo" = tosti_guo_ml_model$rf_model_whole,
                   "baseline_model" = Moffitt_baseline_model$rf_model_reduced))
