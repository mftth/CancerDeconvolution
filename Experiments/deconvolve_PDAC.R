## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
## scRNA-seq von PDAC werden noch gesucht
## Deconvolve Yang, PAAD, Moffitt (Seq + array), Guo, Kirby, Janky
## nreps = 1000, ncores = 15
## survival analysis
## correlation analysis
## ML anaylsis

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/visualization.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/survival_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/correlation_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_Deconvolution.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_MachineLearning.R")

## bulk RNA-seq datasets
## 183 samples; survival; tumor grading; RNA-seq
PAAD_bulk <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_raw_bulk.RDS")
PAAD_meta <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_metadata.RDS")
## 69 samples; survival; tumor grading; microarray
Yang_bulk <- readRDS("~/Masterthesis/Data/Bulk/Yang/Yang_bulk.RDS")
Yang_meta <- readRDS("~/Masterthesis/Data/Bulk/Yang/Yang_metadata.RDS")
## 62 samples; survival; tumor subtype (basal, classical, hybrid); RNA-seq
Guo_bulk <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_bulk.RDS")
Guo_meta <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_metadata.RDS")
## 131 samples; microarray
# Janky_bulk <- readRDS("~/Masterthesis/Data/Bulk/Janky/Janky_bulk.RDS")
# Janky_meta <- readRDS("~/Masterthesis/Data/Bulk/Janky/Janky_metadata.RDS")
# ## 51 samples; survival; RNA-seq
# Kirby_bulk <- readRDS("~/Masterthesis/Data/Bulk/Kirby/Kirby_bulk.RDS")
# Kirby_meta <- readRDS("~/Masterthesis/Data/Bulk/Kirby/Kirby_metadata.RDS")
## 61 samples; RNA-seq
## 15 primary tumors, 37 pancreatic cancer patient-derived xenografts (PDXs), 
## 3 PDAC cell lines and 6 cancer-associated fibroblast (CAF)
# Moffitt_seq_bulk <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_seq_bulk.RDS")
# Moffitt_seq_meta <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_seq_metadata.RDS")
## 357 samples; survival; tumor subtype; microarray
Moffitt_array_bulk <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_bulk.RDS")
Moffitt_array_meta <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_metadata.RDS")

bulk_list <- list("PAAD" = PAAD_bulk,
                  "Yang" = Yang_bulk,
                  "Guo" = Guo_bulk,
                  #"Janky" = Janky_bulk,
                  #"Kirby" = Kirby_bulk,
                  #"Moffitt_seq" = Moffitt_seq_bulk,
                  "Moffitt_array" = Moffitt_array_bulk,
                  "Flowers" = Flowers_bulk)
bulk_meta_list <- list("PAAD" = PAAD_meta,
                       "Yang" = Yang_meta,
                       "Guo" = Guo_meta,
                       #"Janky" = Janky_meta,
                       #"Kirby" = Kirby_meta,
                       #"Moffitt_seq" = Moffitt_seq_meta,
                       "Moffitt_array" = Moffitt_array_meta,
                       "Flowers" = Flowers_meta)

## single-cell RNA-seq datasets
qc_baron_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_baron_exo.RDS")
qc_segerstolpe_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/Segerstolpe_qc_exo.RDS")
qc_lawlor_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/Lawlor_qc_exo.RDS")
qc_tosti_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_tosti.RDS")
sc_list <- list("Baron" = qc_baron_sc$sc.eset.qc, 
                "Segerstolpe" = qc_segerstolpe_sc$sc.eset.qc, 
                "Lawlor" = qc_lawlor_sc$sc.eset.qc)

## perform deconvolution
reps <- 1000
ncores <- 15
cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")
cts_tosti <- c("sacinar", "racinar", "iacinar", "ductal", "mductal", "alpha", "beta", "gamma", "delta")

###
res_path_baron <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Baron"
decon_baron <- lapply(1:length(bulk_list), function(x) {
  decon_baron_x <- Calculate_pvalue(nrep = reps, ncores = ncores, silent = FALSE, 
                                    bulk_data = bulk_list[[x]], bulk_meta = bulk_meta_list[[x]],
                                    sc_data = qc_baron_sc$sc.eset.qc, cell_types = cts,
                                    ensemble = FALSE, multiple_donors = TRUE)
  saveRDS(decon_baron_x, file = paste(res_path_baron, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})

decon_baron <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_baron, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_baron) <- names(bulk_list)
###
###
res_path_tosti <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti"
decon_tosti <- lapply(1:length(bulk_list), function(x) {
  decon_tosti_x <- Calculate_pvalue(nrep = reps, ncores = ncores, silent = FALSE, 
                                    bulk_data = bulk_list[[x]], bulk_meta = bulk_meta_list[[x]],
                                    sc_data = qc_tosti_sc$sc.eset.qc, cell_types = cts_tosti,
                                    ensemble = FALSE, multiple_donors = TRUE)
  saveRDS(decon_tosti_x, file = paste(res_path_tosti, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})

decon_tosti <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_tosti, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_tosti) <- names(bulk_list)
###
###
res_path_ensemble <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Ensemble"
decon_ensemble <- lapply(1:length(bulk_list), function(x) {
  decon_ensemble_x <- Calculate_pvalue(nrep = reps, ncores = ncores, silent = FALSE, 
                                    bulk_data = bulk_list[[x]], bulk_meta = bulk_meta_list[[x]],
                                    sc_data = sc_list, cell_types = cts,
                                    ensemble = TRUE, multiple_donors = c(TRUE, TRUE, FALSE))
  saveRDS(decon_ensemble_x, file = paste(res_path_ensemble, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})

decon_ensemble <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_ensemble, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_ensemble) <- names(bulk_list)
###
###
## remove cell line samples from Moffitseq
all(rownames(decon_baron$Moffitt_seq$p_value_per_sample) == rownames(Moffitt_seq_meta))
all(rownames(decon_baron$Moffitt_seq$decon_res$prop.est.mvw) == rownames(Moffitt_seq_meta))
all(rownames(decon_tosti$Moffitt_seq$p_value_per_sample) == rownames(Moffitt_seq_meta))
all(rownames(decon_tosti$Moffitt_seq$decon_res$prop.est.mvw) == rownames(Moffitt_seq_meta))
Moffitt_seq_primary <- which(Moffitt_seq_meta$Primary == 1)
decon_baron$Moffitt_seq$p_value_per_sample <- decon_baron$Moffitt_seq$p_value_per_sample[Moffitt_seq_primary,]
decon_baron$Moffitt_seq$decon_res$prop.est.mvw <- decon_baron$Moffitt_seq$decon_res$prop.est.mvw[Moffitt_seq_primary,]
decon_tosti$Moffitt_seq$p_value_per_sample <- decon_tosti$Moffitt_seq$p_value_per_sample[Moffitt_seq_primary,]
decon_tosti$Moffitt_seq$decon_res$prop.est.mvw <- decon_tosti$Moffitt_seq$decon_res$prop.est.mvw[Moffitt_seq_primary,]


technology = c("RNA-seq", "microarray", "RNA-seq", "microarray", "RNA-seq")
baron_pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = decon_baron,
                                              pvalue_type = "Spearman", technology = technology) + ggtitle("Reference: Baron et al.") 
tosti_pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = decon_tosti,
                                              pvalue_type = "Spearman", technology = technology) + ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_spearman, tosti_pval_boxplot_spearman, common.legend = TRUE) 

baron_pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = decon_baron,
                                              pvalue_type = "Pearson", technology = technology) + ggtitle("Reference: Baron et al.")
tosti_pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = decon_tosti,
                                              pvalue_type = "Pearson", technology = technology) + ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_pearson, tosti_pval_boxplot_pearson, common.legend = TRUE)

baron_pval_boxplot_mad <- boxplot_pvalue(decon_output_list = decon_baron,
                                              pvalue_type = "mAD", technology = technology) + ggtitle("Reference: Baron et al.")
tosti_pval_boxplot_mad <- boxplot_pvalue(decon_output_list = decon_tosti,
                                              pvalue_type = "mAD", technology = technology) + ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_mad, tosti_pval_boxplot_mad, common.legend = TRUE)

baron_pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_baron,
                                         pvalue_type = "RMSD", technology = technology) + 
                           ggtitle("Reference: Baron et al.")
tosti_pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_tosti,
                                         pvalue_type = "RMSD", technology = technology) + 
                           ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_rmsd, tosti_pval_boxplot_rmsd, common.legend = TRUE) 

## grading
baron_yang_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Yang,
                                               clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`, 
                                                                                     #stage = Yang_meta$`Stage:ch1`, 
                                                                                     row.names = rownames(Yang_meta)))
tosti_yang_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Yang,
                                               clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`, 
                                                                                     #stage = Yang_meta$`Stage:ch1`, 
                                                                                     row.names = rownames(Yang_meta)))

baron_yang_prop_bar <- barplot_proportions(decon_output = decon_baron$Yang,
                                               clinical_characteristics = Yang_meta$`grading:ch1`)
tosti_yang_prop_bar <- barplot_proportions(decon_output = decon_tosti$Yang,
                                               clinical_characteristics = Yang_meta$`grading:ch1`)

baron_PAAD_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$PAAD,
                                               clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                     row.names = rownames(PAAD_meta)))
tosti_PAAD_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$PAAD,
                                               clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                     row.names = rownames(PAAD_meta)), 
                                               clustering_method = "ward.D2")

baron_PAAD_prop_bar <- barplot_proportions(decon_output = decon_baron$PAAD,
                                           clinical_characteristics = PAAD_meta$neoplasm_histologic_grade)
tosti_PAAD_prop_bar <- barplot_proportions(decon_output = decon_tosti$PAAD,
                                               clinical_characteristics = PAAD_meta$neoplasm_histologic_grade)

## tumor subtype
baron_Guo_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Guo,
                                               clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description, row.names = rownames(Guo_meta)))
tosti_Guo_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Guo,
                                               clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description, row.names = rownames(Guo_meta)))

baron_Guo_prop_bar <- barplot_proportions(decon_output = decon_baron$Guo,
                                           clinical_characteristics = Guo_meta$description)
tosti_Guo_prop_bar <- barplot_proportions(decon_output = decon_tosti$Guo,
                                           clinical_characteristics = Guo_meta$description)

baron_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Moffitt_array,
                                              clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`, row.names = rownames(Moffitt_array_meta)))
tosti_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Moffitt_array,
                                              clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`, row.names = rownames(Moffitt_array_meta)))

baron_Moffitt_array_prop_bar <- barplot_proportions(decon_output = decon_baron$Moffitt_array,
                                          clinical_characteristics = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`)
tosti_Moffitt_array_prop_bar <- barplot_proportions(decon_output = decon_tosti$Moffitt_array,
                                          clinical_characteristics = Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2`)


##################
## survival analysis
Yang_OS <- as.numeric(Yang_meta$`survival months:ch1`)
Yang_censor <- as.numeric(Yang_meta$`survival status:ch1`)
Yang_OS <- Yang_OS[-c(32, 66, 67, 68)]
Yang_censor <- Yang_censor[-c(32, 66, 67, 68)]
PAAD_OS <- rep(NA, nrow(PAAD_meta))
PAAD_OS[which(is.na(PAAD_meta$days_to_death))] <- PAAD_meta$days_to_last_followup[which(is.na(PAAD_meta$days_to_death))]
PAAD_OS[which(is.na(PAAD_meta$days_to_last_followup))] <- PAAD_meta$days_to_death[which(is.na(PAAD_meta$days_to_last_followup))]
PAAD_OS <- as.numeric(PAAD_OS)
PAAD_censor <- rep(NA, nrow(PAAD_meta))
PAAD_censor[which(PAAD_meta$vital_status == "alive")] <- 0
PAAD_censor[which(PAAD_meta$vital_status == "dead")] <- 1
PAAD_OS <- PAAD_OS[-c(13, 17, 79, 100)]
PAAD_censor <- PAAD_censor[-c(13, 17, 79, 100)]

decon_baron_surv <- decon_baron
decon_tosti_surv <- decon_tosti
decon_baron_surv$Yang$decon_res$prop.est.mvw <-decon_baron_surv$Yang$decon_res$prop.est.mvw[-c(32, 66, 67, 68),]
decon_tosti_surv$Yang$decon_res$prop.est.mvw <-decon_tosti_surv$Yang$decon_res$prop.est.mvw[-c(32, 66, 67, 68),]
decon_baron_surv$PAAD$decon_res$prop.est.mvw <-decon_baron_surv$PAAD$decon_res$prop.est.mvw[-c(13, 17, 79, 100),]
decon_tosti_surv$PAAD$decon_res$prop.est.mvw <-decon_tosti_surv$PAAD$decon_res$prop.est.mvw[-c(13, 17, 79, 100),]

baron_yang_survival <- survival_analysis(decon_output = decon_baron_surv$Yang, 
                                         OS = Yang_OS, censor = Yang_censor, 
                                         clinical_characteristics =  data.frame("grading" = Yang_meta$`grading:ch1`[-c(32, 66, 67, 68)], 
                                                                                "mki67" = as.numeric(Yang_bulk["MKI67",])[-c(32, 66, 67, 68)], 
                                                                                row.names = rownames(Yang_meta)[-c(32, 66, 67, 68)]))
tosti_yang_survival <- survival_analysis(decon_output = decon_tosti_surv$Yang, 
                                         OS = Yang_OS, censor = Yang_censor, 
                                         clinical_characteristics =  data.frame("grading" = Yang_meta$`grading:ch1`[-c(32, 66, 67, 68)], 
                                                                                "mki67" = as.numeric(Yang_bulk["MKI67",])[-c(32, 66, 67, 68)],
                                                                                row.names = rownames(Yang_meta)[-c(32, 66, 67, 68)]))


baron_PAAD_survival <- survival_analysis(decon_output = decon_baron_surv$PAAD, 
                                         OS = PAAD_OS, censor = PAAD_censor, cell_types = "acinar",
                                         clinical_characteristics =  data.frame(#"grading" = PAAD_meta$neoplasm_histologic_grade[-c(13, 17, 79, 100)], 
                                                                                "mki67" = as.numeric(PAAD_bulk["MKI67",])[-c(13, 17, 79, 100)],
                                                                                row.names = rownames(PAAD_meta)[-c(13, 17, 79, 100)]),
                                         legend.labs = c("baron_PAAD_acinar_high", "baron_PAAD_acinar_low", "MKi67_high", "MKi67_low"), xlab = "Time in days")

tosti_PAAD_survival <- survival_analysis(decon_output = decon_tosti_surv$PAAD, 
                                         OS = PAAD_OS, censor = PAAD_censor, 
                                         clinical_characteristics =  data.frame("grading" = PAAD_meta$neoplasm_histologic_grade[-c(13, 17, 79, 100)], 
                                                                                "mki67" = as.numeric(PAAD_bulk["MKI67",])[-c(13, 17, 79, 100)],
                                                                                row.names = rownames(PAAD_meta)[-c(13, 17, 79, 100)]))


##################
## correlation analysis
baron_yang_correlation <- correlation_analysis(decon_output = decon_baron$Yang,
                                               clinical_characteristic = Yang_meta$`grading:ch1`)
baron_yang_correlation2 <- correlation_analysis(decon_output = decon_baron$Yang,
                                               clinical_characteristic = as.numeric(Yang_bulk["MKI67",]))

tosti_yang_correlation <- correlation_analysis(decon_output = decon_tosti$Yang,
                                               clinical_characteristic = Yang_meta$`grading:ch1`)
tosti_yang_correlation2 <- correlation_analysis(decon_output = decon_tosti$Yang,
                                                clinical_characteristic = as.numeric(Yang_bulk["MKI67",]))

baron_PAAD_correlation <- correlation_analysis(decon_output = decon_baron$PAAD,
                                               clinical_characteristic = PAAD_meta$neoplasm_histologic_grade)
baron_PAAD_correlation2 <- correlation_analysis(decon_output = decon_baron$PAAD,
                                                clinical_characteristic = as.numeric(PAAD_bulk["MKI67",]))

tosti_PAAD_correlation <- correlation_analysis(decon_output = decon_tosti$PAAD,
                                               clinical_characteristic = PAAD_meta$neoplasm_histologic_grade)
tosti_PAAD_correlation2 <- correlation_analysis(decon_output = decon_tosti$PAAD,
                                                clinical_characteristic = as.numeric(PAAD_bulk["MKI67",]))


##################
## ML analysis
baron_yang_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_baron$Yang, clinical_char = Yang_meta$`grading:ch1`)
baron_yang_prepped <- baron_yang_prepped[-c(32, 66, 67, 68) ,] # remove G1, G4, Gx
baron_yang_prepped$response <- factor(baron_yang_prepped$response, levels = c("G2", "G3"))
baron_yang_trainRowNumbers <- createDataPartition(baron_yang_prepped$response, p = 0.8, list = FALSE)
baron_yang_train <- baron_yang_prepped[baron_yang_trainRowNumbers,]
baron_yang_test <- baron_yang_prepped[-baron_yang_trainRowNumbers,]
baron_yang_ml_model <- train_ML_model(trainData = baron_yang_train)
baron_yang_ml_pred <- test_ML_model(train_output = baron_yang_ml_model, 
                                    testData = baron_yang_test[,-ncol(baron_yang_test)], 
                                    truth_vec = baron_yang_test$response)
baron_yang_whole_ml_model <- train_ML_model(trainData = baron_yang_prepped)
## fasse die zusammen fuer jedes bulk zs mit vergleichsmodell (ki67)

baron_PAAD_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_baron$PAAD, clinical_char = PAAD_meta$neoplasm_histologic_grade)
baron_PAAD_prepped <- baron_PAAD_prepped[-c(13, 17, 79, 100) ,] # remove G4, Gx
baron_PAAD_prepped$response <- factor(toupper(baron_PAAD_prepped$response), levels = c("G1", "G2", "G3"))
baron_PAAD_trainRowNumbers <- createDataPartition(baron_PAAD_prepped$response, p = 0.8, list = FALSE)
baron_PAAD_train <- baron_PAAD_prepped[baron_PAAD_trainRowNumbers,]
baron_PAAD_test <- baron_PAAD_prepped[-baron_PAAD_trainRowNumbers,]
baron_PAAD_ml_model <- train_ML_model(trainData = baron_PAAD_train)
#baron_PAAD_ml_pred <- test_ML_model(train_output = baron_PAAD_ml_model, 
#                                    testData = baron_PAAD_test[,-ncol(baron_PAAD_test)], 
#                                    truth_vec = baron_PAAD_test$response)
baron_PAAD_prepped2 <- baron_PAAD_prepped
baron_PAAD_prepped2$response <- as.character(baron_PAAD_prepped2$response)
baron_PAAD_prepped2$response[baron_PAAD_prepped2$response == "G1"] <- "G1_G2"
baron_PAAD_prepped2$response[baron_PAAD_prepped2$response == "G2"] <- "G1_G2"
baron_PAAD_prepped2$response <- factor(baron_PAAD_prepped2$response, levels = c("G1_G2", "G3"))
baron_PAAD_whole_ml_model <- train_ML_model(trainData = baron_PAAD_prepped2)


tosti_yang_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$Yang, clinical_char = Yang_meta$`grading:ch1`)
tosti_yang_prepped <- tosti_yang_prepped[-c(32, 66, 67, 68) ,] # remove G1, G4, Gx
tosti_yang_prepped$response <- factor(tosti_yang_prepped$response, levels = c("G2", "G3"))
tosti_yang_trainRowNumbers <- createDataPartition(tosti_yang_prepped$response, p = 0.8, list = FALSE)
tosti_yang_train <- tosti_yang_prepped[tosti_yang_trainRowNumbers,]
tosti_yang_test <- tosti_yang_prepped[-tosti_yang_trainRowNumbers,]
tosti_yang_ml_model <- train_ML_model(trainData = tosti_yang_train)
tosti_yang_ml_pred <- test_ML_model(train_output = tosti_yang_ml_model, 
                                    testData = tosti_yang_test[,-ncol(tosti_yang_test)], 
                                    truth_vec = tosti_yang_test$response)
tosti_yang_whole_ml_model <- train_ML_model(trainData = tosti_yang_prepped)

tosti_PAAD_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$PAAD, clinical_char = PAAD_meta$neoplasm_histologic_grade)
tosti_PAAD_prepped <- tosti_PAAD_prepped[-c(13, 17, 79, 100) ,] # remove G4, Gx
tosti_PAAD_prepped$response <- factor(toupper(tosti_PAAD_prepped$response), levels = c("G1", "G2", "G3"))
tosti_PAAD_trainRowNumbers <- createDataPartition(tosti_PAAD_prepped$response, p = 0.8, list = FALSE)
tosti_PAAD_train <- tosti_PAAD_prepped[tosti_PAAD_trainRowNumbers,]
tosti_PAAD_test <- tosti_PAAD_prepped[-tosti_PAAD_trainRowNumbers,]
tosti_PAAD_ml_model <- train_ML_model(trainData = tosti_PAAD_train)
tosti_PAAD_ml_pred <- test_ML_model(train_output = tosti_PAAD_ml_model, 
                                    testData = tosti_PAAD_test[,-ncol(tosti_PAAD_test)], 
                                    truth_vec = tosti_PAAD_test$response)
tosti_PAAD_prepped2 <- tosti_PAAD_prepped
tosti_PAAD_prepped2$response <- as.character(tosti_PAAD_prepped2$response)
tosti_PAAD_prepped2$response[tosti_PAAD_prepped2$response == "G1"] <- "G1_G2"
tosti_PAAD_prepped2$response[tosti_PAAD_prepped2$response == "G2"] <- "G1_G2"
tosti_PAAD_prepped2$response <- factor(tosti_PAAD_prepped2$response, levels = c("G1_G2", "G3"))
tosti_PAAD_whole_ml_model <- train_ML_model(trainData = tosti_PAAD_prepped2)


##################
## visualize ML results
baron_yang_roc <- roc_curve(labels = baron_yang_test$response, 
                            predictions = baron_yang_ml_pred$predicted_reduced, levels = c("G2", "G3"))
#baron_PAAD_roc <- roc_curve(labels = baron_PAAD_test$response, 
#                            predictions = baron_PAAD_ml_pred$predicted_reduced, levels = c("G1", "G2", "G3"))
tosti_yang_roc <- roc_curve(labels = tosti_yang_test$response, 
                            predictions = tosti_yang_ml_pred$predicted_reduced, levels = c("G2", "G3"))
tosti_PAAD_roc <- roc_curve(labels = tosti_PAAD_test$response, 
                            predictions = tosti_PAAD_ml_pred$predicted_reduced, levels = c("G1", "G2", "G3"))

ml_model_eval_yang <- list("baron_yang" = baron_yang_ml_pred$evaluation_reduced,
                           "tosti_yang" = tosti_yang_ml_pred$evaluation_reduced)
ml_model_eval_PAAD <- list(#"baron_PAAD" = baron_PAAD_ml_pred$evaluation_reduced,
                           "tosti_PAAD" = tosti_PAAD_ml_pred$evaluation_reduced)
ml_model_eval_yang_plot <- barplot_ML_evaluation(ml_model_eval_yang)
ml_model_eval_PAAD_plot <- barplot_ML_evaluation(ml_model_eval_PAAD)


ml_model_yang <- list(#"baron_yang_whole" = baron_yang_whole_ml_model$rf_model_whole,
                      "baron_yang" = baron_yang_whole_ml_model$rf_model_reduced,
                      #"tosti_yang_whole" = tosti_yang_whole_ml_model$rf_model_whole,
                      "tosti_yang" = tosti_yang_whole_ml_model$rf_model_reduced)
ml_model_PAAD <- list(#"baron_PAAD_whole" = baron_PAAD_whole_ml_model$rf_model_whole,
                      "baron_PAAD" = baron_PAAD_whole_ml_model$rf_model_reduced,
                      #"tosti_PAAD_whole" = tosti_PAAD_whole_ml_model$rf_model_whole,
                      "tosti_PAAD" = tosti_PAAD_whole_ml_model$rf_model_reduced)
ml_model_yang_plot <- boxplot_ML_sd(ml_model_yang)
ml_model_PAAD_plot <- boxplot_ML_sd(ml_model_PAAD)
boxplot_ML_sd(c(ml_model_yang, ml_model_PAAD))

######################
## Mki-67 comparison model
fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  sampling = "down",
  savePred = TRUE
)

yang_mki67_df <- data.frame("mki67" = as.numeric(Yang_bulk["MKI67",]),
                            "response" = Yang_meta$`grading:ch1`, 
                            row.names = rownames(Yang_meta))
yang_mki67_df <- yang_mki67_df[-c(32, 66, 67, 68) ,] # remove G1, G4, Gx
yang_mki67_df$response <- factor(yang_mki67_df$response, levels = c("G2", "G3"))
yang_mki67_model <- train(x = data.frame("mki67" = yang_mki67_df$mki67, row.names = rownames(yang_mki67_df)), 
                          y = yang_mki67_df$response, 
                          method = "rf", metric = "Accuracy", trControl = fitControl,
                          type = "Classification", ntree = 500)
PAAD_mki67_df <- data.frame("mki67" = as.numeric(PAAD_bulk["MKI67",]),
                            "response" = PAAD_meta$neoplasm_histologic_grade, 
                            row.names = rownames(PAAD_meta))
PAAD_mki67_df <- PAAD_mki67_df[-c(13, 17, 79, 100) ,] # remove G4, Gx
PAAD_mki67_df$response <- toupper(PAAD_mki67_df$response)
PAAD_mki67_df$response[PAAD_mki67_df$response == "G1"] <- "G1_G2"
PAAD_mki67_df$response[PAAD_mki67_df$response == "G2"] <- "G1_G2"
PAAD_mki67_df$response <- factor(PAAD_mki67_df$response, levels = c("G1_G2", "G3"))
PAAD_mki67_model <- train(x = data.frame("mki67" = PAAD_mki67_df$mki67, row.names = rownames(PAAD_mki67_df)), 
                          y = PAAD_mki67_df$response, 
                          method = "rf", metric = "Accuracy", trControl = fitControl,
                          type = "Classification", ntree = 500)
mki67_ct_models <- list("baron_yang" = baron_yang_whole_ml_model$rf_model_reduced,
                        "tosti_yang" = tosti_yang_whole_ml_model$rf_model_reduced,
                        "mki67_yang" = yang_mki67_model)
boxplot_ML_sd(mki67_ct_models)
mki67_ct_models2 <- list("baron_PAAD" = baron_PAAD_whole_ml_model$rf_model_reduced,
                        "tosti_PAAD" = tosti_PAAD_whole_ml_model$rf_model_reduced,
                        "mki67_PAAD" = PAAD_mki67_model )
boxplot_ML_sd(mki67_ct_models2)
boxplot_ML_sd(c(mki67_ct_models, mki67_ct_models2))




############################
basal_classical_signature <- data.frame(t(read.table("~/Masterthesis/Data/Bulk/Moffitt/basal_classical_signature.txt",
                                                     header = FALSE, row.names = 1, sep = "\t")))
## import Guo signature
guo_signature <- data.frame(read.table("~/Masterthesis/Data/Bulk/Guo/signature_genes.txt",
                                       header = FALSE, sep = "\t"))
## import supplementary file from Hayashi et al. 2020
hayashi_PAAD_meta <- read.table("~/Masterthesis/Data/Bulk/PAAD/hayashi_suppl.txt", 
                                sep = "\t", header = TRUE)
hayashi_PAAD_meta$TCGA_ID <- hayashi_PAAD_meta$Tumor.Sample.ID
hayashi_PAAD_meta$TCGA_ID <- gsub("-", ".", hayashi_PAAD_meta$TCGA_ID)
hayashi_PAAD_meta$TCGA_ID <- sapply(hayashi_PAAD_meta$TCGA_ID, function(x) substr(x, 1, nchar(x)-4))
rownames(hayashi_PAAD_meta) <- hayashi_PAAD_meta$TCGA_ID
hayashi_idx <- match(rownames(PAAD_meta), rownames(hayashi_PAAD_meta))
hayashi_PAAD_meta$tumor_subtype <- hayashi_PAAD_meta$mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical
hayashi_PAAD_meta$tumor_subtype[hayashi_PAAD_meta$tumor_subtype == 1] <- "Basal"
hayashi_PAAD_meta$tumor_subtype[hayashi_PAAD_meta$tumor_subtype == 2] <- "Classical"
hayashi_PAAD_meta$tumor_collisson <- hayashi_PAAD_meta$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM
hayashi_PAAD_meta$tumor_collisson[hayashi_PAAD_meta$tumor_collisson == 1] <- "Classical"
hayashi_PAAD_meta$tumor_collisson[hayashi_PAAD_meta$tumor_collisson == 2] <- "Exocrine"
hayashi_PAAD_meta$tumor_collisson[hayashi_PAAD_meta$tumor_collisson == 3] <- "QM"
hayashi_PAAD_meta$tumor_bailey <- hayashi_PAAD_meta$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 1] <- "Squamous"
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 2] <- "Immunogenic"
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 3] <- "Progenitor"
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 4] <- "ADEX"

PAAD_meta$tumor_subtype <- hayashi_PAAD_meta$tumor_subtype[hayashi_idx]
PAAD_meta$tumor_collisson <- hayashi_PAAD_meta$tumor_collisson[hayashi_idx]
PAAD_meta$tumor_bailey <- hayashi_PAAD_meta$tumor_bailey[hayashi_idx]

PAAD_mki67 <- continuous_to_discrete(PAAD_mki67_df$mki67, "MKI67")
PAAD_mki67_thirds <- quantile(PAAD_mki67_df$mki67, probs = seq(0, 1, 1/3))
PAAD_mki67.2 <- as.numeric(PAAD_mki67_df$mki67)
PAAD_mki67.2[sapply(as.numeric(PAAD_mki67_df$mki67), function(x) x <= PAAD_mki67_thirds[2])] <- "MKI67_low"
PAAD_mki67.2[sapply(as.numeric(PAAD_mki67_df$mki67), function(x) x > PAAD_mki67_thirds[2] && x <= PAAD_mki67_thirds[3])] <- "MKI67_medium"
PAAD_mki67.2[sapply(as.numeric(PAAD_mki67_df$mki67), function(x) x > PAAD_mki67_thirds[3])] <- "MKI67_high"

annot_colors <- list(collisson = c(Exocrine = "#f683ad", Classical = "#f8fc88", QM ="#6eacf2"),
                     bailey = c(ADEX = "#77f387", Immunogenic = "#f19e5b", Progenitor ="#6eacf2", Squamous = "#f8fc88"),
                     moffitt = c(Basal = "#77f387", Classical = "#f19e5b"),
                     grading = c(g1 = "#f683ad", g2 = "#77f387", g3 = "#6eacf2"),
                     MKI67 = c(MKI67_low = "#77f387", MKI67_medium ="#f19e5b",MKI67_high = "#d34545")) 
tosti_PAAD_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti_surv$PAAD,
                                               clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade[-c(13, 17, 79, 100)],
                                                                                     "moffitt" = PAAD_meta$tumor_subtype[-c(13, 17, 79, 100)],
                                                                                     "collisson" = PAAD_meta$tumor_collisson[-c(13, 17, 79, 100)],
                                                                                     "bailey" = PAAD_meta$tumor_bailey[-c(13, 17, 79, 100)],
                                                                                     "MKI67" = PAAD_mki67.2,
                                                                                     row.names = rownames(PAAD_meta)[-c(13, 17, 79, 100)]), 
                                               clustering_method = "ward.D2", annotation_colors = annot_colors)


tosti_PAAD_subtype_heatmap <- heatmap_corr_genes(decon_tosti$PAAD, cell_types = c("sacinar", "mductal"),
                                                 bulk_data = PAAD_bulk, 
                                                 clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                       "moffitt" = PAAD_meta$tumor_subtype,
                                                                                       "collisson" = PAAD_meta$tumor_collisson,
                                                                                       "bailey" = PAAD_meta$tumor_bailey,
                                                                                       row.names = rownames(PAAD_meta)),
                                                 marker_genes = c(basal_classical_signature$Basal, basal_classical_signature$Classical))
# heatmap_label_order <- rownames(PAAD_meta)[tosti_PAAD_subtype_heatmap$tree_col$order]
# which(heatmap_label_order == "TCGA.XD.AAUG") #79 until here all are classical
# which(heatmap_label_order == "TCGA.2J.AABH") #149 until here all are hybrid
# which(heatmap_label_order == "TCGA.IB.A7LX") #158 from here until
# which(heatmap_label_order == "TCGA.XN.A8T3") #183 here all are basal
# heatmap_label_order_type <- heatmap_label_order
# heatmap_label_order_type[1:79] <- "Classical"
# heatmap_label_order_type[80:149] <- "Hybrid"
# heatmap_label_order_type[150:157] <- NA
# heatmap_label_order_type[158:183] <- "Basal"
# PAAD_meta$tumor_subtype2 <- heatmap_label_order_type

tosti_PAAD_correlation3 <- correlation_analysis(decon_output = decon_tosti$PAAD, clinical_characteristic = PAAD_meta$tumor_subtype)
tosti_PAAD_correlation4 <- correlation_analysis(decon_output = decon_tosti$PAAD, clinical_characteristic = PAAD_meta$tumor_collisson)
tosti_PAAD_correlation4 <- correlation_analysis(decon_output = decon_tosti$PAAD, clinical_characteristic = PAAD_meta$tumor_bailey)
tosti_PAAD_survival2 <- survival_analysis(decon_output = decon_tosti_surv$PAAD, OS = PAAD_OS, censor = PAAD_censor, 
                                          clinical_characteristics =  data.frame("grading" = PAAD_meta$neoplasm_histologic_grade[-c(13, 17, 79, 100)], 
                                                                                 "moffitt" = PAAD_meta$tumor_subtype[-c(13, 17, 79, 100)],
                                                                                 "collisson" = PAAD_meta$tumor_collisson[-c(13, 17, 79, 100)],
                                                                                 "bailey" = PAAD_meta$tumor_bailey[-c(13, 17, 79, 100)],
                                                                                 "MKI67" = PAAD_mki67.2,
                                                                                 row.names = rownames(PAAD_meta)[-c(13, 17, 79, 100)]))
