## Mastherthesis, Melanie Fattohi
## test previously trained Random Forest model regarding predicting tumor grading
## model saved in /home/fattohim/Praktikum/ML_RF_SC2.RData
## bulk RNA-seq data to be tested saved in /home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/*
## evaluate prediction by means of accuracy, sensitivity and specificity
## perform feature importance 


library(dplyr)
library(pROC)
library(caret)
library(stringr)
library(e1071)
library(SCDC)
set.seed(1)

## load workspace of internship containing the RF ML model
load("/home/fattohim/Masterthesis/Workspaces/Prak_RF_ML_model.RData")

## retrain the RF models for exocrine model
## (models were not saved in internship)  
## repeat for all 4 scRNA-seq references
fitControl <- trainControl(
  method = "cv",
  number = 10,
  sampling = "down",
  savePred=T
)

res_exo_new <- list()
d_exo_new <- list()
fImp_exo_new <- list()
exo_models <- list() ## will contain the models

for (i_alg in 1:length(cell_type_prop_exocrine_repset)) {
  train_mat <- cell_type_prop_exocrine_repset[[i_alg]]
  colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "RMSD", "mAD", "spearman")
  #truth_vec <- repset_meta$Grading
  
  model_exo <- caret::train(
    x = train_mat,
    y = truth_vec,
    method = "rf",
    norm.votes = T,
    type = "Classification",
    metric = "Accuracy",
    ntree = 500,
    trControl = fitControl
  )
  
  exo_models[[i_alg]] <- model_exo
  #truth_vec = factor(truth_vec, levels = c("G1","G2","G3"))
  prediction_ml = predict(model_exo, train_mat)
  con_mat_exo = confusionMatrix(prediction_ml, truth_vec, positive = "G3")
  
  res_exo_new[[i_alg]] <- con_mat_exo
  d_exo_new[[i_alg]] <- res_exo[[i_alg]]$byClass
  
  # feature importance
  fImp_exo_new[[i_alg]] <- varImp(model_exo, scale = T) 
  
}

names(exo_models) <- names(cell_type_prop_exocrine_repset)
names(fImp_exo_new) <- names(cell_type_prop_exocrine_repset)
names(res_exo_new) <- names(cell_type_prop_exocrine_repset)
names(d_exo_new) <- names(cell_type_prop_exocrine_repset)


## write function for testing of all four ML models
## input: bulk RNA-seq data, grading
## output: evaluation of prediction performance (confusion_matrix)

test_int_ML_model <- function(decon_bulk, grading_vec, ensemble) {
	if(ensemble) {
	  mean_weights <- sapply(1:3, function(x) mean(Yang_ensemble$w_table[1:5,1:3][,x]))
	  decon_bulk_prop <- Yang_ensemble$prop.only[[which.max(mean_weights)]]
	  train_mat <- cbind(
		       decon_bulk_prop,
		       t(decon_bulk$prop.list[[which.max(mean_weights)]]$yeval$RMSDy.sample.table),
		       t(decon_bulk$prop.list[[which.max(mean_weights)]]$yeval$mADy.sample.table),
		       t(decon_bulk$prop.list[[which.max(mean_weights)]]$yeval$spearmany.sample.table))
	  colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "RMSD", "mAD", "spearman")
	  
	} else { 
	   train_mat <- cbind(
			decon_bulk$prop.est.mvw, 
			t(decon_bulk$yeval$RMSDy.sample.table),
                        t(decon_bulk$yeval$mADy.sample.table),
                        t(decon_bulk$yeval$spearmany.sample.table))
           colnames(train_mat) <- c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "RMSD", "mAD", "spearman")
	}
	
        prediction_ml <- sapply(names(cell_type_prop_exocrine_repset), function(x) predict(exo_models[[x]], train_mat))
        confusion_matrix <- apply(prediction_ml, 2,
			    function(x) confusionMatrix(factor(x, levels = c("G1", "G2", "G3")), 
			                grading_vec, positive = "G3"))
        
        
        
        RF_performance <- matrix(data = NA, nrow = length(confusion_matrix)*nrow(confusion_matrix[[1]]$byClass), 
			         ncol = 2+ncol(confusion_matrix[[1]]$byClass))
	colnames(RF_performance) <- c("SCRef_Model", "Class", colnames(confusion_matrix[[1]]$byClass))
	RF_performance[,1] <- rep(names(confusion_matrix), each = nrow(confusion_matrix[[1]]$byClass))
	RF_performance[,2] <- rep(c("G1", "G2", "G3"), length(confusion_matrix))

	j = 1
	for (i in 1:length(confusion_matrix)) {
		RF_performance[j:(j+2), 3:ncol(RF_performance)] <- confusion_matrix[[i]]$byClass
		j = j+3
	}
        
        
        return(list("confusion_matrix" = confusion_matrix, summary = RF_performance))     
}



## Yang, 2016, GSE62452, microarray, 69 samples
## read bulk data +  meta data of Yang
Yang_bulk <- readRDS(file = "/home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/Yang/Yang_bulk.RDS")
Yang_metadata <- readRDS(file = "/home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/Yang/Yang_metadata.RDS")

## perform deconvolution of Yang with SCDC(Baron)/ENSEMBLE 
## remove Gx and G4 samples
Yang_bulk <- Yang_bulk[,-which(Yang_metadata$'grading:ch1'== 'Gx')]
Yang_bulk <- Yang_bulk[,-which(Yang_metadata$'grading:ch1'== 'G4')]
Yang_metadata <- Yang_metadata[-which(Yang_metadata$'grading:ch1'== 'Gx'),]
Yang_metadata <- Yang_metadata[-which(Yang_metadata$'grading:ch1'== 'G4'),]
all(colnames(Yang_bulk) == rownames(Yang_metadata))
## create Expression Set
fdata_Yang <- rownames(Yang_bulk)
pdata_Yang <- cbind(cellname = rownames(Yang_metadata), subjects = Yang_metadata$'grading:ch1')
eset_bulk <- getESET(Yang_bulk, fdata = fdata_Yang, pdata = pdata_Yang)
## SCDC deconvolution
qc_baron_exo <- readRDS("/home/fattohim/Praktikum/Deko_SCDC/Training_Data/Baron_qc_exo.RDS")

Yang_scdc_baron <- SCDC_prop(bulk.eset = eset_bulk, sc.eset = qc_baron_exo$sc.eset.qc, 
                                      ct.varname = "cluster", sample = "sample", 
                                      ct.sub = c("alpha","beta","delta","gamma", "acinar", "ductal"))                               
## ENSEMBLE deconvolution
qc_lawlor_exo <- readRDS("/home/fattohim/Praktikum/Deko_SCDC/Training_Data/Lawlor_qc_exo.RDS")
qc_segerstolpe_exo <- readRDS("/home/fattohim/Praktikum/Deko_SCDC/Training_Data/Segerstolpe_qc_exo.RDS")
pancreas_ref_exo <- list(baron = qc_baron_exo$sc.eset.qc,
                          seger = qc_segerstolpe_exo$sc.eset.qc,
                          lawlor = qc_lawlor_exo$sc.eset.qc)
Yang_ensemble <- SCDC_ENSEMBLE(bulk.eset = eset_bulk, sc.eset.list = pancreas_ref_exo, ct.varname = "cluster",
                                    sample = "sample", ct.sub =  c(#"alpha","beta","delta","gamma", 
                                    "acinar", "ductal"), search.length = 0.01)
Yang_ensemble$prop.only$lawlor <- cbind("alpha" = numeric(67), "beta" = numeric(67), "delta" = numeric(67), 
                                        "gamma" = numeric(67), Yang_ensemble$prop.only$lawlor)
Yang_ensemble$prop.only$seger <- cbind("alpha" = numeric(67), "beta" = numeric(67), "delta" = numeric(67), 
                                        "gamma" = numeric(67), Yang_ensemble$prop.only$seger)

## apply test_int_ML_model() to the data
Yang_scdc_baron_pred <- test_int_ML_model(decon_bulk = Yang_scdc_baron, 
                                          grading_vec = eset_bulk@phenoData$subjects, 
                                          ensemble = FALSE)
Yang_ensemble_pred <- test_int_ML_model(decon_bulk = Yang_ensemble, 
                                        grading_vec = eset_bulk@phenoData$subjects, 
                                        ensemble = TRUE)

## save results
write.table(Yang_scdc_baron_pred$summary, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Results/Prak_grading_prediction/Yang_scdc_baron_pred.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(Yang_scdc_baron_pred, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Results/Prak_grading_prediction/Yang_scdc_baron_pred.RDS")

write.table(Yang_ensemble_pred$summary, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Results/Prak_grading_prediction/Yang_ensemble_pred.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(Yang_ensemble_pred, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Results/Prak_grading_prediction/Yang_ensemble_pred.RDS")


