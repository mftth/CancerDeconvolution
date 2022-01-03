## Mastherthesis, Melanie Fattohi
## test functions of Execute_MachineLearning.R

source("~/Masterthesis/CancerDeconvolution/Scripts/Quality_control.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_Deconvolution.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_MachineLearning.R")

library(parallel)
library(robustbase)
library(caret)
library(skimr)

set.seed(4)

#Repset_scdc_baron <- readRDS("~/Masterthesis/Deconvolution_results/Repset_scdc_baron.RDS")
repset <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")
#all(rownames(repset_meta) == rownames(Repset_scdc_baron$decon_res$prop.est.mvw))

qc_baron <- readRDS(file = "~/Praktikum/Data/Baron/qc_baron_exo.RDS")

## create basis
# baron_basis <- Create_basis(multiple_donors = TRUE, sc_data = qc_baron$sc.eset.qc, 
#                             cell_types = c("alpha", "beta", "gamma", "delta", "acinar", "ductal"))
# 
# ## deconvolve once
# repset_baron_scdc <- Deconvolve_SCDC(bulk_data = repset, bulk_meta = repset_meta, sc_data = qc_baron$sc.eset.qc,
#                                      sc_basis = baron_basis, cell_types = c("alpha", "beta", "gamma", "delta", "acinar", "ductal"),
#                                      ensemble = FALSE, multiple_donors = TRUE)

## perform deconvolution and calc p-value
repset_decon <- Calculate_pvalue(nrep = 1000, ncores = 10, silent = TRUE, bulk_data = repset,
                                 sc_data = qc_baron$sc.eset.qc, 
                                 cell_types =  c("alpha", "beta", "gamma", "delta", "acinar", "ductal"),
                                 ensemble = FALSE, multiple_donors = TRUE, bulk_meta = repset_meta)

prepped_repset_decon <- prepare_decon_res(decon_res = repset_decon, clinical_char = repset_meta$Grading, p_value = TRUE)

#decon_res <- prepare_decon_res(decon_res = Repset_scdc_baron, clinical_char = repset_meta$Grading, p_value = TRUE)
trainRowNumbers <- createDataPartition(prepped_repset_decon$response, p = 0.8, list = FALSE)
trainData <- prepped_repset_decon[trainRowNumbers,]
testData <- prepped_repset_decon[-trainRowNumbers,]
train_output <- train_ML_model(trainData = trainData, preprocess = TRUE)
#rfProfile
#predictors(rfProfile)
#plot(varimp_rf)
#plot(varimp_rf_sel)
test_output <- test_ML_model(train_output = train_output, testData = testData[, -ncol(testData)], 
                             true_labels = testData$response, preprocess = TRUE)




library(plotROC)
selectedIndices <- train_output$rf_model_whole$pred$mtry == train_output$rf_model_whole$bestTune[[1]]
g <- ggplot(train_output$rf_model_whole$pred[selectedIndices, ], aes(m=pred, d=factor(obs, levels = c("R", "M")))) +  # bei levels gib die klassen deiner response variable an
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc()

g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))











#x = trainData[, - ncol(trainData)]
#y = trainData$response

## analyze feature importance
# featurePlot(x = trainData[, - ncol(trainData)], 
#             y = factor(trainData$response), 
#             plot = "box",
#             strip=strip.custom(par.strip.text=list(cex=.7)),
#             scales = list(x = list(relation="free"), 
#                           y = list(relation="free")))
# 
# featurePlot(x = trainData[, - ncol(trainData)], 
#             y = factor(trainData$response), 
#             plot = "density",
#             strip=strip.custom(par.strip.text=list(cex=.7)),
#             scales = list(x = list(relation="free"), 
#                           y = list(relation="free")))

## ensemble predictions
# library(caretEnsemble)
# trainControl <- trainControl(method="repeatedcv", 
#                             number=10, 
#                             repeats=3,
#                             savePredictions=TRUE, 
#                             classProbs=TRUE)
# algorithmList <- c('rf', 'svmRadial')
# models <- caretList(response~ ., data=trainData, trControl=trainControl, methodList=algorithmList) 
# results <- resamples(models)
# summary(results)