## Mastherthesis, Melanie Fattohi
## machine learning
## input
### decon_res +  reconstruction error + pvalue
### meta data / clinical data of bulk RNA-seq
##output
### ML model trained on decon_res + recon error + pval (optional)
### prediction of a clinical characteristic of bulk samples based on a newly trained or existing ML model
### evaluation of prediction

library(caret)
library(skimr)

# Repset_scdc_baron <- readRDS("~/Masterthesis/Deconvolution_results/Repset_scdc_baron.RDS")
# 
# repset <- read.table(file = "~/Praktikum/Data/RepSet.S57.tsv", header = TRUE, 
#                      sep = "\t",)
# pannen_meta <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", 
#                           header = TRUE, sep = "\t")
# scarpa_meta <- pannen_meta[which(pannen_meta$Study == "Scarpa"),]
# riemer_meta <- pannen_meta[which(pannen_meta$Study == "Riemer"),]
# riemer_meta <- riemer_meta[which(riemer_meta$Subtype == "Cancer"),]
# repset_meta <- rbind(scarpa_meta, riemer_meta)
# repset_meta$Name <- as.character(repset_meta$Name)
# repset_meta$Name[1:55] <- paste("X", repset_meta$Name[1:55], sep = "")
# repset_meta <- repset_meta[match(colnames(repset), repset_meta$Name),]
# all(repset_meta$Name == colnames(repset))
# rownames(repset_meta) <- repset_meta$Name
# all(rownames(repset_meta) == rownames(Repset_scdc_baron$decon_res$prop.est.mvw))



train_ML_model <- function(trainData, preprocess = FALSE, preprocess_method = "scale"){
  ## descriptive statistics (necessary? kinda)
  skimmed <- skim(trainData)
  
  ## impute missing values if there are any
  if (any(skimmed$n_missing > 0)){
    preProcess_missingdata_model <- preProcess(trainData, method = "knnImpute")
    trainData <- predict(preProcess_missingdata_model, newdata = trainData)
  }
  
  ## preprocess data --> which method?
  if (preprocess){
    preProcess_method_model <- preProcess(trainData, method = preprocess_method)
    trainData <- predict(preProcess_method_model, newdata = trainData)
  }
  
  ## feature selection
  subsets <- c(1:8)
  ctrl <- rfeControl(functions = rfFuncs,
                     method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     verbose = FALSE)
  
  rfProfile <- rfe(x = trainData[,-ncol(trainData)], 
                   y = trainData$response,
                   sizes = subsets,
                   rfeControl = ctrl,
                   metric = "Accuracy",)
  #rfProfile
  #predictors(rfProfile)
  
  ## train RF model
  fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    sampling = "down",
    savePred=T
  )
  
  model_rf <- train(x = trainData[, - ncol(trainData)], 
                    y = trainData$response, 
                    method = 'rf', 
                    metric = "Accuracy",
                    trControl = fitControl, 
                    type = "Classification",
                    ntree = 500)
  #fitted <- predict(model_rf)
  varimp_rf <- varImp(model_rf)
  #plot(varimp_rf)
  
  ## train on selected features
  trainData_sel <- trainData[, predictors(rfProfile)[1:5]]
  model_rf_sel <- train(x = trainData_sel, 
                        y = trainData$response, 
                        method = 'rf', 
                        metric = "Accuracy",
                        trControl = fitControl, 
                        type = "Classification",
                        ntree = 500)
  fitted_sel <- predict(model_rf_sel)
  varimp_rf_sel <- varImp(model_rf_sel)
  #plot(varimp_rf_sel)
  
  train_output <- list("rf_model_whole" = model_rf, "rf_model_sel" = model_rf_sel, 
                       "varimp_whole" = varimp_rf, "varimp_sel" = varimp_rf_sel)
  return(train_output)
}


test_ML_model <- function(train_output, testData, true_labels,  preprocess = FALSE, preprocess_method = "scale"){ ## test data, aka hold out set
  ## descriptive statistics
  skimmed <- skim(testData)
  
  ## impute missing values if there are any
  if (any(skimmed$n_missing > 0)){
    preProcess_missingdata_model <- preProcess(testData, method = "knnImpute")
    testData <- predict(preProcess_missingdata_model, newdata = testData)
  }
  
  ## preprocess data 
  if (preprocess){
    preProcess_method_model <- preProcess(testData, method = preprocess_method)
    testData <- predict(preProcess_method_model, newdata = testData)
  }
  
  predicted_whole <- predict(train_output$rf_model_whole, testData)
  predicted_sel <- predict(train_output$rf_model_sel, testData)
  con_mat_whole <- caret::confusionMatrix(data = predicted_whole, reference = true_labels, mode = "everything")
  con_mat_sel <- caret::confusionMatrix(data = predicted_sel, reference = true_labels, mode = "everything")
  
  test_output <- list("predicted_whole" = predicted_whole, "con_mat_whole" = con_mat_whole,
                      "predicted_sel" = predicted_sel, "con_mat_sel" = con_mat_sel)
  return(test_output)
}




## split into training and test set (make it optional?)
set.seed(4)
decon_res <- prepare_decon_res(decon_res = Repset_scdc_baron, clinical_char = repset_meta$Grading, p_value = TRUE)
trainRowNumbers <- createDataPartition(decon_res$response, p = 0.8, list = FALSE)
trainData <- decon_res[trainRowNumbers,]
testData <- decon_res[-trainRowNumbers,]
train_output <- train_ML_model(trainData = trainData, preprocess = TRUE)
test_output <- test_ML_model(train_output = train_output, testData = testData[, -ncol(testData)], true_labels = testData$response, preprocess = TRUE)

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
