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


train_ML_model <- function(trainData, preprocess = FALSE, preprocess_method = "scale"){
  ## descriptive statistics (necessary? kinda)
  skimmed <- skim(trainData)
  
  ## impute missing values if there are any
  if (any(skimmed$n_missing > 0)){
    message("Imputing missing values ..")
    preProcess_missingdata_model <- preProcess(trainData, method = "knnImpute")
    trainData <- predict(preProcess_missingdata_model, newdata = trainData)
  }
  
  ## preprocess data 
  if (preprocess){
    message("Preprocessing data ..")
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
  
  message("Performing feature selection ..")
  rfProfile <- rfe(x = trainData[,-ncol(trainData)], 
                   y = trainData$response,
                   sizes = subsets,
                   rfeControl = ctrl,
                   metric = "Accuracy",)
  
  fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    sampling = "down", # not specified for regression
    savePred = TRUE
  )
  
  ## train RF model on all features
  message("Training a model on all features")
  model_rf <- train(x = trainData[, - ncol(trainData)], 
                    y = trainData$response, 
                    method = 'rf', 
                    metric = "Accuracy", # RMSE for regression
                    trControl = fitControl, 
                    type = "Classification",
                    ntree = 500)
  varimp_rf <- varImp(model_rf)
  
  ## train on RF model selected features
  trainData_sel <- trainData[, predictors(rfProfile)[1:5]]
  message("Training a model on selected features")
  model_rf_sel <- train(x = trainData_sel, 
                        y = trainData$response, 
                        method = 'rf', 
                        metric = "Accuracy",
                        trControl = fitControl, 
                        type = "Classification",
                        ntree = 500)
  varimp_rf_sel <- varImp(model_rf_sel)
  
  train_output <- list("rf_model_whole" = model_rf, "rf_model_sel" = model_rf_sel, 
                       "varimp_whole" = varimp_rf, "varimp_sel" = varimp_rf_sel)
  return(train_output)
}


test_ML_model <- function(train_output, testData, true_labels,  preprocess = FALSE, preprocess_method = "scale"){ ## test data, aka hold out set
  ## descriptive statistics
  skimmed <- skim(testData)
  
  ## impute missing values if there are any
  if (any(skimmed$n_missing > 0)){
    message("Imputing missing values ..")
    preProcess_missingdata_model <- preProcess(testData, method = "knnImpute")
    testData <- predict(preProcess_missingdata_model, newdata = testData)
  }
  
  ## preprocess data 
  if (preprocess){
    message("Preprocessing data ..")
    preProcess_method_model <- preProcess(testData, method = preprocess_method)
    testData <- predict(preProcess_method_model, newdata = testData)
  }
  
  message("Applying models to test data ..")
  predicted_whole <- predict(train_output$rf_model_whole, testData)
  predicted_sel <- predict(train_output$rf_model_sel, testData)
  con_mat_whole <- caret::confusionMatrix(data = predicted_whole, reference = true_labels, mode = "everything")
  con_mat_sel <- caret::confusionMatrix(data = predicted_sel, reference = true_labels, mode = "everything")
  
  test_output <- list("predicted_whole" = predicted_whole, "con_mat_whole" = con_mat_whole,
                      "predicted_sel" = predicted_sel, "con_mat_sel" = con_mat_sel)
  return(test_output)
}




