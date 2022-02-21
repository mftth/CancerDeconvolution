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


train_ML_model <- function(trainData, preprocess = FALSE, preprocess_method = "scale", classification = TRUE){
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
  
  if(classification){
    my_metric <- "Accuracy"
    my_type <- "Classification"
    my_sampling <- "down"
  } else{
    my_metric <- "RMSE"
    my_type <- "Regression"
    my_sampling <- NULL
  }
  
  ## feature selection
  subsets <- c(1:8)
  # creates a control object, i.e. a list of options for specifying the model and method
  ctrl <- rfeControl(functions = rfFuncs, # rfFuncs = random forest
                     method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     verbose = FALSE)
  
  message("Performing feature selection ..")
  rfProfile <- rfe(x = trainData[,-ncol(trainData)], # rfe is feature selection algorithm
                   y = trainData$response,
                   sizes = subsets,
                   rfeControl = ctrl,
                   metric = my_metric)
  
  fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    sampling = my_sampling,
    savePred = TRUE
  )
  
  ## train RF model on all features
  message("Training a model on all features ..")
  model_rf <- train(x = trainData[, - ncol(trainData)], 
                    y = trainData$response, 
                    method = 'rf', 
                    metric = my_metric,
                    trControl = fitControl, 
                    type = my_type,
                    ntree = 500)
  varimp_rf <- varImp(model_rf)
  
  ## train on RF model selected features
  if(length(predictors(rfProfile)) == length(trainData[, - ncol(trainData)])){
    predictors_red <- predictors(rfProfile)[1:5]
  } else {
    predictors_red <- predictors(rfProfile)
  }
  trainData_red <- trainData[, predictors_red]
  message("Training a model on selected features ..")
  model_rf_red <- train(x = trainData_red, 
                        y = trainData$response, 
                        method = 'rf', 
                        metric = my_metric,
                        trControl = fitControl, 
                        type = my_type,
                        ntree = 500)
  varimp_rf_red <- varImp(model_rf_red)
  
  train_output <- list("rf_model_whole" = model_rf, "rf_model_reduced" = model_rf_red, 
                       "varimp_whole" = varimp_rf, "varimp_reduced" = varimp_rf_red)
  return(train_output)
}


test_ML_model <- function(train_output, testData, truth_vec,  preprocess = FALSE, preprocess_method = "scale",
                          classification = TRUE){ ## test data, aka hold out set
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
  predicted_red <- predict(train_output$rf_model_reduced, testData)
  
  if(classification){
    evaluation_whole <- caret::confusionMatrix(data = predicted_whole, reference = truth_vec, mode = "everything")
    evaluation_red <- caret::confusionMatrix(data = predicted_red, reference = truth_vec, mode = "everything")
  } else {
    evaluation_whole <- postResample(pred = predicted_whole, obs = truth_vec)
    evaluation_red <- postResample(pred = predicted_red, obs = truth_vec)  
  }
  
  test_output <- list("predicted_whole" = predicted_whole, "evaluation_whole" = evaluation_whole,
                      "predicted_reduced" = predicted_red, "evaluation_reduced" = evaluation_red)
  return(test_output)
}




