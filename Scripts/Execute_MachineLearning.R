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
library(RANN)

Repset_scdc_baron <- readRDS("~/Masterthesis/Deconvolution_results/Repset_scdc_baron.RDS")

repset <- read.table(file = "~/Praktikum/Data/RepSet.S57.tsv", header = TRUE, 
                     sep = "\t",)
pannen_meta <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", 
                          header = TRUE, sep = "\t")
scarpa_meta <- pannen_meta[which(pannen_meta$Study == "Scarpa"),]
riemer_meta <- pannen_meta[which(pannen_meta$Study == "Riemer"),]
riemer_meta <- riemer_meta[which(riemer_meta$Subtype == "Cancer"),]
repset_meta <- rbind(scarpa_meta, riemer_meta)
repset_meta$Name <- as.character(repset_meta$Name)
repset_meta$Name[1:55] <- paste("X", repset_meta$Name[1:55], sep = "")
repset_meta <- repset_meta[match(colnames(repset), repset_meta$Name),]
all(repset_meta$Name == colnames(repset))
rownames(repset_meta) <- repset_meta$Name
all(rownames(repset_meta) == rownames(Repset_scdc_baron$decon_res$prop.est.mvw))

## form a dataframe out of the deco output --> maybe do this already in Execute_Deconvolution
decon_res <- as.data.frame(cbind(Repset_scdc_baron$decon_res$prop.est.mvw, t(Repset_scdc_baron$decon_res$yeval$spearmany.sample.table), 
                   t(Repset_scdc_baron$decon_res$yeval$mADy.sample.table), t(Repset_scdc_baron$decon_res$yeval$RMSDy.sample.table),
                   repset_meta$Grading))
colnames(decon_res) <- c(colnames(Repset_scdc_baron$decon_res$prop.est.mvw), "spearman", "mad", "rmsd", "response") ## p-value
decon_res[,1:(ncol(Repset_scdc_baron$decon_res$prop.est.mvw)+3)] <- sapply(1:(ncol(Repset_scdc_baron$decon_res$prop.est.mvw)+3), 
                                                                           function(x) as.numeric(decon_res[,x]))
str(decon_res)


## descriptive statistics (necessary? kinda)
skimmed <- skim(decon_res)

## impute missing values if there are any
if (any(skimmed$n_missing > 0)){
  preProcess_missingdata_model <- preProcess(decon_res, method = "knnImpute")
  decon_res <- predict(preProcess_missingdata_model, newdata = decon_res)
}

## preprocess data --> which method?
preProcess_range_model <- preProcess(decon_res, method = "scale")
decon_res <- predict(preProcess_range_model, newdata = decon_res)

## split into training and test set (make it optional?)
set.seed(4)
trainRowNumbers <- createDataPartition(decon_res$response, p = 0.8, list = FALSE)
trainData <- decon_res[trainRowNumbers,]
testData <- decon_res[-trainRowNumbers,]
x = trainData[, - ncol(trainData)]
y = trainData$response

## analyze feature importance
featurePlot(x = trainData[, - ncol(trainData)], 
            y = factor(trainData$response), 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

featurePlot(x = trainData[, - ncol(trainData)], 
            y = factor(trainData$response), 
            plot = "scatter",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

## train RF model and predict on training data itself
fitControl <- trainControl(
  method = "cv",
  number = 10,
  sampling = "down",
  savePred=T
)

model_rf <- train(response ~ ., data = trainData, 
                  method = 'rf', 
                  trControl = fitControl, 
                  type = "Classification",
                  metric = "Accuracy",
                  ntree = 500)
fitted <- predict(model_rf)
varimp_rf <- varImp(model_rf)
plot(varimp_rf)

## predict on test data
predicted <- predict(model_rf, testData)
confusionMatrix(reference = factor(testData$response), data = predicted, mode = "everything", positive = "G3")

