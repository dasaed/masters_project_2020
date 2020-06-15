#!/usr/bin/env Rscript
# A) PRE-SETUP ----
# ..1. Clear console and environment ----

clc <- function() cat("\014") ; clc()
remove(list = ls())

# ..2. Install and load libraries ----

#uncomment the following lines if you need to install the libraries
# install.packages("h2o")
# install.packages("lime")
# install.packages("ROCR")
# install.packages("PRROC")
# install.packages("randomForest")

library("lime")          # ML local interpretation
library("h2o")           # ML model building
library("ROCR")          # ML model performance measuring
library("PRROC")         # ML model performance measuring
library("randomForest")  # ML model building

# ..3. Load and preprocess data ----

# Training Data
fullTrainingDataSet <- read.csv("./tests/combinedData.csv", header = TRUE)
fullTrainingDataSet[,"Class"] <- as.logical(fullTrainingDataSet[,"Class"]) # guarantees that
                          # our h2o model trains the random forest (RF) as a classification 
                          # tree and not a regression tree

dataSetTrainX <- fullTrainingDataSet[,-(8:9)] # to be used by orig training RF model later on
dataSetTrainY <- fullTrainingDataSet[,(8:9)] 

# Testing Data
slt2dataPos <- read.csv("./testing_datasets/SLT2_Positives.tsv", sep = "\t", header = TRUE)
slt2dataNeg <- read.csv("./testing_datasets/SLT2_Negatives.tsv", sep = "\t", header = TRUE)
slt2data <- rbind(slt2dataPos,slt2dataNeg)

ludataPos <- read.csv("./testing_datasets/Lu_Positives.tsv", sep = "\t", header = TRUE)
ludataNeg <- read.csv("./testing_datasets/Lu_Negatives.tsv", sep = "\t", header = TRUE)
lu2data <- rbind(ludataPos,ludataNeg)



# ..4. Load original model ----

# The original article, titled "Prioritizing bona fide bacterial small RNAs with machine learning 
# classifiers", and the source materials(training data, testing data, R scripts, .rds file, etc.) 
# can all be found in PeerJ under the following url: https://peerj.com/articles/6304/
# The original .rds file can be found in github under the following link:
# https://github.com/BioinformaticsLabAtMUN/sRNARanking


## ....4.1 Load model from rds file ----

origRF <- readRDS("RF_classifier4sRNA.rds")

## ....4.2 Create functions for compatibility  ----

# The original model was built using the randomForest library found in CRAN. However, even though
# LIME is supposed to be model agnostic, it's current R implementation only supports certain
# models out of the box. To extend support of the LIME library to other models, the following
# functions need to be defined for each model type. 

predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, ...)
  switch(
    type,
    raw = data.frame(Response = ifelse(res[,2] > 0.5, "sRNA", "notSRNA"), 
                     stringsAsFactors = FALSE
    ),
    prob = res 
  )
  print(class(res))
  print(dim(res))
  print(res)
}

model_type.randomForest <- function(x, ...) 'classification'

# B) TRAIN THE NEW H2O MODEL ----
# ..1. Intialize H2O and load training data ----

h2o.init(              # Initialize h2o
  nthreads = -1,       # -1 = use all available threads
  max_mem_size = "4G"  # specify the amount of memory to use
)
h2o.removeAll()        # clean up the system, h2o-wise (ie. kill any running h2o clusters)
#h2o.no_progress() 
trainData <- as.h2o(fullTrainingDataSet) 

# ..2. Build model ----
?h2o.randomForest
rfh2o <- h2o.randomForest( 
             training_frame = trainData,
             x = 1:7,                   # features to use to generate the prediction
             y = 9,                     # Class type -> what we want to predict
             model_id = "rf1_sRNA",     # name of model in h2o
             ntrees = 400,              # max number of trees  
             seed = 1234,               # seed, has to be set WITHIN the h2o function
                                        # and it's supposed to be different from "R's seed", so 
                                        # results might not be exactly the same as orig model, but
                                        # should be similar enough
             mtries = 2,                 # Same as original model 
             max_depth = 30, 
             verbose = TRUE
)
rfh2o # This will give us a preview

# ..3. Preview model predictions ----
# ....3.1 With SLT2 data ----

testData_slt2 <- as.h2o(slt2data)
rfh2o_pred_slt2 <- h2o.predict(object = rfh2o, newdata = testData_slt2 )

dim(rfh2o_pred_slt2)
head(rfh2o_pred_slt2)
tail(rfh2o_pred_slt2)

#rfh2o_pred_slt2_performance <- h2o.performance(rfh2o, testData_slt2)

# ....3.2 With LU data ----

testData_lu <- as.h2o(lu2data)
rfh2o_pred_lu <- h2o.predict(object = rfh2o, newdata = testData_lu )

dim(rfh2o_pred_lu)
head(rfh2o_pred_lu)
tail(rfh2o_pred_lu)

#rfh2o_pred_lu_performance <- h2o.performance(rfh2o, testData_lu)

# C) COMPARE NEW AND ORIG MODEL ----

# ..1. Implement evaluation function ----
# This function was copied directly from the source materials, as it was the one used
# to evaluate the original model

evaluateData <- function(RF, data, labels){
  require(ROCR)
  require(PRROC)
  require(randomForest)
  res <- list()
  #obtain predictions
  res$predD <- predict(RF, data, type = "prob")
  #Evaluate predictions
  res$pred <- prediction(res$predD[,2], as.logical(labels))
  res$PR <- performance(res$pred, measure = "prec", x.measure = "rec")
  res$SS <- performance(res$pred, measure="sens", x.measure="spec")
  res$auc  <- performance(res$pred, measure = "auc")
  res$acc <-  performance(res$pred, measure = "acc")
  res$pr <- pr.curve(scores.class0 = res$predD[labels == 1,2], 
                     scores.class1 = res$predD[labels == 0,2], 
                     curve = T)
  return(res)
}

# ..2. Get evaluations from models ----

origRF_eval_with_trainData <- evaluateData(origRF, # Evaluations from orig RF
                                           dataSetTrainX, 
                                           dataSetTrainY[,2]) 
rfh2o_peformance <- h2o.performance(rfh2o) # Evaluations from h2o RF

### Variable Importance ###
varImp <- importance(origRF, scale = TRUE)
varImpPlot(origRF, sort = TRUE, type = 1, main = "Initial RF Variable Importance")
rfh2o@model$variable_importances
h2o.varimp_plot(rfh2o, main = "H2O RF Variable Importance")

### AUC Comparison ###

origRF_eval_with_trainData$auc@y.name
origRF_eval_with_trainData$auc@y.values # orig RF
h2o.auc(rfh2o)

### Accuracy comparison ###
plot(origRF_eval_with_trainData$acc, main = "Initial RF Accuracy Plot") # orig RF
plot(h2o.accuracy(rfh2o_peformance), main = "H2O RF Accuraccy Plot")

### Precision-Recall ###
plot(origRF_eval_with_trainData$pr, main = "Initial RF PR Plot")
head(origRF_eval_with_trainData$predD, 2)
class(origRF_eval_with_trainData$predD)

pr.curve(scores.class0 = origRF_eval_with_trainData$predD[labels == 1,2], 
         scores.class1 = origRF_eval_with_trainData$predD[labels == 0,2], 
         curve = T)

rfh2o_train_predictions <- h2o.predict(object = rfh2o, newdata = as.h2o(dataSetTrainX) )
rf2h2o_train_predictions_df <- as.data.frame(rfh2o_train_predictions)
pr.curve(scores.class0 = as.matrix(rf2h2o_train_predictions_df[labels == 1,3]), 
         scores.class1 = as.list.data.frame(rf2h2o_train_predictions_df[labels == 0,3]), 
         curve = T)


### Results on training data ###
rfh2o_td_preds <- h2o.predict(object = rfh2o, newdata = trainData[,c(1:7)])
head(training_data_predictions)
dim(training_data_predictions)

origRF_td_preds <- h2o.predict(object = rfh2o, newdata = trainData[,c(1:7)])
head(training_data_predictions)
dim(training_data_predictions)


# D) LIME ----

# ..1. Apply LIME to the new RF models ----
lime_explainer_rfh2o <- lime( as.data.frame(trainData[,c(1:7)]), # original training data
                             rfh2o,
                             bin_continuous = TRUE,
                             quantile_bins = FALSE
                            )
lime_explanations_rfh2o <- explain( as.data.frame(testData_slt2[1,]),   # Data to explain
                              lime_explainer_rf1,     # Explainer to use
                              n_labels = 1, # only 1 type of category
                              n_features = 7, # Number of features we want to use for explanation
                              n_permutations = 250,
                              feature_select = "none"
)
lime_explanations_rfh2o
plot_features(lime_explanations)

# E) Generate PDPs for the RF models ----

h2o.partialPlot(rfh2o, data = trainData, cols = "SS")
h2o.partialPlot(rfh2o, data = trainData, cols = "Pos10wrtsRNAStart")
h2o.partialPlot(rfh2o, data = trainData, cols = "DistTerm")
h2o.partialPlot(rfh2o, data = trainData, cols = "Distance")
h2o.partialPlot(rfh2o, data = trainData, cols = "DownDistance")
h2o.partialPlot(rfh2o, data = trainData, cols = "sameStrand")
h2o.partialPlot(rfh2o, data = trainData, cols = "sameDownStrand")

h2o.partialPlot(rfh2o, data = trainData, cols = c("Distance","DownDistance"))

# F) Generate Shapley Values for the H2O RF ----

?predict_contributions.H2OModel()
predict_contributions.H2OModel(rfh2o, testData_slt2[1,])

# G) Microsoft R API

# install the latest version from CRAN
install.packages("azuremlsdk")
azuremlsdk::install_azureml(envname = 'r-reticulate')
library(azuremlsdk)

