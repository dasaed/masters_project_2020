#!/usr/bin/env Rscript
# A) PRE-SETUP ----
# 1. Clear console and environment ----

clc <- function() cat("\014") ; clc()
remove(list = ls())

# 2. Install and load libraries ----

#uncomment the following lines if you need to install the libraries
# install.packages("h2o")
# install.packages("lime")
# install.packages("ROCR")
# install.packages("PRROC")
# install.packages("tidyverse")
library("lime")       # ML local interpretation
library("h2o")        # ML model building
library("ROCR")       # ML evaluation
library("PRROC")      # ML evaluation
library("tidyverse")

# 3. Load and preprocess data ----

# Training Data
fullTrainingDataSet <- read.csv("./tests/combinedData.csv", header = TRUE)
fullTrainingDataSet[,"Class"] <- as.logical(fullTrainingDataSet[,"Class"]) # guarantees that
                          # our h2o model trains the random forest (RF) as a classification 
                          # tree and not a regression tree

dataSetTrainX <- fullTrainingDataSet[,-(8:9)] # to be used by orig training RF model later on
dataSetTrainY <- fullTrainingDataSet[,(8:9)] 

# Testing Data
slt2dataPos <- read.csv("./testing_datasets/SLT2_Positives.tsv", sep = "\t", header = TRUE)
slt2dataPos$Class <- rep(1,nrow(slt2dataPos))
slt2dataNeg <- read.csv("./testing_datasets/SLT2_Negatives.tsv", sep = "\t", header = TRUE)
slt2dataNeg$Class <- rep(0,nrow(slt2dataNeg))
slt2data <- rbind(slt2dataPos,slt2dataNeg)

ludataPos <- read.csv("./testing_datasets/Lu_Positives.tsv", sep = "\t", header = TRUE)
ludataPos$Class <- rep(1,nrow(ludataPos))
ludataNeg <- read.csv("./testing_datasets/Lu_Negatives.tsv", sep = "\t", header = TRUE)
ludataNeg$Class <- rep(0,nrow(ludataNeg))
ludata <- rbind(ludataPos,ludataNeg)



# 4. Load original model ----

# The original article, titled "Prioritizing bona fide bacterial small RNAs with machine learning 
# classifiers", and the source materials(training data, testing data, R scripts, .rds file, etc.) 
# can all be found in PeerJ under the following url: https://peerj.com/articles/6304/
# The original .rds file can be found in github under the following link:
# https://github.com/BioinformaticsLabAtMUN/sRNARanking


## 4.1 Load model from rds file ----

origRF <- readRDS("RF_classifier4sRNA.rds")

## 4.2 Create functions for compatibility  ----

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
# 1. Intialize H2O and load training data ----

h2o.init(              # Initialize h2o
  nthreads = -1,       # -1 = use all available threads
  max_mem_size = "4G"  # specify the amount of memory to use
)
h2o.removeAll()        # clean up the system, h2o-wise (ie. kill any running h2o clusters)
h2o.no_progress() 
trainData <- as.h2o(fullTrainingDataSet) 

# 2. Build model ----

rfh2o <- h2o.randomForest( 
             training_frame = trainData,
             x = 1:7,                   # features to use to generate the prediction
             y = 9,                     # Class type -> what we want to predict
             model_id = "rf1_sRNA",     # name of model in h2o
             ntrees = 1000,              # max number of trees  
             seed = 1234,               # seed, has to be set WITHIN the h2o function
                                        # and it's supposed to be different from "R's seed", so 
                                        # results might not be exactly the same as orig model, but
                                        # should be similar enough
             mtries = 2,                 # Same as original model 
             max_depth = 30,
             verbose = TRUE
)

# 3. Preview Model's Performance ----

rfh2o@model$variable_importances
h2o.varimp_plot(rfh2o)
rfh2o_peformance <- h2o.performance(rfh2o)
rfh2o_peformance

# C) COMPARE NEW AND ORIG MODEL ----

# 1. Implement Orig RF evaluation function ----
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
  res$pr <- pr.curve(scores.class0 = res$predD[labels == 1,2], scores.class1 = res$predD[labels == 0,2], curve = T)
  return(res)
}

origRF_performance <- evaluateData(origRF, dataSetTrainX, dataSetTrainY[,2])
origRF_performance
origRF_performance$auc@y.values

# 2. Get Predictions from the Models ----

origRF_lu_pred <- predict(origRF, ludata[,-8], type = "prob")
origRF_slt2_pred <- predict(origRF, slt2data[,-8], type = "prob")

rfh2o_slt2_pred <- h2o.predict(object = rfh2o, newdata = as.h2o(slt2data[,-8]) )
rfh2o_lu_pred <- h2o.predict(object = rfh2o, newdata = as.h2o(ludata[,-8]) )

# 3. Build Comparison Tables ----

slt2_predictions <-  cbind(as.data.frame(origRF_slt2_pred),         # Orig RF predictions 
                           as.data.frame(rfh2o_slt2_pred[,c(2:3)]), # Orig RF predictions
                           slt2data[,8]                             # Actual Value
                           )
colnames(slt2_predictions) <- c("OrigRF_F","OrigRF_T","RF_H2O_F","RF_H2O_T","Real_Value")

lu_predictions <-  cbind(as.data.frame(origRF_lu_pred),
                         as.data.frame(rfh2o_lu_pred[,c(2:3)]),
                         ludata[,8]
                         )
colnames(lu_predictions) <- c("OrigRF_F","OrigRF_T","RF_H2O_F","RF_H2O_T","Real_Value")

compareAnswers <- function(a,b,c){
  if ( a == c & b == c ){
    res="BOTH_RIGHT"
    #res="BOTH"
  }
  else if ( a == c & b != c) {
    res="OnlyA"    
  }
  else if ( a != c & b == c) {
    res="OnlyB"    
  }
  else{
    res="BOTH_WRONG"
    #res="BOTH"  
  }
  return(res)
}

slt2_predictions$origPreds <- ifelse( slt2_predictions$OrigRF_T >= 0.5, 1, 0)
slt2_predictions$h2oPreds  <- ifelse( slt2_predictions$RF_H2O_T >= 0.5, 1, 0)


for( i in 1:nrow(slt2_predictions) ){
  slt2_predictions[i,]$Similarities <- compareAnswers(
    slt2_predictions[i,]$origPreds,
    slt2_predictions[i,]$h2oPreds,
    slt2_predictions[i,]$Real_Value
  )
}
head(slt2_predictions, 30)


#write.csv(slt2_predictions, "./slt2_predictions.csv", row.names = FALSE)
differences <-slt2_predictions[slt2_predictions$Comparison != "BOTH", ]

# D) LIME ----
# 1. Apply LIME to the new RF models ----
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



# E) PDP ----
# 1. Generate PDPs for the RF models ----

h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "SS")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "Pos10wrtsRNAStart")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "DistTerm")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "Distance")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "DownDistance")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "sameStrand")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "sameDownStrand")




# F) SHAPley??? ----
# 1. ----
