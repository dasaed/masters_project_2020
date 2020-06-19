#!/usr/bin/env Rscript
# Some terminology
# OrigRF/Orig RF = Original random forest as found in https://github.com/BioinformaticsLabAtMUN/sRNARanking

# A) PRE-SETUP ----
# 1. Clear console and environment ----

clc <- function() cat("\014") ; clc()
remove(list = ls())

# 2. Install and load libraries ----

# uncomment the following lines if you need to install the libraries
# install.packages("lime")
# install.packages("ROCR")
# install.packages("PRROC")
# install.packages("tidyverse")
# install.packages("ggplot2") 
library("lime")         # ML local interpretation
library("ROCR")         # ML evaluation
library("PRROC")        # ML evaluation
library("randomForest") # ML model building
#library("tidyverse")    # Graphing purposes
#library("ggplot2")      # Graphing purposes

## The following is copied from the h2o documentation site, as it's their recommended way for
## downloading and installing h2o for R (http://h2o-release.s3.amazonaws.com/h2o/rel-zahradnik/4/index.html)
## The following two commands remove any previously installed H2O packages for R.
#if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
#if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
## Next, we download packages that H2O depends on.
#pkgs <- c("RCurl","jsonlite")
#for (pkg in pkgs) {
#  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
#}
## Now we download, install and initialize the H2O package for R.
#install.packages("h2o", type="source", repos="http://h2o-release.s3.amazonaws.com/h2o/rel-zahradnik/4/R")
library("h2o")          # ML model building

# 3. Load and preprocess data ----
### Feature Explanations
## "SS" = the free energy of the predicted secondary structure of the sRNA - mostly negative values
## "Pos10wrtsRNAStart" = distance to their closest predicted promoter site
## "DistTerm" = distance to their closest predicted Rho-independent terminator 
## "Distance", = distance to the closest reading frame on the LEFT("upstream") side
## "sameStrand" = boolean, if transcription is going in the same direction as the ORF(Left Open Reading Frame)
#     ORF = genomic sequence that's supposed to code for a protein
## "DownDistance"  = distance to the closest reading frame on the RIGHT("downstream") side
## "sameDownStrand" =  boolean, if transcription is going in the same direction as the RORF(Right Open Reading Frame)
# 
# Load Training Data
fullTrainingDataSet <- read.csv("./tests/combinedData.csv", header = TRUE)
fullTrainingDataSet[,"Class"] <- as.logical(fullTrainingDataSet[,"Class"]) # guarantees that
                          # our h2o model trains the random forest (RF) as a classification 
                          # tree and not a regression tree

#dataSetTrainX <- fullTrainingDataSet[,-(8:9)] # UNCOMMENT THIS TO RETRAIN ORIG RF MODEL
#dataSetTrainY <- fullTrainingDataSet[,(8:9)]  # UNCOMMENT THIS TO RETRAIN ORIG RF MODEL

# Scale Training Data
head(fullTrainingDataSet)
dataSetTrain_scaled <- scale(fullTrainingDataSet[,-c(5,7,8,9)], center = TRUE, scale = TRUE)
dataSetTrain_scaled <- cbind( dataSetTrain_scaled[,(1:4)], 
                              fullTrainingDataSet[,5],
                              dataSetTrain_scaled[,5],
                              fullTrainingDataSet[,c(5,7,8,9)])
head(dataSetTrainX_scaled)
colnames(dataSetTrain_scaled) <- colnames(fullTrainingDataSet)


# Normal distribution? Useful for the LIME settings later in section D
hist(x = as.numeric(dataSetTrainX[,1]))
qqnorm(as.numeric(dataSetTrainX[,1]))
qqline(as.numeric(dataSetTrainX[,1]))

hist(x = as.numeric(dataSetTrainX[,2]))
qqnorm(as.numeric(dataSetTrainX[,2]))
qqline(as.numeric(dataSetTrainX[,2]))


hist(x = as.numeric(dataSetTrainX[,3]))
qqnorm(as.numeric(dataSetTrainX[,3]))
qqline(as.numeric(dataSetTrainX[,3]))

hist(x = as.numeric(dataSetTrainX[,4]))
qqnorm(as.numeric(dataSetTrainX[,4]))
qqline(as.numeric(dataSetTrainX[,4]))

hist(x = as.numeric(dataSetTrainX[,5]))
qqnorm(as.numeric(dataSetTrainX[,5]))
qqline(as.numeric(dataSetTrainX[,5]))

hist(x = as.numeric(dataSetTrainX[,6]))
qqnorm(as.numeric(dataSetTrainX[,6]))
qqline(as.numeric(dataSetTrainX[,6]))

hist(x = as.numeric(dataSetTrainX[,7]))
qqnorm(as.numeric(dataSetTrainX[,7]))
qqline(as.numeric(dataSetTrainX[,7]))


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


## ..4.1 Load model from rds file ----

origRF <- readRDS("RF_classifier4sRNA.rds")

## ..4.2 Create functions for compatibility  ----

# The original model was built using the randomForest library found in CRAN. However, even though
# LIME is supposed to be model agnostic, it's current R implementation only supports certain
# models out of the box. To extend support of the LIME library to other models, the following
# functions need to be defined for each model type. 

predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, type = "prob")
  switch(
    type,
    raw = data.frame(Response = ifelse(res[,2] > 0.5, "sRNA", "notSRNA"), 
                     stringsAsFactors = FALSE
    ),
    prob = as.data.frame(res, check.names = F)
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
             ntrees = 400,              # max number of trees  
             seed = 1234,               # seed, has to be set WITHIN the h2o function
                                        # and it's supposed to be different from "R's seed", so 
                                        # results might not be exactly the same as orig model, but
                                        # should be similar enough
             mtries = 2,                 # Same as original model 
             max_depth = 30
)

# 3. Preview H2O's RF ----
# This is the performance with the OOB error, based on the training process (ie training data)
rfh2o
rfh2o@model$variable_importances
h2o.varimp_plot(rfh2o)
rfh2o_training_peformance <- h2o.performance(rfh2o)
rfh2o_training_peformance

# C) COMPARE MODEL'S PREDICTIONS ----
# In this section we'll get the predictions both models generate on the testing data LU and SLT2. The goal
# of this section is not to prove which model is the best, but to prove that the models are similar enough
# to the point that the analysis of the H2O model will be a valid analogy model for the original RF model. 

# 1. Get Testing Data Predictions from the Models ----

origRF_lu_pred <- predict(origRF, ludata[,-8], type = "prob")
origRF_slt2_pred <- predict(origRF, slt2data[,-8], type = "prob")

rfh2o_slt2_pred <- h2o.predict(object = rfh2o, newdata = as.h2o(slt2data[,-8]) )
rfh2o_lu_pred <- h2o.predict(object = rfh2o, newdata = as.h2o(ludata[,-8]) )

# 2. Create a comparison function ----
# This function is to simplify the comparison between the predictions from the Original RF model and the H2O RF model
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

# 3. SLT2 comparisons ----
slt2_predictions <-  cbind(as.data.frame(slt2data),               # Input
                           as.data.frame(origRF_slt2_pred[,2]),   # Orig RF predictions 
                           as.data.frame(rfh2o_slt2_pred[,3])     # Orig RF predictions
                           )

colnames(slt2_predictions) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", 
                                "sameStrand", "DownDistance", "sameDownStrand", "Class",
                                "OrigRF_T","RF_H2O_T")
slt2_predictions
cor(slt2_predictions[,"OrigRF_T"],slt2_predictions[,"RF_H2O_T"]) # the closest to 1, the more correlated(similar) they are
slt2_predictions$origPreds <- ifelse( slt2_predictions$OrigRF_T >= 0.5, 1, 0) 
slt2_predictions$h2oPreds  <- ifelse( slt2_predictions$RF_H2O_T >= 0.5, 1, 0) 
slt2_predictions$Similarities <- NA

for( i in 1:nrow(slt2_predictions) ){
  # Based on the way we are feeding the values to the compareAnswers function, 
  # Similiarities = OnlyA means that only the Orig RF got the right answer
  # Similiarities = OnlyB means that only the H2O RF got the right answer
  # All other answers will tell us where both models were right("BOTH_RIGHT"), or
  # wrong ("BOTH_WRONG")
  slt2_predictions[i,]$Similarities <- compareAnswers(
    slt2_predictions[i,]$origPreds,  
    slt2_predictions[i,]$h2oPreds,  
    slt2_predictions[i,]$Class 
  )
} 

head(slt2_predictions,2)

slt2_OnlyA <- slt2_predictions[slt2_predictions$Similarities == "OnlyA",] 
slt2_OnlyB <- slt2_predictions[slt2_predictions$Similarities == "OnlyB",]
slt2_BothW <- slt2_predictions[slt2_predictions$Similarities == "BOTH_WRONG",]

# By looking at the difference between the number of predictions the Orig. RF and the H2O RF,
# we can conclude that the models are similar enough in accuracy
nrow(slt2_OnlyA)
nrow(slt2_OnlyB)
nrow(slt2_BothW)

# 4. LU Comparisons ----

lu_predictions <-  cbind(as.data.frame(ludata),               # Input
                           as.data.frame(origRF_lu_pred[,2]),   # Orig RF predictions 
                           as.data.frame(rfh2o_lu_pred[,3])     # Orig RF predictions
)

colnames(lu_predictions) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", 
                                "sameStrand", "DownDistance", "sameDownStrand", "Class",
                                "OrigRF_T","RF_H2O_T")

lu_predictions$origPreds <- ifelse( lu_predictions$OrigRF_T >= 0.5, 1, 0) 
lu_predictions$h2oPreds  <- ifelse( lu_predictions$RF_H2O_T >= 0.5, 1, 0) 
lu_predictions$Similarities <- NA

for( i in 1:nrow(lu_predictions) ){
  # Based on the way we are feeding the values to the compareAnswers function, 
  # Similiarities = OnlyA means that only the Orig RF got the right answer
  # Similiarities = OnlyB means that only the H2O RF got the right answer
  # All other answers will tell us where both models were right("BOTH_RIGHT"), or
  # wrong ("BOTH_WRONG")
  lu_predictions[i,]$Similarities <- compareAnswers(
    lu_predictions[i,]$origPreds,  
    lu_predictions[i,]$h2oPreds,  
    lu_predictions[i,]$Class 
  )
} 

head(lu_predictions,2)

lu_OnlyA <- lu_predictions[lu_predictions$Similarities == "OnlyA",] 
lu_OnlyB <- lu_predictions[lu_predictions$Similarities == "OnlyB",]
lu_BothW <- lu_predictions[lu_predictions$Similarities == "BOTH_WRONG",]

# By looking at the difference between the number of predictions the Orig. RF and the H2O RF,
# we can conclude that the models are similar enough in accuracy
nrow(lu_OnlyA)
nrow(lu_OnlyB)
nrow(lu_BothW)
cor(lu_predictions[,"OrigRF_T"],lu_predictions[,"RF_H2O_T"]) # the closest to 1, the more correlated(similar) they are


# D) COMPARE MODEL'S METRICS ----
# In this section we'll get the model's metrics, and do a side by side comparison of
# the metrics. 
# 1. Get Performance Metrics ----
# This function was copied directly from the source materials, as it was the one used
# to evaluate the original model, and only the comments have been changed
evaluateData <- function(RF, data, labels){
  require(ROCR)
  require(PRROC)
  require(randomForest)
  res <- list()
  res$predD <- predict(RF, data, type = "prob")                        # Predictions as Generated by the model
  res$pred <- prediction(res$predD[,2], as.logical(labels))            # Predictions, but converted to S4 Object          
  res$PR <- performance(res$pred, measure = "prec", x.measure = "rec") # Precision Recall Curve
  res$SS <- performance(res$pred, measure="sens", x.measure="spec")    # Sensitivity vs. Specificity
  res$auc  <- performance(res$pred, measure = "auc")                   # Sensitivity vs. Specificity AUC
  res$acc <-  performance(res$pred, measure = "acc")                   # Accuracy
  res$pr <- pr.curve(scores.class0 = res$predD[labels == 1,2],         # Precision vs. Recall 
                     scores.class1 = res$predD[labels == 0,2], 
                     curve = T)
  return(res)
}

# OrigRF Performance Metrics
origRF_slt2_performance <- evaluateData(origRF, slt2data[,-8], slt2data[,8])
origRF_lu_performance <- evaluateData(origRF, ludata[,-8], ludata[,8])

# H2O RF Performance Metrics
# H2O provides a way to retrieve most of the desired metrics
# http://docs.h2o.ai/h2o/latest-stable/h2o-r/docs/reference/h2o.metric.html

slt2data_h2o <- slt2data
slt2data_h2o[,"Class"] <- as.logical(slt2data_h2o[,"Class"]) 
rfh2o_slt2_performance <- h2o.performance(rfh2o, newdata = as.h2o(slt2data_h2o))

ludata_h2o <- ludata
ludata_h2o[,"Class"] <- as.logical(ludata_h2o[,"Class"]) 
rfh2o_lu_performance <- h2o.performance(rfh2o, newdata = as.h2o(ludata_h2o))

# 2. Compare Metrics  ----
# ..2.1 Accuracy ----

metrics_table <- data.frame("Accuracy" = 1:4)
rownames(metrics_table) <- c("origRF_slt2_perf", "rfh2o_slt2_perf","origRF_lu_perf", "rfh2o_lu_perf")
metrics_table["origRF_slt2_perf","Accuracy"] <- nrow(slt2_predictions[slt2_predictions$Similarities == "BOTH_RIGHT" | 
                                                     slt2_predictions$Similarities == "OnlyA",]
                                    )/nrow(slt2_predictions)
metrics_table["rfh2o_slt2_perf","Accuracy"] <- nrow(slt2_predictions[slt2_predictions$Similarities == "BOTH_RIGHT" | 
                                                       slt2_predictions$Similarities == "OnlyB",]
                                    )/nrow(slt2_predictions)
metrics_table["origRF_lu_perf","Accuracy"] <- nrow(lu_predictions[lu_predictions$Similarities == "BOTH_RIGHT" | 
                                                     lu_predictions$Similarities == "OnlyA",]
                                    )/nrow(lu_predictions)
metrics_table["rfh2o_lu_perf","Accuracy"] <- nrow(lu_predictions[lu_predictions$Similarities == "BOTH_RIGHT" | 
                                                     lu_predictions$Similarities == "OnlyB",]
                                    )/nrow(lu_predictions)

metrics_table

## Accuracy graph based on different thresholds("cutoff")

plot(origRF_slt2_performance$acc, col = "blue") 
lines(h2o.accuracy(rfh2o_slt2_performance), type = "l", col = "red")

plot(origRF_lu_performance$acc, col = "blue")
lines(h2o.accuracy(rfh2o_lu_performance), type = "l", col = "red")

# ..2.2 AUCPR ----

metrics_table["origRF_slt2_perf","AUCPR"] <- origRF_slt2_performance$pr$auc.integral
metrics_table["rfh2o_slt2_perf", "AUCPR"] <- h2o.aucpr(rfh2o_slt2_performance)
metrics_table["origRF_lu_perf",  "AUCPR"] <- origRF_lu_performance$pr$auc.integral
metrics_table["rfh2o_lu_perf",  "AUCPR"] <- h2o.aucpr(rfh2o_lu_performance)

metrics_table

plot(origRF_slt2_performance$PR,  col = "blue", lwd = 2 )
lines(x = h2o.recall(rfh2o_slt2_performance)[,"tpr"],
      y = h2o.precision(rfh2o_slt2_performance)[,"precision"],
      col = "red", type = "l", lwd = 2
)

plot(origRF_lu_performance$PR,  col = "blue", lwd = 2 )
lines(x = h2o.recall(rfh2o_lu_performance)[,"tpr"],
      y = h2o.precision(rfh2o_lu_performance)[,"precision"],
      col = "red", type = "l", lwd = 2
)


# ..2.3 Sensitivy vs Specificity Graph ----

plot(origRF_slt2_performance$SS, col = "blue", lwd = 2)
#rfh2o_slt2_specificity <- h2o.specificity(rfh2o_slt2_performance)[,"tnr"] # 400 samples
#rfh2o_slt2_sensitivity <- h2o.sensitivity(rfh2o_slt2_performance)[,"tpr"] # 400 samples
lines(x = h2o.specificity(rfh2o_slt2_performance)[,"tnr"],
     y = h2o.sensitivity(rfh2o_slt2_performance)[,"tpr"],
     col = "red", type = "l", lwd = 2
     )

plot(origRF_lu_performance$SS, col = "blue", lwd = 2)
lines(x = h2o.specificity(rfh2o_lu_performance)[,"tnr"],
      y = h2o.sensitivity(rfh2o_lu_performance)[,"tpr"],
      col = "red", type = "l", lwd = 2
)


# What's the difference here? Is it just that "pr" is prettier and gives the AUCPR as well? 
origRF_slt2_performance$PR    
plot(origRF_slt2_performance$PR)
origRF_slt2_performance$pr    
plot(origRF_slt2_performance$pr)

# ..2.4 Other H2O Metrics ----
h2o.confusionMatrix(rfh2o_slt2_performance)
plot(h2o.F1(rfh2o_slt2_performance))
plot(rfh2o_slt2_performance, # REMINDER: TPR = Sensitivity, FPR = (1 - specificity)
     type = "roc", 
     col = "red",
     cex = 0.2,
     pch = 10
)

h2o.confusionMatrix(rfh2o_lu_performance)
plot(h2o.F1(rfh2o_lu_performance))
plot(rfh2o_lu_performance,
     type = "roc", 
     col = "red",
     cex = 0.2,
     pch = 10
)


# D) LOCAL METHODS ----
# Obtain Data to Analyze ----

head(slt2_OnlyA)
head(lu_OnlyA)

nrow(slt2_BothW)
nrow(lu_BothW)

# E) LIME ----

# LIME ARGUMENTS
# x	              The training data used for training the model that should be explained.
# model	          The model whose output should be explained
# preprocess	    Function to transform a character vector to the format expected from the model.
# bin_continuous  Should continuous variables be binned when making the explanation
# n_bins	        The number of bins for continuous variables if bin_continuous = TRUE
# quantile_bins	  Should the bins be based on n_bins quantiles or spread evenly over the range of
#                   the training data
# use_density     If bin_continuous = FALSE should continuous data be sampled using a kernel density 
#                 estimation. If not, continuous features are expected to follow a normal distribution.

# 1. Apply LIME to the new RF models ----

lime_explainer_rfh2o <- lime( as.data.frame( trainData[,c(1:7)] ), # original training data
                              rfh2o,
                              bin_continuous = FALSE, # Having this as T generates inconsistent explanations
                                                      # b/c of the mixture of numerical and categorical (T/F) 
                                                      # features
                              quantile_bins = FALSE,
                              use_density = TRUE,
                              n_bins = 10
)
nrow(slt2_OnlyA)
lime_explanations_rfh2o <- explain( as.data.frame(slt2_OnlyA[14,c(1:7)] ),   # Data to explain
                                    lime_explainer_rfh2o,     # Explainer to use
                                    n_labels = 1, # only 1 type of category
                                    n_features = 7, # Number of features we want to use for explanation
                                    n_permutations = 200,
                                    dist_fun = "gower"  # b/c training contains numerical and categorical(T/F) features
#                                    kernel_width = .75,
#                                    feature_select = "highest_weights"
)

plot_features(lime_explanations_rfh2o)
head(lime_explanations_rfh2o)
slt2_OnlyA[1,c(1:8)]
lime_explainer_rfh2o$n_bins


# F) SHAP Values ----
# 1. something ... ----
SHAP_H2O1 <- h2o.predict_contributions(rfh2o, as.h2o(slt2_OnlyA[14,c(1:7)]))
SHAP_H2O2 <- h2o.predict_contributions(rfh2o, as.h2o(slt2_OnlyA[14,c(1:7)]))
SHAP_H2O3 <- h2o.predict_contributions(rfh2o, as.h2o(slt2_OnlyA[14,c(1:7)]))
SHAP_H2O4 <- h2o.predict_contributions(rfh2o, as.h2o(slt2_OnlyA[14,c(1:7)]))
SHAP_H2O5 <- h2o.predict_contributions(rfh2o, as.h2o(slt2_OnlyA[14,c(1:7)]))

class(SHAP_H2O1)



slt2data[2,]


# G) PDP ----
# 1. Generate PDPs for the RF models ----

?partialPlot()
h2o.partialPlot(rfh2o, data = as.h2o(slt2data_h2o), cols = "SS", plot=TRUE, nbins=10)

h2o.partialPlot(rfh2o, data = as.h2o(ludata_h2o), cols = "SS")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "SS")

h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "Pos10wrtsRNAStart")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "DistTerm")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "Distance")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "DownDistance")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "sameStrand")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "sameDownStrand")

#TODO - Create model with scaled data and try LIME
#TODO - PDPs for classification tasks, instead of regression
#TODO - See how LIME and SHAP agree and compare

#### QUESTIONS ----
