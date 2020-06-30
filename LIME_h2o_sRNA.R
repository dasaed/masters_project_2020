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

library("lime")         # ML local interpretation
library("ROCR")         # ML evaluation
library("PRROC")        # ML evaluation
library("randomForest") # ML model building


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
 
# ..3.1 Load Training Data ----
trainDataSet <- read.csv("./tests/combinedData.csv", header = TRUE)
trainDataSet[,"Class"] <- as.logical(trainDataSet[,"Class"]) # guarantees that
                          # our h2o model trains the random forest (RF) as a classification
                          # tree and not a regression tree

# ..3.2 Scale and Center the Training Data ----
# the scale function can scale and center the data for us, but we would still need 
# the training data's standard deviation and means to be able to apply the same 
# normalization to our testing data sets. 

dataSetTrain_scaled <- scale(trainDataSet[,-c(5,7,8,9)], center = TRUE, scale = TRUE)
dataSetTrain_scaled <- cbind( dataSetTrain_scaled[,(1:4)], 
                              trainDataSet[,5],
                              dataSetTrain_scaled[,5],
                              trainDataSet[,c(7,8,9)])
colnames(dataSetTrain_scaled) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", 
                                   "sameStrand", "DownDistance", "sameDownStrand", "ID", "Class")
dataSetTrain_scaled[,"Class"] <- as.logical(dataSetTrain_scaled[,"Class"]) # Just making sure Class is taken a logical value
head(dataSetTrain_scaled,5)
# ..3.3 Normalize the Training Data ----
# Max and mins for regular normalization
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

dataSetTrain_norm <- trainDataSet
dataSetTrain_norm$SS <- normalize(dataSetTrain_norm$SS)
dataSetTrain_norm$Pos10wrtsRNAStart <- normalize(dataSetTrain_norm$Pos10wrtsRNAStart)
dataSetTrain_norm$DistTerm <- normalize(dataSetTrain_norm$DistTerm)
dataSetTrain_norm$Distance <- normalize(dataSetTrain_norm$Distance)
dataSetTrain_norm$DownDistance <- normalize(dataSetTrain_norm$DownDistance)
dataSetTrain_norm[,"Class"] <- as.logical(dataSetTrain_norm[,"Class"]) # Just making sure Class is taken a logical value

# ..3.4 Check for normality in the data ----
# There are certain parameters in the LIME functions that are recommended based on
# the distribution of the data. For best results, it's best to check out the data,
# and make the changes accordingly

dataSetTrainX <- trainDataSet[,-(8:9)] 
dataSetTrainY <- trainDataSet[,(8:9)]  

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


# ..3.5 Load testing Datasets ----
## Load Original Testing Data Sets
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

## Scale Testing Data Sets
### Xi    = An instance 
### Xmean = Mean/Average of all samples =  
### Xsd   = Standard Deviation of our data set 
### Zi    = Scaled Value of Xi
### Zi    = (Xi - Xmean) / Xsd 

# Standard Deviation and Mean for Standarization
scale_td <- function( Xi, Xmean, Xsd) {
  return ( (Xi - Xmean) / Xsd ) # Zi
}
trainData_sd_SS        <- sd(trainDataSet[,"SS"])
trainData_sd_Pos       <- sd(trainDataSet[,"Pos10wrtsRNAStart"])
trainData_sd_DisTerm   <- sd(trainDataSet[,"DistTerm"])
trainData_sd_Dis       <- sd(trainDataSet[,"Distance"])
trainData_sd_DownDis   <- sd(trainDataSet[,"DownDistance"])

trainData_mean_SS      <- mean(trainDataSet[,"SS"])
trainData_mean_Pos     <- mean(trainDataSet[,"Pos10wrtsRNAStart"])
trainData_mean_DisTerm <- mean(trainDataSet[,"DistTerm"])
trainData_mean_Dis     <- mean(trainDataSet[,"Distance"])
trainData_mean_DownDis <- mean(trainDataSet[,"DownDistance"])


slt2data_scaled <- slt2data
slt2data_scaled$SS                <- scale_td(slt2data_scaled$SS, 
                                              trainData_mean_SS, trainData_sd_SS)
slt2data_scaled$Pos10wrtsRNAStart <- scale_td(slt2data_scaled$Pos10wrtsRNAStart,
                                              trainData_mean_Pos, trainData_sd_Pos)
slt2data_scaled$DistTerm          <- scale_td(slt2data_scaled$DistTerm, 
                                              trainData_mean_DisTerm, trainData_sd_DisTerm)
slt2data_scaled$Distance          <- scale_td(slt2data_scaled$Distance,
                                              trainData_mean_Dis, trainData_sd_Dis)
slt2data_scaled$DownDistance      <- scale_td(slt2data_scaled$DownDistance,
                                              trainData_mean_DownDis, trainData_sd_DownDis)
ludata_scaled <- ludata
ludata_scaled$SS                <- scale_td(ludata_scaled$SS, 
                                              trainData_mean_SS, trainData_sd_SS)
ludata_scaled$Pos10wrtsRNAStart <- scale_td(ludata_scaled$Pos10wrtsRNAStart,
                                              trainData_mean_Pos, trainData_sd_Pos)
ludata_scaled$DistTerm          <- scale_td(ludata_scaled$DistTerm, 
                                              trainData_mean_DisTerm, trainData_sd_DisTerm)
ludata_scaled$Distance          <- scale_td(ludata_scaled$Distance,
                                              trainData_mean_Dis, trainData_sd_Dis)
ludata_scaled$DownDistance      <- scale_td(ludata_scaled$DownDistance,
                                              trainData_mean_DownDis, trainData_sd_DownDis)

## Normalized Testing Data sets
normalize_td <- function(x, min, max) {
  return ( (x - min) / (max - min) )
}

trainData_max_SS      <- max(trainDataSet[,"SS"])
trainData_max_Pos     <- max(trainDataSet[,"Pos10wrtsRNAStart"])
trainData_max_DisTerm <- max(trainDataSet[,"DistTerm"])
trainData_max_Dis     <- max(trainDataSet[,"Distance"])
trainData_max_DownDis <- max(trainDataSet[,"DownDistance"])

trainData_min_SS      <- min(trainDataSet[,"SS"])
trainData_min_Pos     <- min(trainDataSet[,"Pos10wrtsRNAStart"])
trainData_min_DisTerm <- min(trainDataSet[,"DistTerm"])
trainData_min_Dis     <- min(trainDataSet[,"Distance"])
trainData_min_DownDis <- min(trainDataSet[,"DownDistance"])

slt2data_norm <- slt2data
slt2data_norm$SS                <- normalize_td(slt2data_norm$SS, 
                                                 trainData_min_SS, trainData_max_SS)
slt2data_norm$Pos10wrtsRNAStart <- normalize_td(slt2data_norm$Pos10wrtsRNAStart,
                                                 trainData_min_Pos, trainData_max_Pos)
slt2data_norm$DistTerm          <- normalize_td(slt2data_norm$DistTerm,
                                        trainData_min_DisTerm, trainData_max_DisTerm)
slt2data_norm$Distance          <- normalize_td(slt2data_norm$Distance,
                                        trainData_min_Dis, trainData_max_Dis)
slt2data_norm$DownDistance      <- normalize_td(slt2data_norm$DownDistance,
                                            trainData_min_DownDis, trainData_max_DownDis)

ludata_norm <- ludata
ludata_norm$SS                <- normalize_td(ludata_norm$SS, 
                                                 trainData_min_SS, trainData_max_SS)
ludata_norm$Pos10wrtsRNAStart <- normalize_td(ludata_norm$Pos10wrtsRNAStart,
                                                 trainData_min_Pos, trainData_max_Pos)
ludata_norm$DistTerm          <- normalize_td(ludata_norm$DistTerm,
                                                 trainData_min_DisTerm, trainData_max_DisTerm)
ludata_norm$Distance          <- normalize_td(ludata_norm$Distance,
                                                 trainData_min_Dis, trainData_max_Dis)
ludata_norm$DownDistance      <- normalize_td(ludata_norm$DownDistance,
                                                 trainData_min_DownDis, trainData_max_DownDis)


# 4. Load original model ----

# The original article, titled "Prioritizing bona fide bacterial small RNAs with machine learning 
# classifiers", and the source materials(training data, testing data, R scripts, .rds file, etc.) 
# can all be found in PeerJ under the following url: https://peerj.com/articles/6304/
# The original .rds file can be found in github under the following link:
# https://github.com/BioinformaticsLabAtMUN/sRNARanking

## ..4.1 Load/Retrain Orig Model ----
# Load original model (save time)
origRF <- readRDS("RF_classifier4sRNA.rds")

tuneRF(combData2[,-8], y = factor(combData2[,8]), ntreeTry = 400, mtryStart = 2)
set.seed(1234)
origRF_retrained <- randomForest(x = dataSetTrainX, y = factor(dataSetTrainY[,9]), mtry = 2, ntree = 400, 
                                 importance = TRUE)


## ..4.2 Create functions for compatibility  ----

# The original model was built using the randomForest library found in CRAN. However, even though LIME is supposed to be model agnostic, it's current R implementation only supports certain models out of the box. To extend support of the LIME library to other models, the following functions need to be defined for each model type. 

predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, type = "prob")
  switch(
    type,
    raw = data.frame(Response = ifelse(as.vector(res[,2]) > 0.5, "1", "0"), stringsAsFactors = FALSE),
    #raw = data.frame(Response = ifelse(res[,2] > 0.5, "sRNA", "notSRNA"), stringsAsFactors = FALSE),
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
  max_mem_size = "6G"  # specify the amount of memory to use
)
h2o.removeAll()        # clean up the system, h2o-wise (ie. kill any running h2o clusters)
h2o.no_progress() 
trainData <- as.h2o(trainDataSet) 
trainData_scaled <- as.h2o(dataSetTrain_scaled) 
trainData_norm <- as.h2o(dataSetTrain_norm) 

# 2. Build model ----
# Since we are going to use the h2o library to perform our machine learning interpretability, we need to train our model, or an equivalent one, using h2o tools

# ..2.1 Build model with Orig Training data ----

rfh2o <- h2o.randomForest( 
             training_frame = trainData, # training using the normal data
             x = 1:7,                    # features to use to generate the prediction
             y = 9,                      # Class type -> what we want to predict
             model_id = "rf1_sRNA",      # name of model in h2o
             ntrees = 400,               # max number of trees  
             seed = 1234,                # seed, has to be set WITHIN the h2o function and it's supposed to be different from "R's seed", so results might not be exactly the same as orig model, but should be similar enough
             mtries = 2,                 # Same as original model 
             max_depth = 30
)

# ..2.2 Build model with Scaled Training data ----
rfh2o_scaled <- h2o.randomForest( 
  training_frame = trainData_scaled, # training using the standarized data
  x = 1:7,                   # features to use to generate the prediction
  y = 9,                     # Class type -> what we want to predict
  model_id = "rf2_sRNA",     # name of model in h2o
  ntrees = 400,              # max number of trees  
  seed = 1234,               # seed, has to be set WITHIN the h2o function and it's supposed to be different from "R's seed", so results might not be exactly the same as orig model, but should be similar enough
  mtries = 2,                 # Same as original model 
  max_depth = 30
)

# ..2.3 Build model with Norm Training data ----
rfh2o_norm <- h2o.randomForest( 
  training_frame = trainData_norm, # Training using the normalized data
  x = 1:7,                   # features to use to generate the prediction
  y = 9,                     # Class type -> what we want to predict
  model_id = "rf3_sRNA",     # name of model in h2o
  ntrees = 400,              # max number of trees  
  seed = 1234,               # seed, has to be set WITHIN the h2o function
  # and it's supposed to be different from "R's seed", so 
  # results might not be exactly the same as orig model, but
  # should be similar enough
  mtries = 2,                 # Same as original model 
  max_depth = 30
)
# ..2.4 BONUS: Build GLM model ----
# As we were doing tests with LIME, we discovered that it yielded irregular responses, and as a result, the explanations could not be trusted. It seems to be a known issue that LIME is not perfect for every model, and that it struggles with data that is highly diverse. In this use case, our data contains mixed numerical and categorical features, significantly different scales and magnitudes in our numercial features, and an inbalanced data sets that favors predictions of sRNAs being false. Since lime creates a local interpretable linear model for its explanations, the GLM constructed here was simply to verify a hypothesis that a linear model would be unable to explain the data. The "unbalanced-ness" seems to create models where the features tend always support a result of "false", while claiming that all the features are against a result of true. 

glmh2o <- h2o.glm( # https://docs.h2o.ai/h2o/latest-stable/h2o-docs/data-science/glm.html 
  training_frame = trainData, # training using the normal data
  x = 1:7,                   # features to use to generate the prediction
  y = 9,                     # Class type -> what we want to predict
  model_id = "glm_sRNA",     # name of model in h2o
  seed = 1234,
  family = "binomial"
)

glmh2o_perf <- h2o.performance(glmh2o)
glmh2o_perf
h2o.coef(glmh2o)
h2o.coef_norm(glmh2o)
glmh2o@model$coefficients_table
h2o.r2(glmh2o) # R^2 value is really low, and having a Max Recall being really low, meaning that LIME is also predicting that everything should be false

# 3. Preview H2O's RF ----
# This is the performance with the OOB error, based on the training process (ie training data)
rfh2o
rfh2o@model$variable_importances
h2o.varimp_plot(rfh2o)
rfh2o_training_peformance <- h2o.performance(rfh2o)
rfh2o_training_peformance

# C) COMPARE MODEL'S PREDICTIONS ----
# In this section we'll get the predictions both models generate on the testing data LU and SLT2. The goal of this section is not to prove which model is the best, but to prove that the models are similar enough to the point that the analysis of the H2O model will be a valid analogy model for the original RF model. 

# 1. Get Testing Data Predictions from the Models ----

origRF_lu_pred <- predict(origRF, ludata[,-8], type = "prob")
origRF_slt2_pred <- predict(origRF, slt2data[,-8], type = "prob")

rfh2o_slt2_pred <- h2o.predict(object = rfh2o, newdata = as.h2o(slt2data[,-8]) )
rfh2o_lu_pred <- h2o.predict(object = rfh2o, newdata = as.h2o(ludata[,-8]) )

rfh2o_scaled_slt2_pred <- h2o.predict(object = rfh2o_scaled, newdata = as.h2o(slt2data_scaled[,-8]) )
rfh2o_scaled_lu_pred <- h2o.predict(object = rfh2o_scaled, newdata = as.h2o(ludata_scaled[,-8]) )

rfh2o_norm_slt2_pred <- h2o.predict(object = rfh2o_norm, newdata = as.h2o(slt2data_norm[,-8]) )
rfh2o_norm_lu_pred <- h2o.predict(object = rfh2o_norm, newdata = as.h2o(ludata_norm[,-8]) )

glm_slt2_pred <- h2o.predict(object = glmh2o, newdata = as.h2o(slt2data[,-8]) )
glm_lu_pred <- h2o.predict(object = glmh2o, newdata = as.h2o(ludata[,-8]) )

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
# ..3.1 Build comparison table ----
slt2_predictions <-  cbind(as.data.frame(slt2data),               # Input
                           as.data.frame(origRF_slt2_pred[,2]),   # Orig RF predictions 
                           as.data.frame(rfh2o_slt2_pred[,3]),        # H2O RF predictions
                           as.data.frame(rfh2o_scaled_slt2_pred[,3]), # Scaled RF predictions
                           as.data.frame(rfh2o_norm_slt2_pred[,3])    # Normalized RF predictions
                           )
head(slt2_predictions, 3)
colnames(slt2_predictions) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", 
                                "sameStrand", "DownDistance", "sameDownStrand", "Class",
                                "OrigRF_T","RF_H2O_T","RFH2O_scaledT","RFH2O_normT")
head(slt2_predictions, 3)

# Check correlations between the OrigRF and the H2O models
# the closest to 1, the more correlated(similar) they are
cor(slt2_predictions[,"OrigRF_T"],slt2_predictions[,"RF_H2O_T"]) 
cor(slt2_predictions[,"OrigRF_T"],slt2_predictions[,"RFH2O_scaledT"]) 
cor(slt2_predictions[,"OrigRF_T"],slt2_predictions[,"RFH2O_normT"])

slt2_predictions$origPreds <- ifelse( slt2_predictions$OrigRF_T >= 0.5, 1, 0) 
slt2_predictions$h2oPreds  <- ifelse( slt2_predictions$RF_H2O_T >= 0.5, 1, 0) 
slt2_predictions$h2oScaledPreds  <- ifelse( slt2_predictions$RFH2O_scaledT >= 0.5, 1, 0) 
slt2_predictions$h2oNormPreds  <- ifelse( slt2_predictions$RFH2O_normT >= 0.5, 1, 0)  

slt2_predictions$origVsH2O <- NA
slt2_predictions$origVsH2OScaled <- NA
slt2_predictions$origVsH2ONorm <- NA



for( i in 1:nrow(slt2_predictions) ){
  # Based on the way we are feeding the values to the compareAnswers function, 
  # Similiarities = OnlyA means that only the Orig RF got the right answer
  # Similiarities = OnlyB means that only the H2O RF got the right answer
  # All other answers will tell us where both models were right("BOTH_RIGHT"), or
  # wrong ("BOTH_WRONG")
  slt2_predictions[i,]$origVsH2O <- compareAnswers(
    slt2_predictions[i,]$origPreds,  
    slt2_predictions[i,]$h2oPreds,  
    slt2_predictions[i,]$Class 
  )
  slt2_predictions[i,]$origVsH2OScaled <- compareAnswers(
    slt2_predictions[i,]$origPreds,  
    slt2_predictions[i,]$h2oScaledPreds,  
    slt2_predictions[i,]$Class 
  )
  slt2_predictions[i,]$origVsH2ONorm <- compareAnswers(
    slt2_predictions[i,]$origPreds,  
    slt2_predictions[i,]$h2oNormPreds,  
    slt2_predictions[i,]$Class 
  )
} 

head(slt2_predictions)
summary(slt2_predictions)
str(slt2_predictions)

# ..3.2 Compare OrigRF vs H2O Model ----

# By looking at the difference between the number of predictions the Orig. RF and the H2O RF,
# we can conclude that the models are similar enough in accuracy
slt2_OnlyA <- slt2_predictions[slt2_predictions$origVsH2O == "OnlyA",] 
slt2_OnlyB <- slt2_predictions[slt2_predictions$origVsH2O == "OnlyB",]
slt2_BothW <- slt2_predictions[slt2_predictions$origVsH2O == "BOTH_WRONG",]
slt2_BothR <- slt2_predictions[slt2_predictions$origVsH2O == "BOTH_RIGHT",]

nrow(slt2_OnlyA)
nrow(slt2_OnlyB)
nrow(slt2_BothW)
nrow(slt2_BothR)
nrow(slt2_OnlyA) + nrow(slt2_OnlyB) + nrow(slt2_BothW) + nrow(slt2_BothR)
nrow(slt2_predictions)

# ..3.3 Compare OrigRF vs H2O Model Scaled ----

# By looking at the difference between the number of predictions the Orig. RF and the H2O Scaled RF,
# we can conclude that the models are similar enough in accuracy
slt2_OnlyA_scaled <- slt2_predictions[slt2_predictions$origVsH2OScaled == "OnlyA",] 
slt2_OnlyB_scaled <- slt2_predictions[slt2_predictions$origVsH2OScaled == "OnlyB",]
slt2_BothW_scaled <- slt2_predictions[slt2_predictions$origVsH2OScaled == "BOTH_WRONG",]
slt2_BothR_scaled <- slt2_predictions[slt2_predictions$origVsH2OScaled == "BOTH_RIGHT",]

nrow(slt2_OnlyA_scaled)
nrow(slt2_OnlyB_scaled)
nrow(slt2_BothW_scaled)
nrow(slt2_BothR_scaled)
nrow(slt2_OnlyA_scaled) + nrow(slt2_OnlyB_scaled) + nrow(slt2_BothW_scaled) + nrow(slt2_BothR_scaled)
nrow(slt2_predictions)

# ..3.4 Compare OrigRF vs H2O Model Normalized ----

# By looking at the difference between the number of predictions the Orig. RF and the H2O RF,
# we can conclude that the models are similar enough in accuracy
slt2_OnlyA_norm <- slt2_predictions[slt2_predictions$origVsH2ONorm == "OnlyA",] 
slt2_OnlyB_norm <- slt2_predictions[slt2_predictions$origVsH2ONorm == "OnlyB",]
slt2_BothW_norm <- slt2_predictions[slt2_predictions$origVsH2ONorm == "BOTH_WRONG",]
slt2_BothR_norm <- slt2_predictions[slt2_predictions$origVsH2ONorm == "BOTH_RIGHT",]

nrow(slt2_OnlyA_norm)
nrow(slt2_OnlyB_norm)
nrow(slt2_BothW_norm)
nrow(slt2_BothR_norm)
nrow(slt2_OnlyA_norm) + nrow(slt2_OnlyB_norm) + nrow(slt2_BothW_norm) + nrow(slt2_BothR_norm)
nrow(slt2_predictions)

# 4. LU Comparisons ----
# ..4.1 Build comparison table ----
lu_predictions <-  cbind(  as.data.frame(ludata),               # Input
                           as.data.frame(origRF_lu_pred[,2]),   # Orig RF predictions 
                           as.data.frame(rfh2o_lu_pred[,3]),        # H2O RF predictions
                           as.data.frame(rfh2o_scaled_lu_pred[,3]), # Scaled RF predictions
                           as.data.frame(rfh2o_norm_lu_pred[,3])    # Normalized RF predictions
)
colnames(lu_predictions) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", 
                                "sameStrand", "DownDistance", "sameDownStrand", "Class",
                                "OrigRF_T","RF_H2O_T","RFH2O_scaledT","RFH2O_normT")

# Check correlations between the OrigRF and the H2O models
# the closest to 1, the more correlated(similar) they are
cor(lu_predictions[,"OrigRF_T"],lu_predictions[,"RF_H2O_T"]) 
cor(lu_predictions[,"OrigRF_T"],lu_predictions[,"RFH2O_scaledT"]) 
cor(lu_predictions[,"OrigRF_T"],lu_predictions[,"RFH2O_normT"])

lu_predictions$origPreds <- ifelse( lu_predictions$OrigRF_T >= 0.5, 1, 0) 
lu_predictions$h2oPreds  <- ifelse( lu_predictions$RF_H2O_T >= 0.5, 1, 0) 
lu_predictions$h2oScaledPreds  <- ifelse( lu_predictions$RFH2O_scaledT >= 0.5, 1, 0) 
lu_predictions$h2oNormPreds  <- ifelse( lu_predictions$RFH2O_normT >= 0.5, 1, 0)  

lu_predictions$origVsH2O <- NA
lu_predictions$origVsH2OScaled <- NA
lu_predictions$origVsH2ONorm <- NA

for( i in 1:nrow(lu_predictions) ){
  # Based on the way we are feeding the values to the compareAnswers function, 
  # Similiarities = OnlyA means that only the Orig RF got the right answer
  # Similiarities = OnlyB means that only the H2O RF got the right answer
  # All other answers will tell us where both models were right("BOTH_RIGHT"), or
  # wrong ("BOTH_WRONG")
  lu_predictions[i,]$origVsH2O <- compareAnswers(
    lu_predictions[i,]$origPreds,  
    lu_predictions[i,]$h2oPreds,  
    lu_predictions[i,]$Class 
  )
  lu_predictions[i,]$origVsH2OScaled <- compareAnswers(
    lu_predictions[i,]$origPreds,  
    lu_predictions[i,]$h2oScaledPreds,  
    lu_predictions[i,]$Class 
  )
  lu_predictions[i,]$origVsH2ONorm <- compareAnswers(
    lu_predictions[i,]$origPreds,  
    lu_predictions[i,]$h2oNormPreds,  
    lu_predictions[i,]$Class 
  )
} 

head(lu_predictions)
summary(lu_predictions)
str(lu_predictions)

# ..4.2 Compare OrigRF vs H2O Model ----
# By looking at the difference between the number of predictions the Orig. RF and the H2O RF,
# we can conclude that the models are similar enough in accuracy
lu_OnlyA <- lu_predictions[lu_predictions$origVsH2O == "OnlyA",] 
lu_OnlyB <- lu_predictions[lu_predictions$origVsH2O == "OnlyB",]
lu_BothW <- lu_predictions[lu_predictions$origVsH2O == "BOTH_WRONG",]
lu_BothR <- lu_predictions[lu_predictions$origVsH2O == "BOTH_RIGHT",]

nrow(lu_OnlyA)
nrow(lu_OnlyB)
nrow(lu_BothW)
nrow(lu_BothR)
nrow(lu_OnlyA) + nrow(lu_OnlyB) + nrow(lu_BothW) + nrow(lu_BothR)
nrow(lu_predictions)

# ..4.3 Compare OrigRF vs H2O Model Scaled ----

lu_OnlyA_scaled <- lu_predictions[lu_predictions$origVsH2OScaled == "OnlyA",] 
lu_OnlyB_scaled <- lu_predictions[lu_predictions$origVsH2OScaled == "OnlyB",]
lu_BothW_scaled <- lu_predictions[lu_predictions$origVsH2OScaled == "BOTH_WRONG",]
lu_BothR_scaled <- lu_predictions[lu_predictions$origVsH2OScaled == "BOTH_RIGHT",]

nrow(lu_OnlyA_scaled)
nrow(lu_OnlyB_scaled) # Scaled actually performed only slightly better here
nrow(lu_BothW_scaled)
nrow(lu_BothR_scaled)
nrow(lu_OnlyA_scaled) + nrow(lu_OnlyB_scaled) + nrow(lu_BothW_scaled) + nrow(lu_BothR_scaled)
nrow(lu_predictions)


# ..4.3 Compare OrigRF vs H2O Model Normalized ----

lu_OnlyA_norm <- lu_predictions[lu_predictions$origVsH2ONorm == "OnlyA",] 
lu_OnlyB_norm <- lu_predictions[lu_predictions$origVsH2ONorm == "OnlyB",]
lu_BothW_norm <- lu_predictions[lu_predictions$origVsH2ONorm == "BOTH_WRONG",]
lu_BothR_norm <- lu_predictions[lu_predictions$origVsH2ONorm == "BOTH_RIGHT",]

nrow(lu_OnlyA_norm)
nrow(lu_OnlyB_norm)  # Normalized actually performed only slightly better here
nrow(lu_BothW_norm)
nrow(lu_BothR_norm)
nrow(lu_OnlyA_norm) + nrow(lu_OnlyB_norm) + nrow(lu_BothW_norm) + nrow(lu_BothR_norm)
nrow(lu_predictions)


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

slt2data_h2o_scaled <- slt2data_scaled
slt2data_h2o_scaled[,"Class"] <- as.logical(slt2data_h2o_scaled[,"Class"]) 
rfh2o_slt2_performance_scaled <- h2o.performance(rfh2o, newdata = as.h2o(slt2data_h2o_scaled))

ludata_h2o_scaled <- ludata_scaled
ludata_h2o_scaled[,"Class"] <- as.logical(ludata_h2o_scaled[,"Class"]) 
rfh2o_lu_performance_scaled <- h2o.performance(rfh2o, newdata = as.h2o(ludata_h2o_scaled))

slt2data_h2o_norm <- slt2data_norm
slt2data_h2o_norm[,"Class"] <- as.logical(slt2data_h2o_norm[,"Class"]) 
rfh2o_slt2_performance_norm <- h2o.performance(rfh2o, newdata = as.h2o(slt2data_h2o_norm))

ludata_h2o_norm <- ludata_norm
ludata_h2o_norm[,"Class"] <- as.logical(ludata_h2o_norm[,"Class"]) 
rfh2o_lu_performance_norm <- h2o.performance(rfh2o, newdata = as.h2o(ludata_h2o_norm))



# 2. Compare Metrics  ----
# ..2.1 Accuracy ----

metrics_table <- data.frame("Accuracy" = 1:8)
rownames(metrics_table) <- c("origRF_slt2_perf", "rfh2o_slt2_perf", "rfh2o_slt2_perf_scaled", "rfh2o_slt2_perf_norm",
                             "origRF_lu_perf", "rfh2o_lu_perf",  "rfh2o_lu_perf_scaled",  "rfh2o_lu_perf_norm")
metrics_table["origRF_slt2_perf","Accuracy"] <- nrow(slt2_predictions[slt2_predictions$origVsH2O == "BOTH_RIGHT" | 
                                                     slt2_predictions$origVsH2O == "OnlyA",]
                                    )/nrow(slt2_predictions)
metrics_table["rfh2o_slt2_perf","Accuracy"] <- nrow(slt2_predictions[slt2_predictions$origVsH2O == "BOTH_RIGHT" | 
                                                       slt2_predictions$origVsH2O == "OnlyB",]
                                    )/nrow(slt2_predictions)
metrics_table["rfh2o_slt2_perf_scaled","Accuracy"] <- nrow(slt2_predictions[slt2_predictions$origVsH2OScaled == "BOTH_RIGHT" | 
                                                       slt2_predictions$origVsH2OScaled == "OnlyB",])/nrow(slt2_predictions)

metrics_table["rfh2o_slt2_perf_norm","Accuracy"] <- nrow(slt2_predictions[slt2_predictions$origVsH2ONorm == "BOTH_RIGHT" | 
                                                        slt2_predictions$origVsH2ONorm == "OnlyB",])/nrow(slt2_predictions)

metrics_table["origRF_lu_perf","Accuracy"] <- nrow(lu_predictions[lu_predictions$origVsH2O == "BOTH_RIGHT" | 
                                                     lu_predictions$origVsH2O == "OnlyA",]
                                    )/nrow(lu_predictions)
metrics_table["rfh2o_lu_perf","Accuracy"] <- nrow(lu_predictions[lu_predictions$origVsH2O == "BOTH_RIGHT" | 
                                                     lu_predictions$origVsH2O == "OnlyB",]
                                    )/nrow(lu_predictions)
metrics_table["rfh2o_lu_perf_scaled","Accuracy"] <- nrow(lu_predictions[lu_predictions$origVsH2OScaled == "BOTH_RIGHT" | 
                                                     lu_predictions$origVsH2OScaled == "OnlyB",]
                                    )/nrow(lu_predictions)

metrics_table["rfh2o_lu_perf_norm","Accuracy"] <- nrow(lu_predictions[lu_predictions$origVsH2ONorm == "BOTH_RIGHT" | 
                                                     lu_predictions$origVsH2ONorm == "OnlyB",]
                                    )/nrow(lu_predictions)


metrics_table

## Accuracy graph based on different thresholds("cutoff")

plot(origRF_slt2_performance$acc, col = "blue") 
lines(h2o.accuracy(rfh2o_slt2_performance), type = "l", col = "red")
lines(h2o.accuracy(rfh2o_slt2_performance_scaled), type = "l", col = "green")
lines(h2o.accuracy(rfh2o_slt2_performance_norm), type = "l", col = "magenta")


plot(origRF_lu_performance$acc, col = "blue")
lines(h2o.accuracy(rfh2o_lu_performance), type = "l", col = "red")
lines(h2o.accuracy(rfh2o_lu_performance_scaled), type = "l", col = "green")
lines(h2o.accuracy(rfh2o_lu_performance_norm), type = "l", col = "magenta")

# ..2.2 AUCPR ----

metrics_table["origRF_slt2_perf","AUCPR"] <- origRF_slt2_performance$pr$auc.integral
metrics_table["rfh2o_slt2_perf", "AUCPR"] <- h2o.aucpr(rfh2o_slt2_performance)
metrics_table["rfh2o_slt2_perf_scaled", "AUCPR"] <- h2o.aucpr(rfh2o_slt2_performance_scaled)
metrics_table["rfh2o_slt2_perf_norm", "AUCPR"] <- h2o.aucpr(rfh2o_slt2_performance_norm)

metrics_table["origRF_lu_perf",  "AUCPR"] <- origRF_lu_performance$pr$auc.integral
metrics_table["rfh2o_lu_perf",  "AUCPR"] <- h2o.aucpr(rfh2o_lu_performance)
metrics_table["rfh2o_lu_perf_scaled", "AUCPR"] <- h2o.aucpr(rfh2o_lu_performance_scaled)
metrics_table["rfh2o_lu_perf_norm", "AUCPR"] <- h2o.aucpr(rfh2o_lu_performance_norm)

metrics_table

plot(origRF_slt2_performance$PR,  col = "blue", lwd = 2 )
lines(x = h2o.recall(rfh2o_slt2_performance)[,"tpr"],
      y = h2o.precision(rfh2o_slt2_performance)[,"precision"],
      col = "red", type = "l", lwd = 2
)
lines(x = h2o.recall(rfh2o_slt2_performance_norm)[,"tpr"],
      y = h2o.precision(rfh2o_slt2_performance_norm)[,"precision"],
      col = "magenta", type = "l", lwd = 2
)
lines(x = h2o.recall(rfh2o_slt2_performance_scaled)[,"tpr"],
      y = h2o.precision(rfh2o_slt2_performance_scaled)[,"precision"],
      col = "green", type = "l", lwd = 2
)


plot(origRF_lu_performance$PR,  col = "blue", lwd = 2 )
lines(x = h2o.recall(rfh2o_lu_performance)[,"tpr"],
      y = h2o.precision(rfh2o_lu_performance)[,"precision"],
      col = "red", type = "l", lwd = 2
)
lines(x = h2o.recall(rfh2o_lu_performance_scaled)[,"tpr"],
      y = h2o.precision(rfh2o_lu_performance_scaled)[,"precision"],
      col = "green", type = "l", lwd = 2
)
lines(x = h2o.recall(rfh2o_lu_performance_norm)[,"tpr"],
      y = h2o.precision(rfh2o_lu_performance_norm)[,"precision"],
      col = "magenta", type = "l", lwd = 2
)


# ..2.3 Sensitivy vs Specificity Graph ----

plot(origRF_slt2_performance$SS, col = "blue", lwd = 2)
#rfh2o_slt2_specificity <- h2o.specificity(rfh2o_slt2_performance)[,"tnr"] # 400 samples
#rfh2o_slt2_sensitivity <- h2o.sensitivity(rfh2o_slt2_performance)[,"tpr"] # 400 samples
lines(x = h2o.specificity(rfh2o_slt2_performance)[,"tnr"],
     y = h2o.sensitivity(rfh2o_slt2_performance)[,"tpr"],
     col = "red", type = "l", lwd = 2
     )

lines(x = h2o.specificity(rfh2o_slt2_performance_scaled)[,"tnr"],
      y = h2o.sensitivity(rfh2o_slt2_performance_scaled)[,"tpr"],
      col = "green", type = "l", lwd = 2
)
lines(x = h2o.specificity(rfh2o_slt2_performance_norm)[,"tnr"],
      y = h2o.sensitivity(rfh2o_slt2_performance_norm)[,"tpr"],
      col = "magenta", type = "l", lwd = 2
)



plot(origRF_lu_performance$SS, col = "blue", lwd = 2)
lines(x = h2o.specificity(rfh2o_lu_performance)[,"tnr"],
      y = h2o.sensitivity(rfh2o_lu_performance)[,"tpr"],
      col = "red", type = "l", lwd = 2
)

lines(x = h2o.specificity(rfh2o_lu_performance_scaled)[,"tnr"],
      y = h2o.sensitivity(rfh2o_lu_performance_scaled)[,"tpr"],
      col = "green", type = "l", lwd = 2
)

lines(x = h2o.specificity(rfh2o_lu_performance_norm)[,"tnr"],
      y = h2o.sensitivity(rfh2o_lu_performance_norm)[,"tpr"],
      col = "magenta", type = "l", lwd = 2
)


# ..2.4 Other H2O Metrics ----
h2o.confusionMatrix(rfh2o_slt2_performance)
h2o.confusionMatrix(rfh2o_slt2_performance_scaled)
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

# Conclusion from the model interpretations: The normalization, scaled, and the GLM created here were solely for the purpose of applying LIME to our orig RF.

# D) LOCAL METHODS ----
# Data to Analyze ----

head(slt2_OnlyA)
head(slt2_OnlyA_scaled)
head(slt2_OnlyA_norm)


# E) LIME ----
#
# LIME FUNCTION ARGUMENTS
# x	              The training data used for training the model that should be explained.
# model	          The model whose output should be explained
# bin_continuous  Should continuous variables be binned when making the explanation
# n_bins	        The number of bins for continuous variables if bin_continuous = TRUE
# quantile_bins	  Should the bins be based on n_bins quantiles or spread evenly over the range of
#                   the training data
# use_density     If bin_continuous = FALSE should continuous data be sampled using a kernel density 
#                 estimation. If not, continuous features are expected to follow a normal distribution.
#
# LIME SETTINGS JUSTIFICATIONS
# bin_continuous = TRUE --> Makes explanation easier, and having this setting as true helped in obtaining more consistant results.
# quantile_bins = FALSE --> setting this to TRUE generates an error b/c of the sameStrand and sameDownStrand don't have enough variance for the algorithm to use quantile bining.  
# use_density = TRUE --> in section A, 3.4, we showed that the features don't follow a normal distribution. Based on the explanation from LIME, we should have this setting as FALSE. 

# 1. LIME in original Model as is ----
# this section is expected to fail b/c the data is not normalized, or scaled. We simply tried applying LIME to the Orig RF hoping it would work like a plug-and-play kind of application. To get LIME to work on the "Original" RF model, you must retrain the model with the data normalized, or scaled. 
#lime_explainer_orig <- lime(as.data.frame(trainData[, c(1:7)]), origRF,bin_continuous = TRUE, quantile_bins=FALSE, use_density = TRUE)
#trainData[1,c(1:7)]
#lime_explanations_orig <- explain(as.data.frame(trainData[1,c(1:7)]), lime_explainer_orig, n_labels = 1, n_features = 7, n_permutations = 2000)
#plot_features(lime_explanations_orig)

# 2. LIME to H2O RF Model ----
# Here, the training and testing data is not normalized or scaled. 
# ..2.1 H2O RF Model, density = T ----
sampleData <- slt2data_h2o[1,]  # A true example
sampleData[c(2:4),] <- slt2data_h2o[1,]
sampleData[5,] <- slt2data_h2o[1986,] # A false example
sampleData[c(6:8),] <- slt2data_h2o[1986,]
sampleData
h2o.predict(object = rfh2o, newdata = as.h2o(sampleData[c(1,8),c(1:7)]) )
lime_explainer_rfh2o <- lime( as.data.frame( trainData[,c(1:7)] ), rfh2o,
                              bin_continuous = FALSE, 
                              quantile_bins = FALSE,
                              use_density = TRUE
)
lime_explanations_rfh2o <- explain( as.data.frame(sampleData[,c(1:7)] ),   # Data to explain
                                    lime_explainer_rfh2o,     # Explainer to use
                                    n_labels = 1, # only 1 type of category
                                    n_features = 7, # Number of features we want to use for explanation
                                    n_permutations = 2000,
                                    dist_fun = "gower"  
)
plot_features(lime_explanations_rfh2o) 

# ..2.2 H2O RF Model, density = F ----

lime_explainer_rfh2o <- lime( as.data.frame( trainData[,c(1:7)] ), rfh2o,
                              bin_continuous = TRUE, quantile_bins = FALSE,
                              use_density = TRUE, n_bins = 10
)
lime_explanations_rfh2o <- explain( as.data.frame(sampleData[,c(1:7)] ),  
                                    lime_explainer_rfh2o, n_labels = 1, 
                                    n_features = 7, n_permutations = 2000,
                                    dist_fun = "gower" 
)
plot_features(lime_explanations_rfh2o) 

# 3. LIME to H2O RF scaled Model ----
# ..3.1 H2O RF scaled Model, density = T ----
sampleData_scaled <- slt2data_h2o_scaled[1,]  # A true example
sampleData_scaled[c(2:4),] <- slt2data_h2o_scaled[1,]
sampleData_scaled[5,] <- slt2data_h2o_scaled[1986,] # A false example
sampleData_scaled[c(6:8),] <- slt2data_h2o_scaled[1986,]
sampleData_scaled
h2o.predict(object = rfh2o_scaled, newdata = as.h2o(sampleData_scaled[c(1,8),c(1:7)]))
lime_explainer_rfh2o_scaled <- lime( as.data.frame( trainData_scaled[,c(1:7)] ), rfh2o_scaled,
                              bin_continuous = FALSE, auantile_bins = FALSE, use_density = TRUE) 
lime_explanations_rfh2o_scaled <- explain( as.data.frame(sampleData_scaled[,c(1:7)]),  
                                    lime_explainer_rfh2o_scaled, n_labels = 1, 
                                    n_features = 7, n_permutations = 2000, dist_fun = "gower" )
plot_features(lime_explanations_rfh2o_scaled)

# ..3.2 H2O RF scaled Model, density = F ----
lime_explainer_rfh2o_scaled <- lime( as.data.frame( trainData_scaled[,c(1:7)] ), rfh2o_scaled,
                                     bin_continuous = TRUE, quantile_bins = FALSE,
                                     use_density = FALSE, n_bins = 10)
lime_explanations_rfh2o_scaled <- explain( as.data.frame(sampleData_scaled[,c(1:7)]),   # Data to explain
                                           lime_explainer_rfh2o_scaled,     # Explainer to use
                                           n_labels = 1, n_features = 7, 
                                           n_permutations = 2000, dist_fun = "gower")
plot_features(lime_explanations_rfh2o_scaled)

# 4. LIME to H2O RF normalized Model ----
# ..4.1 H2O RF Normalized Model, density = T ----
sampleData_norm <- slt2data_h2o_norm[1,]  # A true example
sampleData_norm[c(2:4),] <- slt2data_h2o_norm[1,]
sampleData_norm[5,] <- slt2data_h2o_norm[1986,] # A false example
sampleData_norm[c(6:8),] <- slt2data_h2o_norm[1986,]
sampleData_norm
h2o.predict(object = rfh2o_norm, newdata = as.h2o(sampleData_norm[c(1,8),c(1:7)]) )
lime_explainer_rfh2o_norm <- lime( as.data.frame( trainData_norm[,c(1:7)] ), rfh2o_norm,
                                   bin_continuous = FALSE, quantile_bins = FALSE, use_density = TRUE)
lime_explanations_rfh2o_norm <- explain( as.data.frame(sampleData_norm[,c(1:7)] ),
                                         lime_explainer_rfh2o_norm, n_labels = 1, 
                                         n_features = 7, n_permutations = 2000, dist_fun = "gower" )
plot_features(lime_explanations_rfh2o_norm)

# ..4.2 Applying LIME to Normalized Model ----
lime_explainer_rfh2o_norm <- lime( as.data.frame( trainData_norm[,c(1:7)] ), rfh2o_norm,
                                     bin_continuous = TRUE, quantile_bins = FALSE,
                                     use_density = FALSE,  n_bins = 10)
lime_explanations_rfh2o_norm <- explain( as.data.frame(sampleData_norm[,c(1:7)] ), lime_explainer_rfh2o_norm,   
                                         n_labels = 1,  n_features = 7, 
                                         n_permutations = 2000, dist_fun = "gower"  
                                         )

plot_features(lime_explanations_rfh2o_norm)

# ..1.5 Applying LIME to GLM Model ----

lime_explainer_glmh2o <- lime( as.data.frame( trainData[,c(1:7)] ),glmh2o, 
                               bin_continuous = FALSE, quantile_bins = FALSE, 
                               use_density = FALSE, n_bins = 10 )
lime_explanations_glmh2o <- explain( as.data.frame(sampleData[,c(1:7)] ), lime_explainer_glmh2o,   
                                         n_labels = 1,  n_features = 7, 
                                         n_permutations = 2000, dist_fun = "gower"  )
plot_features(lime_explanations_glmh2o)

# LIME CONCLUSIONS ----
#
# A known problem with LIME is that explanations may be inconsistant from one instance to the
# next, even when the instances are very similar to each. With this in mind, each instance
# we tested was tested at least 4 times (ie we manually ran this section over and over at
# least 4 times), and with a high enough number of permutations,
# we managed to get somewhat consistent results. However, some of the results obtained from
# LIME were heavily flawed as sometimes it would yield explanations that contradicted the result. 
# Depending on the data to be explain, it would sometimes leave out sameStrand and sameDownStrand from the 
# explanation.
# Additionally, a GLM was also trained to attempt to measure the variance of the data, but the
# results from the GLM yielded an R^2 value between 0.2 and 0.3, with a max recall of 0.07. In
# other words, the GLM model was only predicting false instances, and since LIME uses linear 
# models for its explanations, based on our observations, it is safe to conclude that LIME
# was doing the same. Whenever a TRUE instance was analyze, LIME's explanations would say 
# that most factors contradict the prediction, which is obviously counter intuitive and useless.
#
# During the making of this project, several problems were encountered while trying
# to apply LIME to the sRNA RF model. LIME attempts to explain a complex model's behavior
# in a small region corresponding to a particular instance of interest by applying a 
# linear model. However, in this particular use case, the mixture of categorical and
# numerical values, plus the differences in ranges within each feature resulted in various
# experinments without consistent results. Fr this particular use case, LIME was implement 
# in models trained with the original training data, scaled training data, and normalized
# training data. Several attempts were also made at modifying the arguments for the lime() 
# function, such as all the possible permutations between bin_continuous, quantile_bins
# and use_density. In the explain function, the number of permutations was increased from
# 10 to 5000, the Euclidean and Gower distance functions were tested, and kernel_widths 
# ranging from 0.1 to 5 were tried. 0.01 was also tried as a kernel_width, but generated
# errors in the LIME function, and anything above 5 just seemed exaggerated. 


# F) SHAP Values ----

# 1. something ... ----
?partialPlot()
SHAP_H2O1 <- h2o.predict_contributions(rfh2o_scaled, as.h2o(sampleData_scaled[1,c(1:7)]))

SHAP_H2O1


# G) PDP ----

# 1. Generate PDPs for the RF models ----

?partialPlot()
?h2o.partialPlot
h2o.partialPlot(rfh2o_scaled, data = as.h2o(slt2data_h2o_scaled), cols = "SS", plot=TRUE, nbins=2)

h2o.partialPlot(rfh2o, data = as.h2o(ludata_h2o), cols = "SS")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "SS")

h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "Pos10wrtsRNAStart")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "DistTerm")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "Distance")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "DownDistance")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "sameStrand")
h2o.partialPlot(rfh2o, data = as.h2o(trainData), cols = "sameDownStrand")

#TODO - PDPs for classification tasks, instead of regression
#TODO - See how LIME and SHAP agree and compare

#### QUESTIONS ----
# Why are the performance metrics so different when scaling or normalizing the data,
# but the accuracy is still extremely similar?


