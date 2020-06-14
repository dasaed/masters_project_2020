#!/usr/bin/env Rscript
# CLEAR CONSOLE AND ENVIROMENT ----

clc <- function() cat("\014") ; clc()
remove(list = ls())

# INSTANLL AND CHANGE PACKAGES ----
#uncomment the following lines if you need to install the libraries
#install.packages("randomForest")
#install.packages("ROCR")
#install.packages("PRROC")
#install.packages("lime")

library("randomForest") # To train the model
library("lime")         # To evaluate model's results
library("ROCR")         # For evaluation purposes 
library("PRROC")        # For evaluation purposes

# IMPORT DATA ----

# Training Data
fullTrainingDataSet <- read.csv("./tests/combinedData.csv", header = TRUE)
fullTrainingDataSet[,"Class"] <- as.logical(fullTrainingDataSet[,"Class"])

dataSetTrainX <- fullTrainingDataSet[,-(8:9)]
dataSetTrainY <- fullTrainingDataSet[,(8:9)]

# Testing Data
slt2dataPos <- read.csv("./testing_datasets/SLT2_Positives.tsv", sep = "\t", header = TRUE)
slt2dataNeg <- read.csv("./testing_datasets/SLT2_Negatives.tsv", sep = "\t", header = TRUE)
slt2data <- rbind(slt2dataPos,slt2dataNeg)

ludataPos <- read.csv("./testing_datasets/Lu_Positives.tsv", sep = "\t", header = TRUE)
ludataNeg <- read.csv("./testing_datasets/Lu_Negatives.tsv", sep = "\t", header = TRUE)
lu2data <- rbind(ludataPos,ludataNeg)


# TRAIN MODEL ----

# Train Model as Orig but with scaled features ----

# combData <- read.table("./tests/combinedData2.csv", sep =",", header = TRUE)

# dim(combData)
# head(combData)
# combData2 <- combData[,-8]#remove IDs
# head(combData2)



set.seed(1234)
tuneRF(myDataSetTrainX, y = factor(myDataSetTrainY[,"Class"]), ntreeTry = 400, mtryStart = 2) # looks for best parameters
RF_orig <- randomForest(x = myDataSetTrainX, 
                        y = factor(myDataSetTrainY[,"Class"]), 
                        mtry = 4, 
                        ntree = 400, 
                        importance = TRUE)

tuneRF(myDataSetTrainX_scaled, y = factor(myDataSetTrainY[,"Class"]), ntreeTry = 400, mtryStart = 2) # looks for best parameters
RF_stan <- randomForest(x = myDataSetTrainX_scaled, 
                        y = factor(myDataSetTrainY[,"Class"]), 
                        mtry = 4, 
                        ntree = 400, 
                        importance = TRUE)

# APPLY LIME ----
# Create necessary functions for randomForest's forest to work with LIME ----
#From Lime manual Return a data.frame in the case of predict_model(). If type = 'raw' it will contain one column named 'Response' holding the predicted values. 
#If type = 'prob' it will contain a column for each of the possible classes named after the class, each column holding the probability score for class membership. 
#For model_type() a character string. Either 'regression' or 'classification' is currently supported

# CREATE FUNCTIONS FOR LIME TO WORK WITH OUR RANDOM FOREST ----
model_type.randomForest <- function(x,...) 'classification'

predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, type = "prob")
  switch(
    type,
    raw = data.frame(Response = ifelse(as.vector(res[,2]) > 0.5, "1", "0"), stringsAsFactors = FALSE),
    prob = as.data.frame(res, check.names = F) 
  )
}

# Apply LIME to "Orig RF" ----

lime_explainer_scaled <- lime( as.data.frame(myDataSetTrainX_scaled),         # Original training data
                        RF_orig,                    # The model to explain
                        bin_continuous = TRUE, # Should continuous variables be binned 
                        # when making the explanation
                        n_bins = 2,            # The number of bins for continuous variables 
                        # if bin_continuous = TRUE
                        quantile_bins = FALSE,  # Should the bins be based on n_bins quantiles
                        # or spread evenly over the range of the training data
                        use_density = TRUE
)

lime_explainer_orig <- lime( as.data.frame(myDataSetTrainX),         # Original training data
                               RF_orig,                    # The model to explain
                               bin_continuous = TRUE, # Should continuous variables be binned 
                               # when making the explanation
                               n_bins = 2,            # The number of bins for continuous variables 
                               # if bin_continuous = TRUE
                               quantile_bins = FALSE,  # Should the bins be based on n_bins quantiles
                               # or spread evenly over the range of the training data
                               use_density = TRUE
)


# head(predictions)
colnames(myDataSetTrainX_scaled)
myDataSetTrainX_scaled[1:3,]
colnames(myDataSetTrainX)
myDataSetTrainX[1:3,]
lime_explanations_orig <- explain( as.data.frame(myDataSetTrainX[1:3,]),           # Data to explain
                              lime_explainer_orig,     # Explainer to use
                              n_labels = 1,
                              n_features = 7,
                              n_permutations = 100,
                              feature_select = "none"
)

lime_explanations_scaled <- explain( as.data.frame(myDataSetTrainX_scaled[1:3,]),           # Data to explain
                              lime_explainer_scaled,     # Explainer to use
                              n_labels = 1,
                              n_features = 7,
                              n_permutations = 100,
                              feature_select = "none"
)

origplot <- plot_features(lime_explanations_orig)
scaledplot <- plot_features(lime_explanations_scaled)
origplot
scaledplot

# Evaluate model's performance
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



# Clean up  ----
gc() # gc = garbage collector

