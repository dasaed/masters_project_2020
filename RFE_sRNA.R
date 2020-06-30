#!/usr/bin/env Rscript
# RFE = Random Forest Explainer
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
# install.packages("randomForestExplainer")
library("randomForestExplainer")  # ML interpretation
library("ROCR")         # ML evaluation
library("PRROC")        # ML evaluation
library("randomForest") # ML model building




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
dataSetTrain <- read.csv("./tests/combinedData.csv", header = TRUE)
#dataSetTrain[,"Class"] <- as.logical(dataSetTrain[,"Class"]) # helps guarantee that our model trains the random forest (RF) as a classification tree and not a regression tree

dataSetTrainX <- dataSetTrain[,-(8:9)] 
dataSetTrainY <- dataSetTrain[,(8:9)]  



# B) Load RF models ----

# The original article, titled "Prioritizing bona fide bacterial small RNAs with machine learning 
# classifiers", and the source materials(training data, testing data, R scripts, .rds file, etc.) 
# can all be found in PeerJ under the following url: https://peerj.com/articles/6304/
# The original .rds file can be found in github under the following link:
# https://github.com/BioinformaticsLabAtMUN/sRNARanking


## 1. Load/Retrain Orig Model with Scaled and Norm Data ----

set.seed(1234)

tuneRF(dataSetTrain[,c(1:7)], y = factor(dataSetTrain[,9]), ntreeTry = 400, mtryStart = 2) 


origRF_rds <- readRDS("RF_classifier4sRNA.rds") # The original RF model doesn't have a "formula" defined within its functions/variables, and therefore, some of the functionality offered by the library is missing. 

dataSetFull <- dataSetTrain[,-8]
dataSetFull$Class <- as.factor(dataSetTrain$Class)
str(dataSetFull)
head(dataSetFull,2)
origRF <- randomForest( formula = Class ~.,
                        data = dataSetFull, # To retrain original model 
                        #x = dataSetTrain[,c(1:7)], # To retrain original model 
                        #y = factor(dataSetTrain[,9]), 
                        mtry = 2, ntree = 400, importance = TRUE, localImp = TRUE)
?randomForest()
#importance(origRF) # importance as provided the RF library (not RF explainer)

# 2. Generate Predictions ----
origRF_lu_pred <- predict(origRF, ludata[,-8], type = "prob")
origRF_slt2_pred <- predict(origRF, slt2data[,-8], type = "prob")

origRF_scaled_lu_pred <- predict(origRF_scaled, ludata_scaled[,-8], type = "prob")
origRF_scaled_slt2_pred <- predict(origRF_scaled, slt2data_scaled[,-8], type = "prob")

origRF_norm_lu_pred <- predict(origRF_norm, ludata_norm[,-8], type = "prob")
origRF_norm_slt2_pred <- predict(origRF_norm, slt2data_norm[,-8], type = "prob")


# This function is to simplify the comparison between the predictions from the Original RF model and the other models
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

# SLT2 comparisons ----
slt2_predictions <-  cbind(as.data.frame(slt2data),               # Input
                           as.data.frame(origRF_slt2_pred[,2]),   # Orig RF predictions 
                           as.data.frame(origRF_scaled_slt2_pred[,2]), # Scaled RF predictions
                           as.data.frame(origRF_norm_slt2_pred[,2])    # Normalized RF predictions
)
head(slt2_predictions, 2)
colnames(slt2_predictions) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", 
                                "sameStrand", "DownDistance", "sameDownStrand", "Class",
                                "OrigRF_T","origRF_scaled_T","origRF_norm_T")

head(slt2_predictions)
slt2_predictions$origPreds <- ifelse( slt2_predictions$OrigRF_T >= 0.5, 1, 0) 
slt2_predictions$origScaledPreds  <- ifelse( slt2_predictions$origRF_scaled_T >= 0.5, 1, 0) 
slt2_predictions$origNormPreds  <- ifelse( slt2_predictions$origRF_norm_T >= 0.5, 1, 0)  
slt2_predictions$origVsScaled <- NA
slt2_predictions$origVsNorm <- NA

# Check correlations between the OrigRF and the H2O models
# the closest to 1, the more correlated(similar) they are
cor(slt2_predictions[,"OrigRF_T"],slt2_predictions[,"origRF_scaled_T"]) 
cor(slt2_predictions[,"OrigRF_T"],slt2_predictions[,"origRF_norm_T"])


head(slt2_predictions,2)

for( i in 1:nrow(slt2_predictions) ){
  # Based on the way we are feeding the values to the compareAnswers function, 
  # Similiarities = OnlyA means that only the Orig RF got the right answer
  # Similiarities = OnlyB means that only the H2O RF got the right answer
  # All other answers will tell us where both models were right("BOTH_RIGHT"), or
  # wrong ("BOTH_WRONG")

  slt2_predictions[i,]$origVsScaled <- compareAnswers(
    slt2_predictions[i,]$origPreds,  
    slt2_predictions[i,]$origScaledPreds,  
    slt2_predictions[i,]$Class 
  )
  slt2_predictions[i,]$origVsNorm <- compareAnswers(
    slt2_predictions[i,]$origPreds,  
    slt2_predictions[i,]$origNormPreds,  
    slt2_predictions[i,]$Class 
  )
} 

head(slt2_predictions,2)
summary(slt2_predictions)
str(slt2_predictions)

# By looking at the difference between the number of predictions the Orig. RF and the Scaled RF,
# we can conclude that the models are similar enough in accuracy
slt2_OnlyA_scaled <- slt2_predictions[slt2_predictions$origVsScaled == "OnlyA",] 
slt2_OnlyB_scaled <- slt2_predictions[slt2_predictions$origVsScaled == "OnlyB",]
slt2_BothW_scaled <- slt2_predictions[slt2_predictions$origVsScaled == "BOTH_WRONG",]
slt2_BothR_scaled <- slt2_predictions[slt2_predictions$origVsScaled == "BOTH_RIGHT",]

nrow(slt2_OnlyA_scaled)
nrow(slt2_OnlyB_scaled)
nrow(slt2_BothW_scaled)
nrow(slt2_BothR_scaled)
nrow(slt2_OnlyA_scaled) + nrow(slt2_OnlyB_scaled) + nrow(slt2_BothW_scaled) + nrow(slt2_BothR_scaled)
nrow(slt2_predictions)


# By looking at the difference between the number of predictions the Orig. RF and the H2O RF,
# we can conclude that the models are similar enough in accuracy
slt2_OnlyA_norm <- slt2_predictions[slt2_predictions$origVsNorm == "OnlyA",] 
slt2_OnlyB_norm <- slt2_predictions[slt2_predictions$origVsNorm == "OnlyB",]
slt2_BothW_norm <- slt2_predictions[slt2_predictions$origVsNorm == "BOTH_WRONG",]
slt2_BothR_norm <- slt2_predictions[slt2_predictions$origVsNorm == "BOTH_RIGHT",]

nrow(slt2_OnlyA_norm)
nrow(slt2_OnlyB_norm)
nrow(slt2_BothW_norm)
nrow(slt2_BothR_norm)
nrow(slt2_OnlyA_norm) + nrow(slt2_OnlyB_norm) + nrow(slt2_BothW_norm) + nrow(slt2_BothR_norm)
nrow(slt2_predictions)




# C) RandomForest Explainer ----
# https://modeloriented.github.io/randomForestExplainer/articles/randomForestExplainer.html

### Feature Explanations
## "SS" = the free energy of the predicted secondary structure of the sRNA - mostly negative values
## "Pos10wrtsRNAStart" = distance to their closest predicted promoter site
## "DistTerm" = distance to their closest predicted Rho-independent terminator 
## "Distance", = distance to the closest reading frame on the LEFT("upstream") side
## "sameStrand" = boolean, if transcription is going in the same direction as the ORF(Left Open Reading Frame)
#                  --> ORF = genomic sequence that's supposed to code for a protein
## "DownDistance"  = distance to the closest reading frame on the RIGHT("downstream") side
## "sameDownStrand" =  boolean, if transcription is going in the same direction as the RORF(Right Open Reading Frame)
 

# 1. Obtain the distribution of minimal depth ----
#  * Note that the depth of a tree is equal to the length of the longest path from root to leave in this tree. This equals the maximum depth of a variable in this tree plus one, as leaves are by definition not split by any variable.
min_depth_frame_origRF <- min_depth_distribution(origRF) # Obtain the distribution of minimal depth
head(min_depth_frame_origRF,3)

#?min_depth_distribution() # Get minimal depth values for all trees in a random forest
#plot_min_depth_distribution(origRF, mean_sample = "all_trees") # Same as plotting the min depth frame, but slower
plot_min_depth_distribution(min_depth_frame_origRF, 
                            mean_sample = "all_trees" # "relevant_trees", "all_trees", or "top_trees"
                            #k = 7 # num of features to plot, 10 by default, but we only have 7, so this is optional as well
                            #mean_scale = TRUE # False by default
                            )
# * min depth in this package means what is the first time this variable is used to split the tree
# * There is an inverse relationship between times_a_root and mean_min_depth 


# 2. Variable Importance Measures ----
## Identify important variables
importance_frame_origRF <- measure_importance(origRF)
importance_frame_origRF
importance_frame_origRF[order(-importance_frame_origRF$accuracy_decrease), c("variable","accuracy_decrease")] # equivalent to graph generated by varImpPlot(origRF)
importance_frame_origRF[order(-importance_frame_origRF$gini_decrease), c("variable","gini_decrease")] # equivalent to graph generated by varImpPlot(origRF)

importance_frame_origRF[order(importance_frame_origRF$p_value), c("variable","p_value")] # This test tells us whether the observed number of successes (number of nodes in which Xj was used for splitting) exceeds the theoretical number of successes if they were random

#   sameDownStrand and sameStrand have very little importance in the model
# measure_importance() in the RF Explainer package = importance() and varImpPlot(), already include in the RF package. RF Explainer package just seems to present things prettier


# 3. Multi-Way Importance Plots ----
# plot_multi_way_importance(origRF, size_measure = "no_of_nodes") # gives the same result as below but takes longer
# Better presents the results from "importance_frame_origRF <- measure_importance(origRF)"

importance_frame_origRF
# These 2 plots use the same variables, but it different Axis
plot_multi_way_importance(importance_frame_origRF, size_measure = "no_of_nodes") # x = mean_min_depth, y = times_a_root
plot_multi_way_importance(importance_frame_origRF, 
                          x_measure = "mean_min_depth", 
                          y_measure = "no_of_nodes", 
                          size_measure = "times_a_root", 
                          no_of_labels = 7)
plot_multi_way_importance(importance_frame_origRF, 
                          x_measure = "p_value",
                          y_measure = "accuracy_decrease", 
                          size_measure = "gini_decrease"
                          # SS and Pos10wrtsRNAStart are insignificant according to the P-Value
                          )


# plot_multi_way_importance(importance_frame_origRF, size_measure = "p_value") # Generates a warning if p_value is used
plot_multi_way_importance(importance_frame_origRF, 
                          x_measure = "p_value", 
                          y_measure = "accuracy_decrease", 
                          size_measure = "times_a_root", 
                          no_of_labels = 7)

# 4. Measure Comparison using ggpairs  ----
## The diagonal are density plots ("advanced versions of histograms")
# plot_importance_ggpairs(origRF) # gives the same result as below but takes longer
importance_frame_origRF[order(-importance_frame_origRF$times_a_root),]
importance_frame_origRF[order(-importance_frame_origRF$p_value),] # a low p-value means that the feature is important for the prediction. 
# While Distance and DownDistance seem to be the most important variables, SS and Pos10wrtsRNAStart have the lowest p-values.

plot_importance_ggpairs(importance_frame_origRF)
plot_importance_ggpairs(importance_frame_origRF,
                        measures = c("gini_decrease","p_value", "no_of_nodes"))
plot_importance_ggpairs(importance_frame_origRF)


# plot_importance_rankings(origRF) # gives the same result as below but takes longer
plot_importance_rankings(importance_frame_origRF)

# (vars <- important_variables(forest, k = 5, measures = c("mean_min_depth", "no_of_trees"))) # gives the same result as below but takes longer
imp_vars <- important_variables(importance_frame_origRF, k = 5, measures = c("mean_min_depth", "gini_decrease"))
imp_vars
?min_depth_interactions
interactions_frame_origRF <- min_depth_interactions(
      forest = origRF,
      vars = imp_vars,
      mean_sample = "top_trees",
      uncond_mean_sample = "top_trees")
plot_min_depth_interactions(interactions_frame_origRF)
interactions_frame_origRF_2 <- min_depth_interactions(
      forest = origRF,
      vars = imp_vars,
      mean_sample = "relevant_trees",
      uncond_mean_sample = "relevant_trees")
plot_min_depth_interactions(interactions_frame_origRF_2)
interactions_frame_origRF_3 <- min_depth_interactions(
      forest = origRF,
      vars = imp_vars,
      mean_sample = "all_trees",
      uncond_mean_sample = "all_trees")
plot_min_depth_interactions(interactions_frame_origRF_3)


# 5. Plot Feature Interactions ----

head(dataSetTrain,3)
head(dataSetTrain[,c(1:7)]) 
str(dataSetTrainX)
dataSetTrainX$Distance <- as.numeric(dataSetTrainX$Distance)
?plot_predict_interaction
plot_predict_interaction(origRF, dataSetTrainX, "Distance", "SS") 
plot_predict_interaction(origRF, dataSetTrainX, "Distance", "Pos10wrtsRNAStart") 
plot_predict_interaction(origRF, dataSetTrainX, "DownDistance", "Pos10wrtsRNAStart") 
plot_predict_interaction(origRF, dataSetTrainX, "Distance", "SS") 
plot_predict_interaction(origRF, dataSetTrainX, "Distance", "DownDistance") 

?explain_forest()
explain_forest(origRF, data = dataSetTrain[,c(1:7)], interactions = TRUE)



# QUESTIONS ----
# Is the P-Value important in random forests?
# I get warning message saying "Using alpha for a discrete variable is not advised", whenever I use the size_measure = "p_value" in the multi-way importance plot function. I don't understand why? 
# Why is the difference in p-value between Distance and DownDistance so high? 

