#!/usr/bin/env Rscript
# CLEAR CONSOLE AND ENVIROMENT ----

clc <- function() cat("\014") ; clc()
remove(list = ls())

# INSTANLL AND CHANGE PACKAGES ----

install.packages("caret")
install.packages('e1071', dependencies=TRUE)
install.packages("lime")

library("h2o")
library("datasets")
library("randomForest")
library("caret")
library("lime")

# IMPORT DATA ----

myDataSetTrainX <- read.csv("training_combined.csv", header = TRUE)
myDataSetTrainY <- read.csv("training_combined_labels.csv", header = TRUE)
myDataSetTest <- read.csv("FeatureTable.tsv", sep="\t", header = TRUE)

summary(myDataSetTrainX)
dim(myDataSetTrainX)
head(myDataSetTrainX,5)

summary(myDataSetTrainY)
dim(myDataSetTrainY)
head(myDataSetTrainY,5)

summary(myDataSetTest)
dim(myDataSetTest)
head(myDataSetTest,5)

# DATA PRE-PROCESSING ----

myDataSetTrainX_scaled <-  scale(myDataSetTrainX[,-c(5,7)], center = TRUE, scale = TRUE)
myDataSetTrainX_scaled  <- cbind(myDataSetTrainX_scaled[,(1:4)], 
                                 myDataSetTrainX[,5],
                                 myDataSetTrainX_scaled[,5],
                                 myDataSetTrainX[,7])

head(myDataSetTrainX)
class(myDataSetTrainX)

head(myDataSetTrainX_scaled)
dim(myDataSetTrainX_scaled)
colnames(myDataSetTrainX_scaled) <- colnames(myDataSetTrainX)

# sameStrand -> not scale
# sameDownStrand -> not scale

# STATS CLASS BY JAVIER LAZARO ----
# scaled and centered
# centered = encontrar la media, y restarla a todos los datos, 
#           y divide todo en la desviacion estandar
# 
#summary(myDataSetTrainX)
#?scale
#myDataSetTrainXCentered <- scale(myDataSetTrainX, center = TRUE, scale = TRUE)
#head(myDataSetTrainXCentered)
#head(myDataSetTrainX)


# TRAIN MODEL ----

# Train Model as Orig ----

combData <- read.table("./tests/combinedData2.csv", sep =",", header = TRUE)

dim(combData)
head(combData)
combData2 <- combData[,-8]#remove IDs
head(combData2)

library(randomForest)
myDataSetTrainY
tuneRF(myDataSetTrainX_scaled, y = factor(myDataSetTrainY[,"Class"]), ntreeTry = 400, mtryStart = 2) # looks for best parameters

# oscillates btwn 2 and 4. Let's use 2 as in Python classifier
set.seed(1234)
RF_orig <- randomForest(x = myDataSetTrainX_scaled, y = factor(myDataSetTrainY[,"Class"]), mtry = 4, ntree = 400, importance = TRUE)

#x_scaled   <- scale(x_orig, center = TRUE, scale = TRUE)
#head(x_scaled)
#RF_scalced <- randomForest(x = x_scaled, y = factor(combData2[,8]), mtry = 2, ntree = 400, importance = TRUE)

# Train model with caret ----

# Train model with h2o ----



# Create model with default paramters
control <- trainControl(method="repeatedcv", number=10, repeats=3)
metric <- "Accuracy"
set.seed(seed)
mtry <- 2
tunegrid <- expand.grid(.mtry=mtry)
?expand.grid(.mtry=mtry)
rf_default <- train(Class~., data=dataset, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
?train()
getModelInfo(model="rf")
print(rf_default)

model <- train(myDataSetTrainX, myDataSetTrainY[,"Class"], method='rf' )
model


# APPLY LIME ----
# Create necessary functions for randomForest's forest to work with LIME ----
predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, ...)
  switch(
    type,
    raw = data.frame(Response = ifelse(res[,2] > 0.5, "sRNA", "notSRNA"), 
                     stringsAsFactors = FALSE
                     ),
    prob = res # as.data.frame(res, check.names = FALSE)
  )
  print(class(res))
  print(dim(res))
  print(res)
}

model_type.randomForest <- function(x, ...) 'classification'
# Apply LIME to "Orig RF" ----
?lime
?explain
lime_explainer <- lime( as.data.frame(myDataSetTrainX_scaled),         # Original training data
                        RF_orig,                    # The model to explain
                        bin_continuous = TRUE, # Should continuous variables be binned 
                        # when making the explanation
                        n_bins = 2,            # The number of bins for continuous variables 
                        # if bin_continuous = TRUE
                        quantile_bins = FALSE,  # Should the bins be based on n_bins quantiles
                        # or spread evenly over the range of the training data
                        use_density = TRUE
)

test3 <- myDataSetTrainX_scaled[1:3,]
dim(test3)
test3
head(myDataSetTrainX_scaled)
predictions <- predict(RF_orig, test3, type="prob")
predict_model.randomForest(RF_orig, test3, type="prob")
ifelse(predictions[,2] > 0.5, "sRNA", "notSRNA")

head(predictions)


?model_type
lime_explanations <- explain( as.data.frame(myDataSetTrainX_scaled),           # Data to explain
                              lime_explainer,     # Explainer to use
                              n_labels = 1,
                              n_features = 7,
                              n_permutations = 100,
                              feature_select = "none"
)

?explain

lime_explanations

plot_features(explanation)

#explainer <- lime(myDataSet_train, model)
#explainer
#explanation <-  explain(myDataSet_test, explainer, n_labels = 1, n_features = 2)
#explanation
#plot_features(explanation)
