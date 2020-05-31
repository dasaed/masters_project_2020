library("lime")
library("randomForest")
RF <- readRDS("RF_classifier4sRNA.rds") # Load the model

origTrainingData <- read.csv( "training_combined.csv", header = TRUE, sep = ",") # load Orig Training data

origTrainingDataLabels <- read.csv( "training_combined_labels.csv", header = TRUE, sep = "," ) 
                                                        # load Orig Training data labes
Classification <- origTrainingDataLabels$Class
origTrainingDataWithLabels <- cbind(origTrainingData, Classification)

# instances to explain ----
inputFile <- "FeatureTable.tsv"
testData <- read.table( inputFile, sep = "\t", header = TRUE)
class(testData)

testDataPredictions <- predict(RF, testData, type="prob")
testDataPre
# randomForest
# RF <- readRDS("RF_classifier4sRNA.rds")
# pred <- predict(RF, data, type = "prob")

predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, ...)
  switch(
    type,
    raw = data.frame(Response = res$class, stringsAsFactors = FALSE),
    prob = as.data.frame(res["posterior"], check.names = FALSE)
  )
}

model_type.randomForest <- function(x, ...) 'classification'

?lime()
lime_explainer <- lime( origTrainingData,      # Original training data
                        RF,                    # The model to explain
                        bin_continuous = TRUE, # Should continuous variables be binned 
                                               # when making the explanation
                        n_bins = 5,           # The number of bins for continuous variables 
                                               # if bin_continuous = TRUE
                        quantile_bins = FALSE  # Should the bins be based on n_bins quantiles
                                               # or spread evenly over the range of the training data
                        )
lime_explanations <- explain( testData,           # Data to explain
                              lime_explainer,     # Explainer to use
                              n_labels = 7,
                              n_features = 7,
                              n_permutations = 10,
                              feature_select = "none"
                            )
lime_explanations
