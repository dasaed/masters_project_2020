# Installing the required libraries ----
## Required Libraries: randomForest, h2o

install.packages("randomForest")
# Installing h2o, according to their website
# "https://docs.h2o.ai/h2o/latest-stable/h2o-docs/downloading.html#install-in-r"
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))
# test the h2o library
library(h2o)
localH2O = h2o.init() 
demo(h2o.kmeans)

# initialize the libraries ----
library("randomForest")
library("h2o")


# ORIGINAL SCRIPT ----
#!/usr/bin/env Rscript
library("optparse")
library(randomForest)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input data (feature table obtained by sRNACharP)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="outPredictions.txt", 
              help="Filename to output predictions [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least the input data filename must be provided", call.=FALSE)
}

RF <- readRDS("RF_classifier4sRNA.rds")
data <- read.table(opt$input, sep = "\t", header = TRUE)

pred <- predict(RF, data, type = "prob")
colnames(pred) <- c("No sRNA", "sRNA")
write.table(pred, file = opt$out, sep = "\t", col.names = TRUE, row.names = TRUE)

# Running the Model from within script ----
inputFile <- "FeatureTable.tsv" # the sample file
outputFile <- paste( "output_from_",inputFile, sep="" )
RF <- readRDS("RF_classifier4sRNA.rds")
data <- read.table( inputFile, sep = "\t", header = TRUE)
data
class(data)
pred <- predict(RF, data, type = "prob")
colnames(pred) <- c("No sRNA", "sRNA") # adding header columns to the prediction output
pred
class(pred)
write.table(pred, file = outputFile, sep = "\t", col.names = TRUE, row.names = TRUE)

# Running the Model with Single Data Point ----


singlePoint <- "singlePoint.tsv"
headers <- c("SS","Pos10wrtsRNAStart","DistTerm","Distance","sameStrand","DownDistance","sameDownStrand")
write.csv(headers, singlePoint, eol="\t", row.names = FALSE)
cat("\n", file=singlePoint, append=TRUE)

#instanceToTest <- list("sRNA00822",-43.8,-36,1000,-38,1,153,0)
instanceToTest <- list(-43.8,-36,1000,-38,1,153,0) 
instanceToTestDf <- as.data.frame(instanceToTest)
write.table(instanceToTestDf, singlePoint, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

myData <- read.table(singlePoint, sep="\t", header = TRUE)
myData
singlePred <- predict(RF, myData, type = "prob")

# Applying the model with a single Data Point - Works, but you lose your labels

# headers <- c("Pos10wrtsRNAStart","DistTerm","Distance","sameStrand","DownDistance","sameDownStrand")
# instanceToTest <- list("sRNA00822",-43.8,-36,1000,-38,1,153,0)

RF <- readRDS("RF_classifier4sRNA.rds") # Load the model

instanceToTest <- list(-43.8,-36,1000,-38,1,153,0) # Load the instance
predictSingleInstance <- function( singleInstance) {
  singlePred <- predict(RF, instanceToTest, type = "prob") # Make Prediction
  colnames(singlePred) <- c("No sRNA", "sRNA") # adding header columns to the prediction output
  rownames(singlePred) <- c("sRNA00822")
  return(singlePred)
}
myPrediction <- predictSingleInstance(instanceToTest)
myPrediction


# LIME Test ----
install.packages("lime")
library("lime")
library("randomForest")
RF <- readRDS("RF_classifier4sRNA.rds") # Load the model

origTrainingData <- read.csv( "training_combined.csv", header = TRUE, sep = ",") # load Orig Training data
class(origTrainingData)
head(origTrainingData)
length(origTrainingData)

origTrainingDataLabels <- read.csv( "training_combined_labels.csv", header = TRUE, sep = "," ) # load Orig Training data labes
class(origTrainingDataLabels)
head(origTrainingDataLabels)
classification <- origTrainingDataLabels$Class
trainingDataWithLabels <- cbind(origTrainingData, classification)

# instances to explain ----
inputFile <- "FeatureTable.tsv"
testData <- read.table( inputFile, sep = "\t", header = TRUE)
class(testData)
head(testData)
length(testData)
testDataPredictions <- predict(RF, testData, type="prob")
testDataPredictions

?lime()

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

origTrainingData
lime_explainer <- lime( origTrainingData,             # Original training data
                        RF,                       # The model to explain
                        bin_continuous = TRUE,    # Should continuous variables be binned 
                                                  # when making the explanation
                        quantile_bins = FALSE     # Should the bins be based on n_bins quantiles
                                                  # or spread evenly over the range of the training data
                        )
lime_explanations <- explain( testData,           # Data to explain
                              lime_explainer,     # Explainer to use
                              n_labels = 1,       # The number of labels to explain
                              n_features = 4     # Number of permutations to use for each explanation
                              )



# Applying LIME using h2o----
library("h2o")
h2o.init()
#h2o.no_progress() 


# Clean up  ----
gc() # gc = garbage collector
