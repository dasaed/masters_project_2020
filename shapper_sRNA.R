# A SHAP implementation exist for R, but in reality, it uses the python library. As a presetup for this script to work, you must run 'pip3 install shap' or 'pip install shap', to first fetch the python library.

#install.packages("shapper")

library("shapper")
library("randomForest")

trainDataSet <- read.csv("./DataSets/combinedData.csv", header = TRUE)
trainDataSet[,"Class"] <- as.logical(trainDataSet[,"Class"]) 

slt2dataPos <- read.csv("./DataSets/SLT2_Positives.tsv", sep = "\t", header = TRUE)
slt2dataPos$Class <- rep(1,nrow(slt2dataPos))
slt2dataNeg <- read.csv("./DataSets/SLT2_Negatives.tsv", sep = "\t", header = TRUE)
slt2dataNeg$Class <- rep(0,nrow(slt2dataNeg))
slt2data <- rbind(slt2dataPos,slt2dataNeg)

ludataPos <- read.csv("./DataSets/Lu_Positives.tsv", sep = "\t", header = TRUE)
ludataPos$Class <- rep(1,nrow(ludataPos))
ludataNeg <- read.csv("./DataSets/Lu_Negatives.tsv", sep = "\t", header = TRUE)
ludataNeg$Class <- rep(0,nrow(ludataNeg))
ludata <- rbind(ludataPos,ludataNeg)

sampleData <- slt2data[1,]  # A true sample
sampleData[c(2:4),] <- slt2data[1,]
sampleData[5,] <- slt2data[1986,] # A false sample
sampleData[c(6:8),] <- slt2data[1986,]
sampleData

origRF <- readRDS("RF_classifier4sRNA.rds")

p_function <- function(model, data) predict(model, newdata = data, type = "prob")

ive_rf <- individual_variable_effect(origRF, data = trainDataSet[,c(1:7)], predict_function = p_function,
                                     new_observation = sampleData[1,-8], nsamples = 50)

ive_rf
plot(ive_rf, bar_width = 4)
