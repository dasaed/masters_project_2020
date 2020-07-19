# CITATIONS ----
clc <- function() cat("\014") ; clc()
remove(list = ls())

library("h2o")          # ML model building
library("ROCR")         # ML evaluation
library("PRROC")        # ML evaluation
library("randomForest") # ML model building
library("randomForestExplainer") # ML Global interpretation
library("lime")         # ML local interpretation
library("shapper")      # ML local interpretation
library("png") 

citation("h2o")          # ML model building
citation("ROCR")         # ML evaluation
citation("PRROC")        # ML evaluation
citation("randomForest") # ML model building
citation("randomForestExplainer") # ML Global interpretation
citation("lime")         # ML local interpretation
citation("shapper")
citation("png")

