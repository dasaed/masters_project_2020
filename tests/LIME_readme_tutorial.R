#!/usr/bin/env Rscript
clc <- function() cat("\014") ; clc()
remove(list = ls())

install.packages("caret")
install.packages('e1071', dependencies=TRUE)
install.packages("lime")
library("caret")
library("lime")
library("datasets")
data(iris)
summary(iris)
head(iris)
View(iris)
class(iris)

iris_test <- iris[1:5, 1:4] # only the first 5 rows without the last column -> row x col
iris_test
iris_train <- iris[-(1:5), 1:4] # all data except first 5 rows and last column -> row x col
head(iris_train)
iris_lab <- iris[[5]][-(1:5)] # only the fifth column, without the first 5 rows 
                              #   -> extract data from data frame -> col x rows
head(iris_lab)
?train
model <- train(iris_train, iris_lab, method='rf')
model
explainer <- lime(iris_train, model)
explainer
class(iris_test)
explanation <-  explain(iris_test, explainer, n_labels = 1, n_features = 2)
explanation
plot_features(explanation)
?lime
