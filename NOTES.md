# NOTES
These are just personal notes, for myself to use as guidance should they be needed.
I would not read into them as they are pieces of code I used at one point or another,
and no longer need them in my final code. However, I want to keep them as reference
should I need to come back to them at some point.

## DATA PRE-PROCESSING ----
## Optional --> this is just to standarize the data, but it could drastically change results 
myDataSetTrainX_scaled <-  scale(myDataSetTrainX[,-c(5,7)], center = TRUE, scale = TRUE)
myDataSetTrainX_scaled  <- cbind(myDataSetTrainX_scaled[,(1:4)], 
                                 myDataSetTrainX[,5],
myDataSetTrainX_scaled[,5],
                                 myDataSetTrainX[,7])
colnames(myDataSetTrainX_scaled) <- colnames(myDataSetTrainX)

\# sameStrand -> not scale
\# sameDownStrand -> not scale


------------------
# Train Model as Orig but with scaled features ----

# combData <- read.table("./tests/combinedData2.csv", sep =",", header = TRUE)

 dim(combData)
 head(combData)
 combData2 <- combData[,-8]#remove IDs
 head(combData2)
# Train Model as Orig but with scaled features ----

 combData <- read.table("./tests/combinedData2.csv", sep =",", header = TRUE)

 dim(combData)
 head(combData)
 combData2 <- combData[,-8]#remove IDs
 head(combData2)

## DATA PRE-PROCESSING ----
# We need to standarize(scale and center) our features, except for sameStrand and sameDownStrand
# to be able to use LIME in our models
# sameStrand -> not scale
# sameDownStrand -> not scale



# data preprocessing with h2o   
dataSetTrainX_scaled <-  scale(dataSetTrainX[,-c(5,7)], center = TRUE, scale = TRUE)
dataSetTrainX_scaled  <- cbind(dataSetTrainX_scaled[,(1:4)], 
                               dataSetTrainX[,5],
                               dataSetTrainX_scaled[,5], 
                               dataSetTrainX[,7])
colnames(dataSetTrainX_scaled) <- colnames(dataSetTrainX)

dataSetTrainX_for_h2o <- cbind(dataSetTrainX_scaled,dataSetTrainY[,2])
colnames(dataSetTrainX_for_h2o)[8] <- "Class"
head(dataSetTrainX_for_h2o,4)


# rf2 > uses the normal data, without standarization
rf2 <- h2o.randomForest( 
  training_frame = trainData_orig,
  x = 1:7,                   # features to use to generate the prediction
  y = 9,                     # Class type -> what we want to predict
  model_id = "rf1_sRNA",     # name of model in h2o
  ntrees = 400,              # max number of trees  
  seed = 1234,               # seed, has to be set WITHIN the h2o function
  # and it's supposed to be different from "R's seed", so 
  # results might not be exactly the same as orig model, but
  # should be similar enough
  mtries = 2,                 # Same as original model 
  max_depth = 20
)

rf2@model$variable_importances
h2o.varimp_plot(rf2)
h2o.performance(rf2)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

> rfh2o_peformance
H2OBinomialMetrics: drf
** Reported on training data. **
** Metrics reported on Out-Of-Bag training samples **

MSE:  0.08009398
RMSE:  0.2830088
LogLoss:  0.2681074
Mean Per-Class Error:  0.1267894
AUC:  0.9391835
AUCPR:  0.8614903
Gini:  0.878367
R^2:  0.5728321

Confusion Matrix (vertical: actual; across: predicted) for F1-optimal threshold:
       FALSE TRUE    Error     Rate
FALSE    443   46 0.094070  =46/489
TRUE      26  137 0.159509  =26/163
Totals   469  183 0.110429  =72/652

Maximum Metrics: Maximum metrics at their respective thresholds
                        metric threshold      value idx
1                       max f1  0.435400   0.791908 164
2                       max f2  0.185289   0.845384 223
3                 max f0point5  0.540218   0.809906 131
4                 max accuracy  0.540218   0.897239 131
5                max precision  0.997085   1.000000   0
6                   max recall  0.000530   1.000000 397
7              max specificity  0.997085   1.000000   0
8             max absolute_mcc  0.435400   0.719317 164
9   max min_per_class_accuracy  0.303686   0.871166 186
10 max mean_per_class_accuracy  0.354028   0.876278 177
11                     max tns  0.997085 489.000000   0
12                     max fns  0.997085 162.000000   0
13                     max fps  0.000085 489.000000 399
14                     max tps  0.000530 163.000000 397
15                     max tnr  0.997085   1.000000   0
16                     max fnr  0.997085   0.993865   0
17                     max fpr  0.000085   1.000000 399
18                     max tpr  0.000530   1.000000 397

Gains/Lift Table: Extract with `h2o.gainsLift(<model>, <data>)` or `h2o.gainsLift(<model>, valid=<T/F>, xval=<T/F>)`
> 



rfh2o@model$variable_importances
h2o.varimp_plot(rfh2o)
rfh2o_peformance <- h2o.performance(rfh2o)
plot(h2o.accuracy(rfh2o_peformance))
