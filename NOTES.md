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

SHOULD I WORRY?
> library("tidyverse")
── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
✔ ggplot2 3.3.1     ✔ purrr   0.3.4
✔ tibble  3.0.1     ✔ dplyr   1.0.0
✔ tidyr   1.1.0     ✔ stringr 1.4.0
✔ readr   1.3.1     ✔ forcats 0.4.0
── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::combine()  masks randomForest::combine()
✖ dplyr::explain()  masks lime::explain()
✖ dplyr::filter()   masks stats::filter()
✖ dplyr::lag()      masks stats::lag()
✖ ggplot2::margin() masks randomForest::margin()



slt2data[500,c(1:7)] # This is a good example to show LIME's inconsistencies










# V) RANDOM FOREST EXPLAINER

Random Forest Explainer (RFE) is a library that allows the user to understand the structure of a random forest by exploring its nodes, feature importance, and finding relationships between the main features. 
Most of the functions RFE has to generate graphs have the option to use the RF as an input, or a dataframe that already contains the values the function would have gotten from the RF. 
You can save the dataframe and simply load it whenever you want to generate the graphs to save some time. 
The methodology behind RFE is based around identifying the importance of the different features and their interactions so that a general idea of what the model is doing can be obtained. 

A lot of information and explanations can be found in the following [link](https://modeloriented.github.io/randomForestExplainer/articles/randomForestExplainer.html).

*REMINDER: Feature Explanations*
  * "SS" = the free energy of the predicted secondary structure of the sRNA - mostly negative values
  * "Pos10wrtsRNAStart" = distance to their closest predicted promoter site
  * "DistTerm" = distance to their closest predicted Rho-independent terminator 
  * "Distance", = distance to the closest reading frame on the LEFT("upstream") side
  * "sameStrand" = boolean, if transcription is going in the same direction as the ORF(Left Open Reading Frame)
     * ORF = genomic sequence that's supposed to code for a protein
  * "DownDistance"  = distance to the closest reading frame on the RIGHT("downstream") side
  * "sameDownStrand" =  boolean, if transcription is going in the same direction as the RORF(Right Open Reading Frame)


## 14. Minimal Depth Distribution

* "Note that the depth of a tree is equal to the length of the longest path from root to leave in this tree. This equals the maximum depth of a variable in this tree plus one, as leaves are by definition not split by any variable." [-CRAN Repository](https://cran.r-project.org/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html)

```{r}

min_depth_frame_origRF <- min_depth_distribution(origRF) # Obtain the distribution of minimal depth
min_depth_frame_origRF
#plot_min_depth_distribution(origRF, mean_sample = "all_trees") # Same as plotting the min depth frame, but slower
plot_min_depth_distribution(min_depth_frame_origRF, 
                            mean_sample = "all_trees" # available options: "relevant_trees", "all_trees", or "top_trees"
                            #k = 7 # num of features to plot, 10 by default, but we only have 7, so this is optional as well
                            #mean_scale = TRUE # False by default
                            )
# * Min depth in this package seems to mean what is the first time this variable is used to split the tree
# * There is an inverse relationship between times_a_root and mean_min_depth 

```

## 15. Variable Importance Measures
There is more than one way to rank feature importance in a RF. 
Thankfully, RFE provides most, if not all of the feature importance measures available for RF, which makes analyzing the interactions between the different features easier to see and understand. 
15.1 provides a table of the features, but 15.2 allows us better visualize the results by providing a graphical view of the features.   

### 15.1 Variable Importance Table
The *variable importance table* provides seven importance measures by which to rank each of the features used in the RF:
  * mean_min_depth 
  * no_of_nodes
  * accuracy_decrease
  * gini_decrease
  * no_of_trees
  * times_a_root
  * p_value
    * A low p-value generally means that the feature is important for the prediction. 
    * While Distance and DownDistance seem to be the most important variables, SS and Pos10wrtsRNAStart have the lowest p-values.

```{r}
importance_frame_origRF <- measure_importance(origRF)
importance_frame_origRF

# The following lines are simply the importance_frame_origRF organized in a different way.
#importance_frame_origRF[order(-importance_frame_origRF$accuracy_decrease), c("variable","accuracy_decrease")] # equivalent to the graph generated by varImpPlot(origRF)
#importance_frame_origRF[order(-importance_frame_origRF$gini_decrease), c("variable","gini_decrease")] # equivalent to the  graph generated by varImpPlot(origRF)
# importance_frame_origRF[order(importance_frame_origRF$p_value), c("variable","p_value")] # This test tells us whether the observed number of successes (number of nodes in which Xj was used for splitting) exceeds the theoretical number of successes if they were random

```

* sameDownStrand and sameStrand have very little importance in the model
* measure_importance() in the RF Explainer package = importance() and varImpPlot(), already include in the RF package. RF Explainer package just seems to present things prettier

### 15.2 Variable Importance Graphs

```{r}
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
                          )

```

## 16. Measure Comparison using ggpairs 

The diagonal are density plots ("advanced versions of histograms")
```{r}
# plot_importance_ggpairs(origRF) # gives the same result as below but takes longer
importance_frame_origRF[order(-importance_frame_origRF$times_a_root),]

plot_importance_ggpairs(importance_frame_origRF,
                        measures = c("gini_decrease","p_value", "no_of_nodes"))
plot_importance_ggpairs(importance_frame_origRF)

```

```{r}
# plot_importance_rankings(origRF) # gives the same result as below but takes longer
plot_importance_rankings(importance_frame_origRF)
```


```{r}
# (vars <- important_variables(forest, k = 5, measures = c("mean_min_depth", "no_of_trees"))) # gives the same result as below but takes longer
imp_vars <- important_variables(importance_frame_origRF, k = 5, measures = c("mean_min_depth", "gini_decrease"))
imp_vars
#?min_depth_interactions = splits appearing in maximal subtrees with respect to one of the variables selected. 

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

```

