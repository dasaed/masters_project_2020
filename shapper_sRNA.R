
# A SHAP implementation exist for R, but in reality, it uses the python library. As a presetup for this script to work, you must run 'pip3 install shap' or 'pip install shap', to first fetch the python library.

#install.packages("shapper")
#install.packages("randomForest")

library("shapper")
library("randomForest")
library(data.table)

# Load the data ----
trainDataSet <- read.csv("./DataSets/combinedData.csv", header = TRUE)
trainDataSet[,"Class"] <- as.logical(trainDataSet[,"Class"]) 

trainDataSet_norm <- read.csv("./DataSets/combinedData_norm.csv", header = TRUE)
trainDataSet_norm[,"Class"] <- as.logical(trainDataSet_norm[,"Class"]) 

slt2_sampleData_norm <- read.csv("./DataSets/slt2_data_to_analyze_normalized.csv", header = TRUE, row.names = 1)
lu_sampleData_norm <- read.csv("./DataSets/lu_data_to_analyze_normalized.csv", header = TRUE, row.names = 1)
sampleData <- rbind(slt2_sampleData_norm, lu_sampleData_norm) # This is the data that will ultimately be used by the IM

# Load the RF models ----
origRF <- readRDS("./RFModels/RF_classifier4sRNA.rds")
origRF_norm <- readRDS("./RFModels/OrigRF_norm.rds")

# Create table to save all the results ----
shap_values_table <- data.frame("SS" = NULL, "Pos10wrtsRNAStart"= NULL, "DistTerm"= NULL, "Distance"= NULL, "sameStrand"= NULL, "DownDistance"= NULL, "sameDownStrand"= NULL, "Prediction"= NULL)

# Create Prediction function for shapper library ----
p_function <- function(model, data) predict(model, newdata = data, type = "prob")

# Run shap against each instance 4 times ----
for (instance in 1:nrow(sampleData) ) {
    instance = 1
    instance_name <-rownames(sampleData[instance,])
    ive_rf_a <- individual_variable_effect(origRF_norm, data = trainDataSet_norm[,c(1:7)], predict_function = p_function, new_observation = sampleData[instance,-8])
    ive_rf_b <- individual_variable_effect(origRF_norm, data = trainDataSet_norm[,c(1:7)], predict_function = p_function, new_observation = sampleData[instance,-8])
    ive_rf_c <- individual_variable_effect(origRF_norm, data = trainDataSet_norm[,c(1:7)], predict_function = p_function, new_observation = sampleData[instance,-8])
    ive_rf_d <- individual_variable_effect(origRF_norm, data = trainDataSet_norm[,c(1:7)], predict_function = p_function, new_observation = sampleData[instance,-8])
    ive_rf_b
    if ( all.equal(ive_rf_a,ive_rf_b) & all.equal(ive_rf_a,ive_rf_c) & all.equal(ive_rf_a,ive_rf_d) ){
      ive_rf <- ive_rf_a
    }
    
    shap_values_summary <- ive_rf[ive_rf$`_ylevel_` == 1,c("_vname_","_attribution_")]
    shap_values_temp_table <- transpose(shap_values_summary)
    shap_values_temp_table_with_pred <- cbind( shap_values_temp_table[2,] , predict(origRF, sampleData[instance,-8], type = "prob")[1,"1"] )
    colnames(shap_values_temp_table_with_pred) <- c("SS", "Pos10wrtsRNAStart", "DistTerm", "Distance", "sameStrand", "DownDistance", "sameDownStrand", "Prediction")
    shap_values_table <- rbind(shap_values_table, shap_values_temp_table_with_pred)
    rownames(shap_values_table)[nrow(shap_values_table)] <- instance_name
}

# Print and save the results ----
shap_values_table
write.csv(shap_values_table, "./Results/SHAP_explanations.csv", row.names = TRUE)

# plot(ive_rf[ive_rf$`_ylevel_` == 1,], bar_width = 10 )
