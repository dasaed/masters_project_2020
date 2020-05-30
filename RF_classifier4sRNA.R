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

#!/usr/bin/env Rscript

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


# Importing the Model ----

inputFile <- "FeatureTable.tsv" # the sample file
outputFile <- "outputfile.txt"
RF <- readRDS("RF_classifier4sRNA.rds")
data <- read.table( inputFile, sep = "\t", header = TRUE)
pred <- predict(RF, data, type = "prob")
pred
colnames(pred) <- c("No sRNA", "sRNA") # adding header columns to the prediction output
pred
write.table(pred, file = outputFile, sep = "\t", col.names = TRUE, row.names = TRUE)
# Rscript --vanilla RF_classifier4sRNA.R -i FeatureTable.tsv -o outFile.txt


singleDataPoint <- c("sRNA00822","-43.8","-36","1000","-38","1","153","0")
singlePred <- predict(RF, singleDataPoint, type = "prob")


# Applying LIME ----



