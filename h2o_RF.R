# http://uc-r.github.io/lime
# required packages
# install vip from github repo: devtools::install_github("koalaverse/vip")
library(lime)       # ML local interpretation
install.packages("vip")
library(vip)        # ML global interpretation
install.packages("pdp")
library(pdp)        # ML global interpretation
library(ggplot2)    # visualization pkg leveraged by above packages
library(caret)      # ML model building
library(h2o)        # ML model building
h2o.init()
h2o.no_progress()
??rsample
