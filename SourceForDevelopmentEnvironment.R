# Load the required libraries
library(devtools)
# # Switch to developer's mode
# devPath <- path.expand("~/R-DevTools/flippant-Development")
# if(!(devPath %in% .libPaths())){
#   dev_mode(path=devPath)
#   update.packages(checkBuilt=TRUE)
#   source("http://bioconductor.org/biocLite.R")
#   biocLite(ask=FALSE)
# }
library(testthat)
# Load the worked on package
library(flippant)
# Use packrat
packrat::on()