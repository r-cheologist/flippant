# Load the required libraries
library(devtools)
# Switch to developer's mode
devPath <- path.expand("~/R-DevTools/mdaa-Development")
if(!(devPath %in% .libPaths())){
  dev_mode(path=devPath)
  update.packages(checkBuilt=TRUE)
}
library(testthat)
# Ensure correct RCFPD package
library(RCFPD,lib.loc="~/R-DevTools/mdaa-Development")
# Load the worked on package
library(MDAA)