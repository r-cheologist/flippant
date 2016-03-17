#' @rdname scramblaseAssayPlot
#' @importFrom plyr rbind.fill
#' @export
scramblaseAssayStats <- function(
  x,
  scaleTo = c("model","data"),
  forceThroughOrigin = FALSE,
  formulaGeneration = c(2, 1),
  splitByExperiment = TRUE){
  UseMethod("scramblaseAssayStats",x)
}
#' @export
scramblaseAssayStats.data.frame <- function(x, ...){
  baseFunctionScramblaseAssayStats(x, ...)
}
#' @export
scramblaseAssayStats.character <- function(x, ...){
  parsedInputFile <- readScramblaseInputFile(x)
  baseFunctionScramblaseAssayStats(x=parsedInputFile, ...)
}
baseFunctionScramblaseAssayStats <- function(
  x,
  scaleTo = c("model","data"),
  forceThroughOrigin = FALSE,
  splitByExperiment = TRUE){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblaseAssayInputValidation(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin,
    formulaGeneration = formulaGeneration,
    splitByExperiment = splitByExperiment)
  x <- validatedParams[["x"]]
  scaleTo <- validatedParams[["scaleTo"]]
  forceThroughOrigin <- validatedParams[["forceThroughOrigin"]]
  formulaGeneration <- validatedParams[["formulaGeneration"]]
  splitByExperiment <- validatedParams[["splitByExperiment"]]

# Processing --------------------------------------------------------------
  processedListFromX <- scramblaseAssayCalculations(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin,
    splitByExperiment = splitByExperiment)
  processedListFromX <- lapply(processedListFromX,function(y){y[["Raw"]]})

# Assemble the output -----------------------------------------
  output <- plyr::rbind.fill(processedListFromX)
  if(splitByExperiment){
   output <- output[!duplicated(output$CombinedId),]
  } else {
    output <- output[!duplicated(output$`Experimental Series`),]
  }
  if(splitByExperiment){
    columns <- c("Fit Constant (a)","Experimental Series","Experiment")
  } else {
    columns <- c("Fit Constant (a)","Experimental Series")
  }
  output <- output[which(names(output) %in% columns)]
  output["Fit Constant (a)"] <- round(output["Fit Constant (a)"],digits=2)
  names(output)[which(names(output) == "Fit Constant (a)")] <- "Fit Constant"
# Return ------------------------------------------------------------------
  return(output)
}