#' @rdname scramblaseAssayPlot
#' @importFrom plyr rbind.fill
#' @export
scramblaseAssayStats <- function(x, scaleTo=c("model","data"), forceThroughOrigin=FALSE){
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
baseFunctionScramblaseAssayStats <- function(x,scaleTo=c("model","data"),forceThroughOrigin=FALSE){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblaseAssayInputValidation(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin)
  x <- validatedParams[["x"]]
  scaleTo <- validatedParams[["scaleTo"]]
  forceThroughOrigin <- validatedParams[["forceThroughOrigin"]]

# Processing --------------------------------------------------------------
  processedListFromX <- scramblaseAssayCalculations(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin)
  processedListFromX <- lapply(processedListFromX,function(y){y[["Raw"]]})

# Assemble the output -----------------------------------------
  output <- plyr::rbind.fill(processedListFromX)
  output <- output[!duplicated(output$CombinedId),]
  output <- output[which(names(output) %in% c("Fit Constant (a)","Experimental Series","Experiment"))]
  output["Fit Constant (a)"] <- round(output["Fit Constant (a)"],digits=2)
  names(output)[which(names(output) == "Fit Constant (a)")] <- "Fit Constant"
# Return ------------------------------------------------------------------
  return(output)
}