#' @rdname scramblase_assay_plot
#' @importFrom plyr rbind.fill
#' @export
scramblase_assay_stats <- function(
  x,
  scaleTo = c("model","data"),
  forceThroughOrigin = FALSE,
  generation_of_algorithm = c(2, 1),
  splitByExperiment = TRUE){
  UseMethod("scramblase_assay_stats",x)
}
#' @export
scramblase_assay_stats.data.frame <- function(x, ...){
  base_function_scramblase_assay_stats(x, ...)
}
#' @export
scramblase_assay_stats.character <- function(x, ...){
  parsedInputFile <- readScramblaseInputFile(x)
  base_function_scramblase_assay_stats(x=parsedInputFile, ...)
}
base_function_scramblase_assay_stats <- function(
  x,
  scaleTo = c("model","data"),
  forceThroughOrigin = FALSE,
  generation_of_algorithm = c(2, 1),
  splitByExperiment = TRUE){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblaseAssayInputValidation(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin,
    generation_of_algorithm = generation_of_algorithm,
    splitByExperiment = splitByExperiment)
  x <- validatedParams[["x"]]
  scaleTo <- validatedParams[["scaleTo"]]
  forceThroughOrigin <- validatedParams[["forceThroughOrigin"]]
  generation_of_algorithm <- validatedParams[["generation_of_algorithm"]]
  splitByExperiment <- validatedParams[["splitByExperiment"]]

# Processing --------------------------------------------------------------
  processedListFromX <- scramblaseAssayCalculations(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin,
    generation_of_algorithm = generation_of_algorithm,
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
  output["Fit Constant (a)"] <- signif(output["Fit Constant (a)"], digits = 2)
  names(output)[which(names(output) == "Fit Constant (a)")] <- "Fit Constant"
# Return ------------------------------------------------------------------
  return(output)
}