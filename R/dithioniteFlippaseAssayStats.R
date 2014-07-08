#' @rdname dithioniteFlippaseAssayPlot
#' @importFrom plyr rbind.fill
#' @export
dithioniteFlippaseAssayStats <- function(x,scaleTo=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- dithioniteFlippaseAssayInputValidation(x=x,scaleTo=scaleTo)
  x <- validatedParams[["x"]]
  scaleTo <- validatedParams[["scaleTo"]]

# Processing --------------------------------------------------------------
  processedListFromX <- dithioniteFlippaseAssayCalculations(x=x,scaleTo=scaleTo)
  processedListFromX <- lapply(processedListFromX,function(y){y[["Raw"]]})

# Assemble the output -----------------------------------------
  output <- plyr::rbind.fill(processedListFromX)
  output <- output[!duplicated(output$CombinedId),]
  output <- output[which(names(output) %in% c("Fit Constant (a)","PPR at P = 0.5","Experimental Series","Experiment"))]
  output[c("Fit Constant (a)","PPR at P = 0.5")] <- round(output[c("Fit Constant (a)","PPR at P = 0.5")],digits=2)

# Return ------------------------------------------------------------------
  return(output)
}