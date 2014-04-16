#' @rdname dithioniteFlippaseAssayPlot
#' @export
dithioniteFlippaseAssayStats <- function(x,scale_to=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validated_params <- dithioniteFlippaseAssayInputValidation(x=x,scale_to=scale_to)
  x <- validated_params[["x"]]
  scale_to <- validated_params[["scale_to"]]

# Processing --------------------------------------------------------------
  processed_list_from_x <- dithioniteFlippaseAssayCalculations(x=x,scale_to=scale_to)
  processed_list_from_x <- lapply(processed_list_from_x,function(y){y[["Raw"]]})

# Assemble the output -----------------------------------------
  output <- rbind.fill(processed_list_from_x)
  output <- output[!duplicated(output$CombinedId),]
  output <- output[which(names(output) %in% c("Fit Constant (a)","PPR at P = 0.5","Experimental Series","Experiment"))]

# Return ------------------------------------------------------------------
  return(output)
}