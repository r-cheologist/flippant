#' @rdname dithioniteFlippaseAssayPlot
#' @export
dithioniteFlippaseAssayTraces <- function(x,scaleTo=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- dithioniteFlippaseAssayInputValidation(x=x,scaleTo=scaleTo)
  x <- validatedParams[["x"]]
  scaleTo <- validatedParams[["scaleTo"]]

# Processing --------------------------------------------------------------
  # Perform assay calculations to retrive PPR
  processedListFromX <- dithioniteFlippaseAssayCalculations(x=x,scaleTo=scaleTo)
  trimmedProcessedListFromX <- rbind.fill(
      lapply(
        names(processedListFromX),
        function(y){
          processedListFromX[[y]][["Raw"]][c("Path","Experimental Series","Experiment","Protein per Phospholipid (mg/mmol)")]
        }))
  # Parse the fluorometer data and whip it into shape
  rawFlourometerOutput <- lapply(x$Path,parseFluorometerOutput)
  names(rawFlourometerOutput) <- x$Path
  dataFromRawFlourometerOutput <- rbind.fill(
    lapply(
      names(rawFlourometerOutput),
      function(y){data.frame(
        Path=y,
        Time.in.sec=rawFlourometerOutput[[y]][["Data"]][["Time.in.sec"]],
        Fluorescense.Intensity=rawFlourometerOutput[[y]][["Data"]][["Fluorescense.Intensity"]])})
  )
  # Merge spectral data and analysis
  mergedData <- merge(x=dataFromRawFlourometerOutput,y=trimmedProcessedListFromX,by="Path")
  names(mergedData) <- make.names(names(mergedData))
  # Use corresponding PPR as path
  mergedData$Path <- round(mergedData$Protein.per.Phospholipid..mg.mmol.,2)
# Assemble the output -----------------------------------------
  # Groundwork
  plotOutput <- ggplot(
    data=mergedData,
    aes_string(
      x="Time.in.sec",
      y="Fluorescense.Intensity",
      group="Path",
      colour="Path"))
  # Plot traces with lines
  plotOutput <- plotOutput + geom_line()
  # Prettify
  plotOutput <- plotOutput +
    labs(
      x="Time (s)",
      y="Relative Fluorescense Intensity")
  # Facetting
  hasExperiment <- any(!is.na(mergedData$Experiment))
  hasSeries <- any(!is.na(mergedData$Experimental.Series))
  if(hasExperiment && hasSeries){
    plotOutput <- plotOutput + facet_grid(Experimental.Series~Experiment)
  } else if(hasExperiment){
    plotOutput <- plotOutput + facet_wrap(~Experiment)
  } else if(hasSeries){
    plotOutput <- plotOutput + facet_wrap(~Experimental.Series)
  }
  return(plotOutput)
}