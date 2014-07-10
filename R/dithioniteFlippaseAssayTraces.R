#' @rdname dithioniteFlippaseAssayPlot
#' @import ggplot2
#' @importFrom assertive assert_is_a_number
#' @importFrom plyr rbind.fill
#' @export
dithioniteFlippaseAssayTraces <- function(x,timeMax=NA_real_){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- dithioniteFlippaseAssayInputValidation(x=x,scaleTo="data",forceThroughOrigin=TRUE)
  x <- validatedParams[["x"]]
  assert_is_a_number(timeMax)

# Processing --------------------------------------------------------------
  # Perform assay calculations to retrive PPR
  processedX <- x
  processedX$"Protein per Phospholipid (mg/mmol)" <- calculatePpr(x)
  processedX <- processedX[c("Path","Experimental Series","Experiment","Protein per Phospholipid (mg/mmol)")]
  # Parse the fluorometer data and whip it into shape
  rawFlourometerOutput <- lapply(processedX$Path,parseFluorometerOutput)
  names(rawFlourometerOutput) <- processedX$Path
  dataFromRawFlourometerOutput <- plyr::rbind.fill(
    lapply(
      names(rawFlourometerOutput),
      function(y){
        time.in.sec <- rawFlourometerOutput[[y]][["Data"]][["Time.in.sec"]]
        fluorescenseIntensity <- rawFlourometerOutput[[y]][["Data"]][["Fluorescense.Intensity"]]
        fluorescenseIntensity <- fluorescenseIntensity/max(fluorescenseIntensity,na.rm=TRUE)
        return(
          data.frame(
            Path = y,
            Time.in.sec = time.in.sec,
            Fluorescense.Intensity=fluorescenseIntensity))
    }))
  # Merge spectral data and analysis
  mergedData <- merge(x=dataFromRawFlourometerOutput,y=processedX,by="Path")
  names(mergedData) <- make.names(names(mergedData))
  timeMin <- min(mergedData$Time.in.sec, na.rm=TRUE)
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
  # Restrict x axis as requested
  if(!is.na(timeMax)){
    plotOutput <- plotOutput + xlim(timeMin,timeMax)
  }
  # Prettify
  plotOutput <- plotOutput +
    labs(
      x="Time (s)",
      y="Relative Fluorescense Intensity",
      colour=expression("PPR "*bgroup("(",frac("mg","mmol"),")")))
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