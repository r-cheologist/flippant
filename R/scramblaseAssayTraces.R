#' @rdname scramblaseAssayPlot
#' @import ggplot2
#' @importFrom assertive assert_is_a_number
#' @importFrom plyr rbind.fill
#' @export
scramblaseAssayTraces <- function(x, timeMin=NA_real_, timeMax=NA_real_){
  UseMethod("scramblaseAssayTraces",x)
}
#' @export
scramblaseAssayTraces.data.frame <- function(x, ...){
  baseFunctionScramblaseAssayTraces(x, ...)
}
#' @export
scramblaseAssayTraces.character <- function(x, ...){
  parsedInputFile <- readScramblaseInputFile(x)
  baseFunctionScramblaseAssayTraces(x=parsedInputFile, ...)
}
baseFunctionScramblaseAssayTraces <- function(x,timeMin=NA_real_,timeMax=NA_real_){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblaseAssayInputValidation(
    x = x,
    scaleTo = "data",
    forceThroughOrigin = TRUE,
    verbose = FALSE)
  x <- validatedParams[["x"]]
  assert_is_a_number(timeMin)
  assert_is_a_number(timeMax)

# Processing --------------------------------------------------------------
  # Perform assay calculations to retrive PPR
  processedX <- x
  processedX$"Protein per Phospholipid (mg/mmol)" <- calculatePpr(x)
  processedX <- processedX[c("Path","Experimental Series","Experiment","Protein per Phospholipid (mg/mmol)")]
  # Parse the fluorimeter data and whip it into shape
  rawFluorimeterOutput <- lapply(processedX$Path,parseFluorimeterOutput)
  names(rawFluorimeterOutput) <- processedX$Path
  dataFromRawFluorimeterOutput <- plyr::rbind.fill(
    lapply(
      names(rawFluorimeterOutput),
      function(y){
        time.in.sec <- rawFluorimeterOutput[[y]][["Data"]][["Time.in.sec"]]
        fluorescenceIntensity <- rawFluorimeterOutput[[y]][["Data"]][["Fluorescence.Intensity"]]
        fluorescenceIntensity <- fluorescenceIntensity/max(fluorescenceIntensity,na.rm=TRUE)
        return(
          data.frame(
            Path = y,
            Time.in.sec = time.in.sec,
            Fluorescence.Intensity=fluorescenceIntensity))
    }))
  # Merge spectral data and analysis
  mergedData <- merge(x=dataFromRawFluorimeterOutput,y=processedX,by="Path")
  names(mergedData) <- make.names(names(mergedData))
  # Use corresponding PPR as path
  mergedData$Path <- round(mergedData$Protein.per.Phospholipid..mg.mmol.,2)
# Assemble the output -----------------------------------------
  # Groundwork
  plotOutput <- ggplot(
    data=mergedData,
    aes_string(
      x="Time.in.sec",
      y="Fluorescence.Intensity",
      group="Path",
      colour="Path"))
  # Plot traces with lines
  plotOutput <- plotOutput + geom_line()
  # Restrict x axis as requested
  if(any(!is.na(c(timeMin,timeMax)))){
    xRange <- range(mergedData$Time.in.sec, na.rm=TRUE)
    if(!is.na(timeMin)){xRange[1] <- timeMin}
    if(!is.na(timeMax)){xRange[2] <- timeMax}
    plotOutput <- plotOutput + xlim(xRange)
  }
  # Prettify
  plotOutput <- plotOutput +
    labs(
      x="Time (s)",
      y="Relative Fluorescence Intensity",
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