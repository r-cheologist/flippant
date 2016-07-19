#' @rdname scramblase_assay_plot
#' @export
scramblase_assay_traces <- function(
  x,
  ppr_scale_factor = 0.65,
  time_min_sec=NA_real_,
  time_max_sec=NA_real_,
  adjust = TRUE){
  UseMethod("scramblase_assay_traces", x)
}
#' @export
scramblase_assay_traces.data.frame <- function(
  x,
  ...){
  base_function_scramblase_assay_traces(x, ...)
}

#' @export
scramblase_assay_traces.character <- function(x, ...){
  parsedInputFile <- read_scramblase_input_file(x)
  withr::with_dir(
    dirname(x),
    base_function_scramblase_assay_traces(x=parsedInputFile, ...)
  )
}
base_function_scramblase_assay_traces <- function(
  x,
  ppr_scale_factor = 0.65,
  time_min_sec=NA_real_,
  time_max_sec=NA_real_,
  adjust = TRUE){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblase_assay_input_validation(
    x = x,
    scale_to = "data",
    ppr_scale_factor = ppr_scale_factor,
    force_through_origin = TRUE,
    generation_of_algorithm = 2,
    split_by_experiment = TRUE,
    r_bar = 88,
    sigma_r_bar = 28,
    verbose = FALSE)
  x <- validatedParams[["x"]]
  assertive.types::assert_is_a_number(time_min_sec)
  assertive.types::assert_is_a_number(time_max_sec)

# Processing --------------------------------------------------------------
  # Perform assay calculations to retrive PPR
  processedX <- x
  processedX$"Protein per Phospholipid (mg/mmol)" <- calculate_ppr(x, ppr_scale_factor = validatedParams[["ppr_scale_factor"]])
  processedX <- processedX[c("Path","Experimental Series","Experiment","Protein per Phospholipid (mg/mmol)")]
  # Parse the fluorimeter data and whip it into shape
  rawFluorimeterOutput <- lapply(processedX$Path,parse_fluorimeter_output,adjust = adjust)
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
  if(any(!is.na(c(time_min_sec,time_max_sec)))){
    xRange <- range(mergedData$Time.in.sec, na.rm=TRUE)
    if(!is.na(time_min_sec)){xRange[1] <- time_min_sec}
    if(!is.na(time_max_sec)){xRange[2] <- time_max_sec}
    plotOutput <- plotOutput + xlim(xRange)
  }
  # Render color scale friendly to the color blind
  plotOutput <- plotOutput + 
    scale_color_continuous(low="#0072B2",high="#E69F00")
  # Prettify
  plotOutput <- plotOutput +
    labs(
      x="Acquisition Time (s)",
      y="Relative Fluorescence Intensity")
  if(is.null(ppr_scale_factor)){
    plotOutput <- plotOutput +
      labs(
        colour=expression("PPR "*bgroup("(",frac("mg","mmol"),")")))
  } else {
    plotOutput <- plotOutput +
      labs(
        colour=expression("adj. PPR "*bgroup("(",frac("mg","mmol"),")")))
  }
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