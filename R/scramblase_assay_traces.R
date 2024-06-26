#' @rdname scramblase_assay_plot
#' @export
scramblase_assay_traces <- function(
  x,
  ppr_scale_factor = 0.65,
  time_min_sec = NA_real_,
  time_max_sec = NA_real_,
  adjust = TRUE,
  timepoint_of_measurement = 400,
  n_averaging = 10,
  annotate_traces = FALSE) {
  UseMethod("scramblase_assay_traces", x)
}
#' @export
scramblase_assay_traces.data.frame <- function(
  x,
  ...) {
  base_function_scramblase_assay_traces(x, ...)
}

#' @export
scramblase_assay_traces.character <- function(x, ...) {
  parsedInputFile <- read_scramblase_input_file(x)
  withr::with_dir(
    dirname(x),
    base_function_scramblase_assay_traces(x = parsedInputFile, ...)
  )
}
base_function_scramblase_assay_traces <- function(
  x,
  ppr_scale_factor = 0.65,
  time_min_sec = NA_real_,
  time_max_sec = NA_real_,
  adjust = TRUE,
  timepoint_of_measurement = 400,
  n_averaging = 10,
  annotate_traces = FALSE) {
# Check Prerequisites -----------------------------------------------------
  protein_is_factor = class(x$`Protein Reconstituted (mg)`) %in%
    c("character", "factor")
  validatedParams <- scramblase_assay_input_validation(
    x = x,
    scale_to = "data",
    ppr_scale_factor = ppr_scale_factor,
    force_through_origin = TRUE,
    generation_of_algorithm = 2,
    split_by_experiment = TRUE,
    r_bar = 88,
    sigma_r_bar = 28,
    verbose = FALSE,
    n_averaging = n_averaging, 
    protein_is_factor = protein_is_factor)
  x <- validatedParams[["x"]]
  assert_is_a_number(time_min_sec)
  assert_is_a_number(time_max_sec)
  assert_is_a_number(timepoint_of_measurement)
  assert_is_a_number(n_averaging)
  assert_all_are_whole_numbers(n_averaging)
  assert_all_are_positive(n_averaging)
  assert_is_a_bool(annotate_traces)

# Processing --------------------------------------------------------------
  columns_retained <- c("Path", "Experimental Series", "Experiment")
  # Perform assay calculations to retrive PPR
  processedX <- x
  if (!protein_is_factor) {
    processedX$"Protein per Phospholipid (mg/mmol)" <- calculate_ppr(
      x, ppr_scale_factor = validatedParams[["ppr_scale_factor"]])
    columns_retained %<>%
      c("Protein per Phospholipid (mg/mmol)")
  } else {
    columns_retained %<>%
      c("Protein Reconstituted (mg)")
  }
  processedX <- processedX %>% magrittr::extract(columns_retained)
  # Parse the fluorimeter data and whip it into shape
  rawFluorimeterOutput <- lapply(
    processedX$Path, parse_fluorimeter_output,
      adjust = adjust, timepoint_of_measurement = timepoint_of_measurement,
      n_averaging = n_averaging)
  names(rawFluorimeterOutput) <- processedX$Path
  dataFromRawFluorimeterOutput <- plyr::rbind.fill(
    lapply(
      names(rawFluorimeterOutput),
      function(y) {
        time.in.sec <- rawFluorimeterOutput[[y]][["Time.in.sec"]]
        fluorescenceIntensity <- rawFluorimeterOutput[[y]][["Fluorescence.Intensity"]]
        fluorescenceIntensity <- fluorescenceIntensity/max(fluorescenceIntensity, na.rm = TRUE)
        return(
          data.frame(
            Path = y,
            Time.in.sec = time.in.sec,
            Fluorescence.Intensity = fluorescenceIntensity))
    }))
  # Deal with fluorescence extrema
  fluorescence_extrema <- rawFluorimeterOutput %>%
    names() %>%
    lapply(function(x) {
      attr(rawFluorimeterOutput[[x]], "FluorescenceExtrema") %>%
        magrittr::divide_by(
          rawFluorimeterOutput[[x]][["Fluorescence.Intensity"]] %>%
            max(na.rm = TRUE)) %>%
        as.list() %>% c(list(Path = x))}) %>%
    lapply(as.data.frame, stringsAsFactors = FALSE) %>%
    plyr::rbind.fill()
  fluorescence_extrema[["Fractional.Fluorescence.Change"]] <- 
                         fluorescence_extrema$Endpoint.Fluorescence %>%
        magrittr::divide_by(fluorescence_extrema$Baseline.Fluorescence)
  # Merge spectral data and analysis
  mergedData <- merge(x = dataFromRawFluorimeterOutput, y = processedX,
                      by = "Path")
  mergedData %<>% merge(y = fluorescence_extrema, by = "Path")
  mergedData$Path %<>% as.character()
  names(mergedData) <- make.names(names(mergedData))
  if (!protein_is_factor) {
    # Use corresponding PPR as a trace identifier
    mergedData$ID <- round(mergedData$Protein.per.Phospholipid..mg.mmol.,2)
  } else {
    mergedData$ID <- mergedData$Protein.Reconstituted..mg.
  }
  
# Assemble the output -----------------------------------------
  # Groundwork
  plotOutput <- ggplot(
    data = mergedData,
    aes_string(
      x      = "Time.in.sec",
      y      = "Fluorescence.Intensity",
      group  = "ID",
      colour = "ID"))
  # Plot traces with lines
  plotOutput <- plotOutput + geom_line()
  # Restrict x axis as requested
  if (any(!is.na(c(time_min_sec,time_max_sec)))) {
    xRange <- range(mergedData$Time.in.sec, na.rm = TRUE)
    if (!is.na(time_min_sec)) {xRange[1] <- time_min_sec}
    if (!is.na(time_max_sec)) {xRange[2] <- time_max_sec}
    plotOutput <- plotOutput + xlim(xRange)
  }
  # Render color scale friendly to the color blind
  # plotOutput <- plotOutput + 
  #   scale_color_continuous(low = "#0072B2",high = "#E69F00")
  # Prettify
  plotOutput <- plotOutput +
    labs(
      x = "Acquisition Time (s)",
      y = "Relative Fluorescence Intensity")
  if (protein_is_factor) {
    plotOutput <- plotOutput + ggplot2::labs(colour = NULL)
  } else {
    if (is.null(ppr_scale_factor)) {
      plotOutput <- plotOutput +
        labs(
          colour = expression("PPR "*bgroup("(",frac("mg","mmol"),")")))
    } else {
      plotOutput <- plotOutput +
        labs(
          colour = expression("adj. PPR "*bgroup("(",frac("mg","mmol"),")")))
    }
  }
  # Annotate as requested
  if (annotate_traces) {
    Fractional.Fluorescence.Change <- Endpoint.Fluorescence <- NULL
    plotOutput <- plotOutput +
      ggplot2::geom_text(
        data = plotOutput$data[!duplicated(plotOutput$data$Path),], 
        mapping = ggplot2::aes(
          label = Fractional.Fluorescence.Change %>%
            magrittr::multiply_by(100) %>%
            round() %>%
            paste0("%"),
          x = Inf, y = Endpoint.Fluorescence),
        position = ggplot2::position_nudge(y = 0.05), vjust = 0, hjust = 1,
        color = "black", size = 3)
  }
  # Facetting
  hasExperiment <- any(!is.na(mergedData$Experiment))
  hasSeries <- any(!is.na(mergedData$Experimental.Series))
  if (hasExperiment && hasSeries) {
    plotOutput <- plotOutput + facet_grid(Experimental.Series~Experiment)
  } else if (hasExperiment) {
    plotOutput <- plotOutput + facet_wrap(~Experiment)
  } else if (hasSeries) {
    plotOutput <- plotOutput + facet_wrap(~Experimental.Series)
  }
  return(plotOutput)
}