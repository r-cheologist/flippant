#' @rdname dithioniteFlippaseAssayPlot
#' @export
dithioniteFlippaseAssayTraces <- function(x,scale_to=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validated_params <- dithioniteFlippaseAssayInputValidation(x=x,scale_to=scale_to)
  x <- validated_params[["x"]]
  scale_to <- validated_params[["scale_to"]]

# Processing --------------------------------------------------------------
  # Perform assay calculations to retrive PPR
  processed_list_from_x <- dithioniteFlippaseAssayCalculations(x=x,scale_to=scale_to)
  trimmed_processed_list_from_x <- rbind.fill(
      lapply(
        names(processed_list_from_x),
        function(y){
          processed_list_from_x[[y]][["Raw"]][c("Path","Experimental Series","Experiment","Protein per Phospholipid (mg/mmol)")]
        }))
  # Parse the fluorometer data and whip it into shape
  raw_florometer_output <- lapply(x$Path,parseFluorometerOutput)
  names(raw_florometer_output) <- x$Path
  data_from_raw_florometer_output <- rbind.fill(
    lapply(
      names(raw_florometer_output),
      function(y){data.frame(
        Path=y,
        Time.in.sec=raw_florometer_output[[y]][["Data"]][["Time.in.sec"]],
        Fluorescense.Intensity=raw_florometer_output[[y]][["Data"]][["Fluorescense.Intensity"]])})
  )
  # Merge spectral data and analysis
  merged_data <- merge(x=data_from_raw_florometer_output,y=trimmed_processed_list_from_x,by="Path")
  names(merged_data) <- make.names(names(merged_data))
  # Use corresponding PPR as path
  merged_data$Path <- round(merged_data$Protein.per.Phospholipid..mg.mmol.,2)
# Assemble the output -----------------------------------------
  # Groundwork
  plot <- ggplot(
    data=merged_data,
    aes_string(
      x="Time.in.sec",
      y="Fluorescense.Intensity",
      group="Path",
      colour="Path"))
  # Plot traces with lines
  plot <- plot + geom_line()
  # Prettify
  plot <- plot +
    labs(
      x="Time (s)",
      y="Relative Fluorescense Intensity")
  # Facetting
  hasExperiment <- any(!is.na(merged_data$Experiment))
  hasSeries <- any(!is.na(merged_data$Experimental.Series))
  if(hasExperiment && hasSeries){
    plot <- plot + facet_grid(Experimental.Series~Experiment)
  } else if(hasExperiment){
    plot <- plot + facet_wrap(~Experiment)
  } else if(hasSeries){
    plot <- plot + facet_wrap(~Experimental.Series)
  }
  return(plot)
}