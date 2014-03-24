#' @title dithionite_flippase_assay_traces
#' @description A function displaying raw traces for dithionite flippase assays
#' @details Input format for and data processing by the function is documented 
#' at \code{\link{dithionite_flippase_assay_plot}}.
#' @param x \code{\link{data.frame}} as described in 
#' \code{\link{dithionite_flippase_assay_plot}}.
#' @param scale_to Defines the source of \code{ymax}, defaulting to 
#' \code{model}. See \code{\link{dithionite_flippase_assay_plot}}.
#' @return Returns a \code{\link{ggplot}} object.
#' @author Johannes Graumann
#' @export
#' @seealso \code{\link{parse_fluorometer_output}},
#' \code{\link{dithionite_flippase_assay_plot}},
#' \code{\link{dithionite_flippase_assay_input_validation}},
#' \code{\link{dithionite_flippase_assay_calculations}}
#' \code{\link{dithionite_flippase_assay_stats}},
#' @keywords methods manip
#' @examples
#' stop("Add citation to Mike's manuscript!")
#' stop("Revisit workflow description.")
#' stop("Add example using actually published data.")
#' # Build input
#' x <- data.frame(
#'  Path = c(
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-plus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-plus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt"),
#'  "Extract Volume (ul)" = c(0,15,40,75,150,0,15,40,75,150,0,15,40,75,150,0,15,40,75,150),
#'  #     "Reaction Volume w/o DT (ul)" = rep(2000,4),
#'  "Reaction Volume with DT (ul)" = rep(2040,20),
#'  "Concentration Egg PC (mM)" = rep(4.5,20),
#'  "Extract Protein Concentration (mg/ml)" = c(rep(0.67,5),rep(1.26,5),rep(0.67,5),rep(1.26,5)),
#'  #     "Timepoint of Measurement (s)",
#'  "Experimental Series"=c(rep("Extract",5),rep("Depleted Extract",5),rep("Extract",5),rep("Depleted Extract",5)),
#'  Experiment=c(rep("Erg1, Replicate 1",10),rep("Erg1, Replicate 2",10)),
#'  check.names=FALSE,
#'  stringsAsFactors=FALSE)
#'  # Run function
#'  DithioniteFlippaseAssay(x)
dithionite_flippase_assay_traces <- function(x,scale_to=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validated_params <- dithionite_flippase_assay_input_validation(x=x,scale_to=scale_to)
  x <- validated_params[["x"]]
  scale_to <- validated_params[["scale_to"]]

# Processing --------------------------------------------------------------
  # Perform assay calculations to retrive PPR
  processed_list_from_x <- dithionite_flippase_assay_calculations(x=x,scale_to=scale_to)
  trimmed_processed_list_from_x <- rbind.fill(
      lapply(
        names(processed_list_from_x),
        function(y){
          processed_list_from_x[[y]][["Raw"]][c("Path","Experimental Series","Experiment","Protein per Phospholipid (mg/mmol)")]
        }))
  # Parse the fluorometer data and whip it into shape
  raw_florometer_output <- lapply(x$Path,parse_fluorometer_output)
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