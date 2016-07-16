#' @title parse_manual_output
#' @description Parse spectra from files provided in a manually assembled 
#' format.
#' @details A helper function to \code{\link{parse_fluorimeter_output}}.
#' The file in question is supposed to contain a tab-delimited table
#' with the columns \code{Time (sec)} and \code{Fluorescense Intensity}.
#' @param x \code{\link{character}} vector resulting from 
#' \code{\link{readLines}} of the corresponding data file.
#' @return See \code{\link{parse_fluorimeter_output}}.
#' @seealso \code{\link{parse_fluorimeter_output}},
#' \code{\link[flippant]{parse_felix_gx_output}}, 
#' \code{\link[flippant]{parse_felix_32_output}}
#' @author Johannes Graumann
#' @keywords manip IO file
parse_manual_output <- function(x = NULL){
  #######################
  # Check Prerequisites #
  #######################
  assertive.properties::assert_is_not_null(x)
  assertive.strings::assert_all_are_non_missing_nor_empty_character(x)
  ##############
  # Processing #
  ##############
  # Extract & process spectral data
  #################################
  output <- list(
    Data = x[2:length(x)])
  # Parse data apart & numerizise
  output$Data <- lapply(
    output$Data,
    function(y){
      as.numeric(unlist(strsplit(x=y,split="\\s+")))
    })
  # Combine into DF & name
  output$Data <- data.frame(
    "Time.in.sec"=sapply(output$Data,function(y){y[1]}),
    "Fluorescence.Intensity"=sapply(output$Data,function(y){y[2]}))
  # Extract/create metadata
  #########################
  # Create data point count
  output$Data.Points <- nrow(output$Data)
  # Create fluorescence max
  output$Max.Fluorescence.Intensity <- max(output$Data$Fluorescence.Intensity,na.rm=TRUE)
  # Create fluorescence min
  output$Min.Fluorescence.Intensity <- min(output$Data$Fluorescence.Intensity,na.rm=TRUE)
  # Create acquisition time  max
  output$Max.Acquisition.Time.in.sec <- max(output$Data$Time.in.sec,na.rm=TRUE)
  # Create acquisition time  min
  output$Min.Acquisition.Time.in.sec <- min(output$Data$Time.in.sec,na.rm=TRUE)
  #################
  # Return result #
  #################
  invisible(output)
}