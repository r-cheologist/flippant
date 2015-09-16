#' @title parseManualOutput
#' @description Parse spectra from files provided in a manually assembled 
#' format.
#' @details A helper function to \code{\link{parseFluorimeterOutput}}.
#' The file in question is supposed to contain a tab-delimited table
#' with the columns \code{Time (sec)} and \code{Fluorescense Intensity}.
#' @param x \code{\link{character}} vector resulting from 
#' \code{\link{readLines}} of the corresponding data file.
#' @return See \code{\link{parseFluorimeterOutput}}.
#' @seealso \code{\link{parseFluorimeterOutput}},
#' \code{\link[flippant]{parseFelixGxOutput}}, 
#' \code{\link[flippant]{parseFelix32Output}}
#' @importFrom assertive assert_any_are_not_missing_nor_empty_characters
#' @importFrom assertive assert_is_not_null
#' @author Johannes Graumann
#' @keywords manip IO file
parseManualOutput <- function(x = NULL){
  #######################
  # Check Prerequisites #
  #######################
  assert_is_not_null(x)
  assert_any_are_not_missing_nor_empty_characters(x)
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