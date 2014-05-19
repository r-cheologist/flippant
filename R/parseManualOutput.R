#' @title parseManualOutput
#' @description Parse spectra from files provided in a manually assembled 
#' format.
#' @details A helper function to \code{\link{parseFluorometerOutput}}.
#' The file in question is supposed to contain a tab-delimited table
#' with the columns \code{Time (sec)} and \code{Blank}.
#' @param x \code{\link{character}} vector resulting from 
#' \code{\link{readLines}} of the corresponding data file.
#' @return See \code{\link{parseFluorometerOutput}}.
#' @seealso \code{\link{parseFluorometerOutput}},
#' \code{\link{parseFelixGxOutput}}, \code{\link{parseFelix32Output}}
#' @author Johannes Graumann
#' @keywords manip IO file
parseManualOutput <- function(x=NA){
  #######################
  # Check Prerequisites #
  #######################
  if(is.na(x[1])){
    stop("'x' must be defined.")
  }
  if(!is.character(x)){
    stop("'x' must be of class 'character'.")
  }
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
    "Fluorescense.Intensity"=sapply(output$Data,function(y){y[2]}))
  # Extract/create metadata
  #########################
  # Create data point count
  output$Data.Points <- nrow(output$Data)
  # Create fluorescense max
  output$Max.Fluorescence.Intensity <- max(output$Data$Fluorescense.Intensity,na.rm=TRUE)
  # Create fluorescense min
  output$Min.Fluorescence.Intensity <- min(output$Data$Fluorescense.Intensity,na.rm=TRUE)
  # Create acquisition time  max
  output$Max.Acquisition.Time.in.sec <- max(output$Data$Time.in.sec,na.rm=TRUE)
  # Create acquisition time  min
  output$Min.Acquisition.Time.in.sec <- min(output$Data$Time.in.sec,na.rm=TRUE)
  #################
  # Return result #
  #################
  invisible(output)
}