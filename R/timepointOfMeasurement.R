#' @title timepointOfMeasurement
#' @description Helper function to establish a measurement time window for which
#' fluorescence intensities are averaged in a dithionite flippase assay.
#' @param x A \code{\link{character}} vector referring to fluorometer output
#' parseable by \code{\link{parseFluorometerOutput}}.
#' @details The function identifies the maximum reported measurement time from
#' all spectral series involved and reports the smallest of them.
#' @return A \code{\link{numeric}} of length \code{1}.
#' @author Johannes Graumann
#' @seealso \code{\link{parseFluorometerOutput}}
#' @keywords manip
#' @examples
#' stop("Function is missing examples!")
timepointOfMeasurement <- function(x=NA){
  #######################
  # Check prerequisites #
  #######################
  if(!is.character(x)){
    stop("'x' must be of class 'character'.")
  }
  if(!all(file.exists(x))){
    stop("All elements in 'x' must refer to existing files.")
  }
  if(any(file.access(x,mode=4) == -1)){
    stop("All elements in 'x' must refer to readable files.")
  }
  ##############
  # Processing #
  ##############
  # Parsing
  #########
  tmpData <- lapply(
    x,
    parseFluorometerOutput)
  # Extract data
  ##############
  # What is the latest data point per series?
  maxAt <- vapply(tmpData,function(y){y$Max.Acquisition.Time.in.sec},1)
  # Identify the shortest one
  minMaxAt <- min(maxAt,na.rm=TRUE)
  ##########
  # Return #
  ##########
  return(minMaxAt)
}