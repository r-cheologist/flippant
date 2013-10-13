#' @export
TimepointOfMeasurement <- function(x=NA){
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
    function(y){ParseQuantMasterData(SpecFile=y)})
  # Extract data
  ##############
  # What is the latest data point per series?
  maxAT <- vapply(tmpData,function(y){y$"Maximal Acquisition Time (s)"},1)
  # Identify the shortest one
  minMaxAT <- min(maxAT,na.rm=TRUE)
  ##########
  # Return #
  ##########
  return(minMaxAT)
}