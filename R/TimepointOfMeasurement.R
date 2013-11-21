#' @export
TimepointOfMeasurement <- function(x=NA,Fluorometer=c("MasterQuant","LS55")){
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
  Fluorometer <- match.arg(arg=Fluorometer,choices=c("QuantMaster","LS55"),several.ok=FALSE)
  ##############
  # Processing #
  ##############
  # Parsing
  #########
  if(Fluorometer=="QuantMaster"){
    tmpData <- lapply(
      x,
      function(y){ParseQuantMasterData(SpecFile=y)})
  } else if(Fluorometer=="LS55"){
    tmpData <- lapply(
      x,
      function(y){ParseLS55Data(SpecFile=y)})
  }
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