timepoint_of_measurement <- function(x=NA){
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
  tmp_data <- lapply(
    x,
    parse_fluorometer_output)
  # Extract data
  ##############
  # What is the latest data point per series?
  maxAT <- vapply(tmp_data,function(y){y$Max.Acquisition.Time.in.sec},1)
  # Identify the shortest one
  minMaxAT <- min(maxAT,na.rm=TRUE)
  ##########
  # Return #
  ##########
  return(minMaxAT)
}