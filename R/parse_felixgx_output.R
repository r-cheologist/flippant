#' @title parse_felixgx_output
#' @description Parse spectra from files provided by a Photon QuantMaster 
#' fluorometer using \code{FelixGX v4.1}
#' @details A helper function to \code{\link{parse_fluorometer_output}}.
#' @param x \code{\link{character}} vector resulting from 
#' \code{\link{readLines}} of the corresponding data file.
#' @return See \code{\link{parse_fluorometer_output}}.
#' @seealso \code{\link{parse_fluorometer_output}}, 
#' \code{\link{parse_felix32_output}},
#' \code{\link{parse_manual_output}}
#' @author Johannes Graumann
#' @keywords manip IO file
parse_felixgx_output <- function(x=NA){
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
  # Extract data
  ##############
  # Where are data blocks starting?
  block_indeces_in_x <- lapply(
    c("^<Trace>","^X\\tY","^</Trace>"),
    function(y){
      which(sapply(x,function(z){grep(pattern=y,x=z)}) == 1)
    })
  # Check for data block counts
  for(y in block_indeces_in_x){
    if(length(y) == 0){
      stop("No discernable '",names(y),"' block present - aborting.")
    } else if(length(y) > 1){
      stop("Multiple '",names(y),"' blocks present - aborting.")
    }
  }
  # Ensure order and append line boundaries
  block_indeces_in_x <- unlist(block_indeces_in_x,use.names=TRUE)
  block_indeces_in_x <- c(
    1,
    block_indeces_in_x[order(unlist(block_indeces_in_x))]+1,
    length(x)+2)
  # Split into blocks
  blocks_in_x <- sapply(
    seq(from=2,to=length(block_indeces_in_x)),
    function(y){
      block <- list(x[(block_indeces_in_x[[y-1]]):(block_indeces_in_x[[y]]-2)])
      names(block) <- names(block_indeces_in_x[y-1])
      return(block)
    })
  names(blocks_in_x) <- sub(pattern="\\s+$",replacement="",names(blocks_in_x))
  blocks_in_x <- blocks_in_x[c("<Trace>","X\tY")]
  # Rough format check
  right_length_blocks_in_x <- sapply(blocks_in_x,length) == c(2,NA)
  right_length_blocks_in_x[is.na(right_length_blocks_in_x)] <- TRUE
  if(!all(right_length_blocks_in_x)){
    stop(
      "Data blocks '",
      paste(names(blocks_in_x[!right_length_blocks_in_x]),collapse=","),
      "' have unexpected member counts.")
  }
  # Process spectral data
  #######################
  output <- list(
    Data = blocks_in_x$"X\tY")
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
  # Extract file name
  output$File.Name <- blocks_in_x[["<Trace>"]][2]
  # Extract/check acquisition time  max
  max_acq_time <- as.numeric(blocks_in_x[["<Trace>"]][1])
  if(max_acq_time != max(output$Data$Time.in.sec,na.rm=TRUE)){
    stop("Reported and home-made acquisition time maxima differ.")
  }
  output$Max.Acquisition.Time.in.sec <- max_acq_time
  # Create acquisition time  min
  output$Min.Acquisition.Time.in.sec <- min(output$Data$Time.in.sec,na.rm=TRUE)
  #################
  # Return result #
  #################
  invisible(output)
}