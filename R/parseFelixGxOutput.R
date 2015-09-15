#' @title parseFelixGxOutput
#' @description Parse spectra from files provided by a Photon QuantMaster 
#' fluorimeter using \code{FelixGX v4.1}
#' @details A helper function to \code{\link{parseFluorimeterOutput}}.
#' @param x \code{\link{character}} vector resulting from 
#' \code{\link{readLines}} of the corresponding data file.
#' @return See \code{\link{parseFluorimeterOutput}}.
#' @seealso \code{\link{parseFluorimeterOutput}}, 
#' \code{\link[flippant]{parseFelix32Output}},
#' \code{\link[flippant]{parseManualOutput}}
#' @importFrom assertive assert_any_are_not_missing_nor_empty_characters
#' @importFrom assertive assert_is_not_null
#' @author Johannes Graumann
#' @keywords manip IO file
parseFelixGxOutput <- function(x=NA){
  #######################
  # Check Prerequisites #
  #######################
  assert_is_not_null(x)
  assert_any_are_not_missing_nor_empty_characters(x)
  ##############
  # Processing #
  ##############
  # Extract data
  ##############
  # Where are data blocks starting?
  blockIndecesInX <- lapply(
    c("^<Trace>","^X\\tY","^</Trace>"),
    function(y){
      which(sapply(x,function(z){grep(pattern=y,x=z)}) == 1)
    })
  # Check for data block counts
  for(y in blockIndecesInX){
    if(length(y) == 0){
      stop("No discernable '",names(y),"' block present - aborting.")
    } else if(length(y) > 1){
      stop("Multiple '",names(y),"' blocks present - aborting.")
    }
  }
  # Ensure order and append line boundaries
  blockIndecesInX <- unlist(blockIndecesInX,use.names=TRUE)
  blockIndecesInX <- c(
    1,
    blockIndecesInX[order(unlist(blockIndecesInX))]+1,
    length(x)+2)
  # Split into blocks
  blocksInX <- sapply(
    seq(from=2,to=length(blockIndecesInX)),
    function(y){
      block <- list(x[(blockIndecesInX[[y-1]]):(blockIndecesInX[[y]]-2)])
      names(block) <- names(blockIndecesInX[y-1])
      return(block)
    })
  names(blocksInX) <- sub(pattern="\\s+$",replacement="",names(blocksInX))
  blocksInX <- blocksInX[c("<Trace>","X\tY")]
  # Rough format check
  rightLengthBlocksInX <- sapply(blocksInX,length) == c(2,NA)
  rightLengthBlocksInX[is.na(rightLengthBlocksInX)] <- TRUE
  if(!all(rightLengthBlocksInX)){
    stop(
      "Data blocks '",
      paste(names(blocksInX[!rightLengthBlocksInX]),collapse=","),
      "' have unexpected member counts.")
  }
  # Process spectral data
  #######################
  output <- list(
    Data = blocksInX$"X\tY")
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
  # Extract file name
  output$File.Name <- blocksInX[["<Trace>"]][2]
  # Extract/check acquisition time  max
  maxAcqTime <- as.numeric(blocksInX[["<Trace>"]][1])
  if(maxAcqTime != max(output$Data$Time.in.sec,na.rm=TRUE)){
    stop("Reported and home-made acquisition time maxima differ.")
  }
  output$Max.Acquisition.Time.in.sec <- maxAcqTime
  # Create acquisition time  min
  output$Min.Acquisition.Time.in.sec <- min(output$Data$Time.in.sec,na.rm=TRUE)
  #################
  # Return result #
  #################
  invisible(output)
}