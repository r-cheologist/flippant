#' @title ParseQuantMasterData
#' @description Parse spectra from files provided by a Photon QuantMaster 
#' fluorometer
#' @details A function to read fluorometer output directly. Intended as a helper
#' function to flippase activity determinations from dithionite assays.
#' @param SpecFile Path to a \file{*.txt} file as a \code{\link{character}} 
#' object.
#' @return Returns a \code{\link{list}} with the follwoing keys:
#' \describe{
#'  \item{\code{Data}}{A \code{\link{data.frame}} representing the actual 
#'    spectrumwith the columns \code{Time (s)} and \code{Fluorescense Intensity}
#'    (all \code{\link{numeric}}).}
#'  \item{\code{Data Points}}{Number of data points in the spectrum as a 
#'    \code{\link{numeric}}. Ecquivalent to \code{\link{nrow}} of the 
#'    \code{link{data.frame}} in \code{Data}.}
#'  \item{\code{Maximal Fluorescence Intensity}}{\code{\link{numeric}} 
#'    representation of the maximal fluorescence intensity from \code{Data}.}
#'  \item{\code{Minimal Fluorescence Intensity}}{\code{\link{numeric}} 
#'    representation of the minimal fluorescence intensity from \code{Data}.}
#'  \item{\code{File Name}}{\code{\link{character}} representation of the 
#'    file name (as saved by the instrument).}
#'  \item{\code{Maximal Acquisition Time (s)}}{\code{\link{numeric}} 
#'    representation of the maximal \code{Time (s)} from \code{Data}.}
#'  \item{\code{Minimal Acquisition Time (s)}}{\code{\link{numeric}} 
#'    representation of the minimal \code{Time(s)} from \code{Data}.}}
#' @seealso \code{\link{DithioniteFlippaseAssayAnalysis}}
#' @author Johannes Graumann
#' @export
#' @keywords manip IO file
#' @examples
#' stop("Function is missing examples!")
ParseQuantMasterData <- function(SpecFile=NA){
  #######################
  # Check Prerequisites #
  #######################
  if(is.na(SpecFile)){
    stop("'SpecFile' must be defined.")
  }
  if(!is.character(SpecFile)){
    stop("'SpecFile' must be of class 'character'.")
  }
  if(length(SpecFile) != 1){
    stop("'SpecFile' must be of length 1.")
  }
  if(!file.exists(SpecFile)){
    stop("'SpecFile' must represent an existing path.")
  }
  if(file.access(names=SpecFile,mode=4) != 0){
    stop("'SpecFile' must be readable.")
  }
  ##############
  # Processing #
  ##############
  tmpOutput <- list()
  # Aspirate the file
  ###################
  tmpLines <- readLines(SpecFile)
  # Extract data
  ##############
  # Where are data blocks starting?
  blockIndeces <- lapply(
    c("^<Trace>","^X\\tY","^</Trace>"),
    function(x){
      which(sapply(tmpLines,function(y){grep(pattern=x,x=y)}) == 1)
    })
  # Check for data block counts
  for(x in blockIndeces){
    if(length(x) == 0){
      stop("No discernable '",names(x),"' block present - aborting.")
    } else if(length(x) > 1){
      stop("Multiple '",names(x),"' blocks present - aborting.")
    }
  }
  # Ensure order and append line boundaries
  blockIndeces <- unlist(blockIndeces,use.names=TRUE)
  blockIndeces <- c(1,blockIndeces[order(unlist(blockIndeces))]+1,length(tmpLines)+2)
  # Split into blocks
  tmpBlocks <- sapply(
    seq(from=2,to=length(blockIndeces)),
    function(x){
      tmpBlock <- list(tmpLines[(blockIndeces[[x-1]]):(blockIndeces[[x]]-2)])
      names(tmpBlock) <- names(blockIndeces[x-1])
      return(tmpBlock)
    })
  names(tmpBlocks) <- sub(pattern="\\s+$",replacement="",names(tmpBlocks))
  tmpBlocks <- tmpBlocks[c("<Trace>","X\tY")]
  # Rough format check
  elementCountTest <- sapply(tmpBlocks,length) == c(2,NA)
  elementCountTest[is.na(elementCountTest)] <- TRUE
  if(!all(elementCountTest)){
    stop(
      "Data blocks '",
      paste(names(tmpBlocks[!elementCountTest]),collapse=","),
      "' have unexpected member counts.")
  }
  # Process spectral data
  #######################
  # Cull empty entries
  tmpOutput$Data <- tmpBlocks$"X\tY"
  # Parse data apart
  tmpOutput$Data <- lapply(tmpOutput$Data,function(x){unlist(strsplit(x=x,split="\\s+"))})
  # Numerizise
  tmpOutput$Data <- lapply(tmpOutput$Data,as.numeric)
  tmpOutput$Data <- data.frame(
    "Time (s)"=sapply(tmpOutput$Data,function(x){x[1]}),
    "Fluorescense Intensity"=sapply(tmpOutput$Data,function(x){x[2]}),
    check.names=FALSE)
  # Extract/create metadata
  #########################
  # Create data point count
  tmpOutput$"Data Points" <- nrow(tmpOutput$Data)
  # Create fluorescense max
  tmpOutput$"Maximal Fluorescence Intensity" <- max(tmpOutput$Data$"Fluorescense Intensity",na.rm=TRUE)
  # Create fluorescense min
  tmpOutput$"Minimal Fluorescence Intensity" <- min(tmpOutput$Data$"Fluorescense Intensity",na.rm=TRUE)
  # Extract file name
  tmpExtract <- tmpBlocks[["<Trace>"]][2]
  tmpOutput$"File Name" <- tmpExtract
  # Extract/check acquisition time  max
  tmpExtract <- as.numeric(tmpBlocks[["<Trace>"]][1])
  if(tmpExtract != max(tmpOutput$Data$"Time (s)",na.rm=TRUE)){
    stop("Reported and home-made acquisition time maxima differ.")
  }
  tmpOutput$"Maximal Acquisition Time (s)" <- tmpExtract
  # Create acquisition time  min
  tmpOutput$"Minimal Acquisition Time (s)" <- min(tmpOutput$Data$"Time (s)",na.rm=TRUE)
  #################
  # Return result #
  #################
  invisible(tmpOutput)
#   library(ggplot2)
#   tmpPlot <- ggplot(data=tmpOutput$Data,aes_string(x="`Time (s)`",y="`Fluorescense Intensity`"))
#   tmpPlot + geom_line() + geom_point()
}