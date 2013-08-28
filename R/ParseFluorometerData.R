#' @title ParseFluorometerData
#' @description Parse spectra from files provided by a Perkin Elmer FS55
#' @details A function to read fluorometer output directly. Intended as a helper
#' function to flippase activity determinations from dithionite assays.
#' @param SpecFile Path to a \file{*.td} file as a \code{\link{character}} 
#' object.
#' @return Returns a \code{\link{list}} with the follwoing keys:
#' \describe{
#'  \item{\code{Data}}{A \code{\link{data.frame}} representing the actual 
#'    spectrumwith the columns \code{Time (s)} and \code{Fluorescense Intensity}
#'    (all \code{\link{numerical}}).}
#'  \item{\code{Data Points}}{Number of data points in the spectrum as a 
#'    \code{\link{numerical}}. Ecquivalent to \code{\link{nrow}} of the 
#'    \code{link{data.frame}} in \code{Data}.}
#'  \item{\code{Maximal Flurescence Intensity}}{\code{\link{numerical}} 
#'    representation of the maximal fluorescence intensity from \code{Data}.}
#'  \item{\code{Minimal Flurescence Intensity}}{\code{\link{numerical}} 
#'    representation of the minimal fluorescence intensity from \code{Data}.}
#'  \item{\code{File Name}}{\code{\link{character}} representation of the 
#'    original file name used by the instrument. Often useless as MS 
#'    Windows-specific truncation applies.}
#'  \item{\code{Acquisition Date (YY/MM/DD)}}{\code{\link{character}} 
#'    representation of the date of spectrum acquisition.}
#'  \item{\code{Acquisition Time (HH:MM:SS)}}{\code{\link{character}} 
#'    representation of the time of spectrum acquisition.}
#'  \item{\code{Access Date (YY/MM/DD)}}{\code{\link{character}} 
#'    representation of the date of file access.}
#'  \item{\code{Access Time (HH:MM:SS)}}{\code{\link{character}} 
#'    representation of the time of file access.}
#'  \item{\code{Maximal Acquisition Time (s)}}{\code{\link{numerical}} 
#'    representation of the maximal \code{Time (s)} from \code{Data}.}
#'  \item{\code{Minimal Acquisition Time (s)}}{\code{\link{numerical}} 
#'    representation of the minimal \code{Time(s)} from \code{Data}.}
#'  \item{\code{Instrument}}{\code{\link{character}} representation of the 
#'    spectrometer used.}
#'  \item{\code{Method File}}{\code{\link{character}} representation of the 
#'    file name of the method used by the instrument. Often useless as MS 
#'    Windows-specific truncation applies.}
#'  \item{\code{Excitation Wavelength (nm)}}{\code{\link{numerical}} 
#'    representation of the excitation wavelength used.}
#'  \item{\code{Emitting Wavelength (nm)}}{\code{\link{numerical}} 
#'    representation of the emission wavelength recorded.}}
#' @author Johannes Graumann
#' @export
#' @keywords manip IO file
#' @examples
#' stop("Function is missing examples!")
ParseFluorometerData <- function(SpecFile=NA){
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
    c("^\\#HDR","^\\#GR","^\\#DATA"),
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
  names(tmpBlocks)[names(tmpBlocks) == ""] <- "Preamble"
  # Rough format check
  elementCountTest <- sapply(tmpBlocks,length) == c(39,2,10,NA)
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
  tmpOutput$Data <- tmpBlocks[["#DATA"]][nchar(tmpBlocks[["#DATA"]]) > 0]
  # Parse data apart
  tmpOutput$Data <- lapply(tmpOutput$Data,function(x){unlist(strsplit(x=x,split="\\s+"))})
  # Numerizise
  tmpOutput$Data <- lapply(tmpOutput$Data,as.numeric)
  tmpOutput$Data <- data.frame(
    "Time (s)"=sapply(tmpOutput$Data,function(x){x[1]}),
    "Fluorescense Intensity"=sapply(tmpOutput$Data,function(x){x[2]}),
    check.names=FALSE)
  # Extract metadata
  ##################
  # Extract/check data point count
  tmpExtract <- as.numeric(tmpBlocks[["#GR"]][7])
  if(tmpExtract != nrow(tmpOutput$Data)){
    stop("Data points of home-extracted spectrum differ from in-file report.")
  }
  tmpOutput$"Data Points" <- tmpExtract
  # Extract/check fluorescense max
  tmpExtract <- as.numeric(tmpBlocks[["#GR"]][9])
  if(tmpExtract != max(tmpOutput$Data$"Fluorescense Intensity",na.rm=TRUE)){
    stop("Reported and home-made fluorescense intentity maxima differ.")
  }
  tmpOutput$"Maximal Flurescence Intensity" <- tmpExtract
  # Extract/check fluorescense min
  tmpExtract <- as.numeric(tmpBlocks[["#GR"]][10])
  if(tmpExtract != min(tmpOutput$Data$"Fluorescense Intensity",na.rm=TRUE)){
    stop("Reported and home-made fluorescense intentity minima differ.")
  }
  tmpOutput$"Minimal Flurescence Intensity" <- tmpExtract
  # Extract file name
  tmpExtract <- tmpBlocks[["Preamble"]][3]
  tmpOutput$"File Name" <- tmpExtract
  # Extract acquisition date/time
  tmpExtract <- tmpBlocks[["Preamble"]][4]
  tmpOutput$"Acquisition Date (YY/MM/DD)" <- tmpExtract
  tmpExtract <- tmpBlocks[["Preamble"]][5]
  tmpOutput$"Acquisition Time (HH:MM:SS)" <- tmpExtract
  # Extract access date/time
  tmpExtract <- tmpBlocks[["Preamble"]][6]
  tmpOutput$"Access Date (YY/MM/DD)" <- tmpExtract
  tmpExtract <- tmpBlocks[["Preamble"]][7]
  tmpOutput$"Access Time (HH:MM:SS)" <- tmpExtract
  # Extract/check acquisition time  max
  tmpExtract <- as.numeric(tmpBlocks[["Preamble"]][10])
  if(tmpExtract != max(tmpOutput$Data$"Time (s)",na.rm=TRUE)){
    stop("Reported and home-made acquisition time maxima differ.")
  }
  tmpOutput$"Maximal Acquisition Time (s)" <- tmpExtract
  # Extract/check acquisition time  min
  tmpExtract <- as.numeric(tmpBlocks[["Preamble"]][11])
  if(tmpExtract != min(tmpOutput$Data$"Time (s)",na.rm=TRUE)){
    stop("Reported and home-made acquisition time minima differ.")
  }
  tmpOutput$"Minimal Acquisition Time (s)" <- tmpExtract
  # Extract instrument type
  tmpExtract <- tmpBlocks[["Preamble"]][21]
  tmpOutput$"Instrument" <- tmpExtract
  # Extract instrument method
  tmpExtract <- tmpBlocks[["Preamble"]][25]
  tmpOutput$"Method File" <- tmpExtract
  # Extract excitation/emitting wave length
  tmpExtract <- tmpBlocks[["Preamble"]][38]
  tmpExtract <- as.numeric(unlist(strsplit(x=tmpExtract,split="\\s+")))
  tmpOutput$"Excitation Wavelength (nm)" <- tmpExtract[1]
  tmpOutput$"Emitting Wavelength (nm)" <- tmpExtract[2]
  #################
  # Return result #
  #################
  invisible(tmpOutput)
#   library(ggplot2)
#   tmpPlot <- ggplot(data=tmpOutput$Data,aes_string(x="`Time (s)`",y="`Fluorescense Intensity`"))
#   tmpPlot + geom_line() + geom_point()
}