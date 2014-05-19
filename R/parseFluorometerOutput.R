#' @title parseFluorometerOutput
#' @description Parse fluorometer spectra
#' @details A function to read fluorometer output directly. Intended as a helper
#' function to flippase activity determinations from dithionite assays.
#' 
#' The function is currently capable to deal with input derived from QuantMaster
#' instruments running software versions \code{FelixGX v4.1} 
#' (see \code{\link{parseFelixGxOutput}}) and \code{Felix32 v1.20} (see 
#' \code{\link{parseFelix32Output}}). The version is devined from the data 
#' structure and appropriate internal parsing functions called.
#' @param specFile Path to a \file{*.txt} file as a \code{\link{character}} 
#' object.
#' @return Returns a \code{\link{list}} with the follwoing keys:
#' \describe{
#'  \item{\code{Data}}{A \code{\link{data.frame}} representing the actual 
#'    spectrum with the columns \code{Time.in.sec} and 
#'    \code{Fluorescense.Intensity} (all \code{\link{numeric}}).}
#'  \item{\code{Data.Points}}{Number of data points in the spectrum as a 
#'    \code{\link{numeric}}. Ecquivalent to \code{\link{nrow}} of the 
#'    \code{link{data.frame}} in \code{Data}.}
#'  \item{\code{Max.Fluorescence.Intensity}}{\code{\link{numeric}} 
#'    representation of the maximal fluorescence intensity from \code{Data}.}
#'  \item{\code{Min.Fluorescence.Intensity}}{\code{\link{numeric}} 
#'    representation of the minimal fluorescence intensity from \code{Data}.}
#'  \item{\code{Max.Acquisition.Time.in.sec}}{\code{\link{numeric}} 
#'    representation of the maximal \code{Time.in.sec} from \code{Data}.}
#'  \item{\code{Min.Acquisition.Time.in.sec}}{\code{\link{numeric}} 
#'    representation of the minimal \code{Time.in.sec} from \code{Data}.}}
#' If provided by the data file parsed, an additional field is present:
#' \describe{
#'  \item{\code{File.Name}}{\code{\link{character}} representation of the 
#'    file name (as saved by the instrument).}
#' }
#' @seealso \code{flippant:::dithioniteFlippaseAssayInputValidation}
#' \code{\link{parseFelixGxOutput}} \code{\link{parseFelix32Output}}
#' \code{\link{parseManualOutput}}
#' @author Johannes Graumann
#' @keywords manip IO file
#' @examples
#' stop("Function is missing examples!")
#' @export
parseFluorometerOutput <- function(specFile=NA){
  #######################
  # Check Prerequisites #
  #######################
  if(is.na(specFile[1])){
    stop("'specFile' must be defined.")
  }
  if(!is.character(specFile)){
    stop("'specFile' must be of class 'character'.")
  }
  if(length(specFile) != 1){
    stop("'specFile' must be of length 1.")
  }
  if(!file.exists(specFile)){
    stop("'specFile' must represent an existing path.")
  }
  if(file.access(names=specFile,mode=4) != 0){
    stop("'specFile' must be readable.")
  }
  ##############
  # Processing #
  ##############
  # Aspirate the file
  ###################
  linesInSpecFile <- readLines(specFile)
  # Divine the output-producing fluorometer
  #########################################
  if(
    grepl(pattern="^<Trace>$",x=linesInSpecFile[1],ignore.case=TRUE) &
      grepl(pattern="^X\tY\t$",x=linesInSpecFile[4],ignore.case=TRUE) &
      grepl(pattern="^</Trace>$",x=tail(linesInSpecFile,n=1),ignore.case=TRUE)){
    formatOfSpecFile <- "FelixGXv4.1.0.3096"
  } else if(
    grepl(pattern="^1$",x=linesInSpecFile[1],ignore.case=TRUE) &
      grepl(pattern="^X\tY$",x=linesInSpecFile[4],ignore.case=TRUE) &
      grepl(pattern="^\\d+\\.{1,1}\\d+\t\\d+$",x=tail(linesInSpecFile,n=1),ignore.case=TRUE)
  ){
    formatOfSpecFile <- "Felix32v1.20"
  } else if(
    grepl(pattern="^Time \\(sec\\)\tBlank$",x=linesInSpecFile[1],ignore.case=TRUE)
  ){
    formatOfSpecFile <- "Manual"
  } else {
    stop("Unsupported data format in ",specFile)
  }
  # Extract data
  ##############
  if(formatOfSpecFile == "FelixGXv4.1.0.3096"){
    output <- parseFelixGxOutput(linesInSpecFile)
  } else if(formatOfSpecFile == "Felix32v1.20"){
    output <- parseFelix32Output(linesInSpecFile)
  } else if(formatOfSpecFile == "Manual"){
    output <- parseManualOutput(linesInSpecFile)
  }
  #################
  # Return result #
  #################
  invisible(output)
}