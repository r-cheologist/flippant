#' @title parse_fluorimeter_output
#' @description Parse fluorimeter spectra
#' @details A function to read fluorimeter output directly. Intended as a helper
#' function to scramblase activity determinations from dithionite assays.
#' 
#' The function is currently capable to deal with input derived from 
#' QuantaMaster instruments (Photon Technology International, Inc., Edison, 
#' New Jersey)running software versions \code{FelixGX v4.1} 
#' (see \code{\link{parseFelixGxOutput}}) and \code{Felix32 v1.20} (see 
#' \code{\link{parseFelix32Output}}). The format used in a given file is devined
#' from the data structure and appropriate internal parsing functions are 
#' called.
#' 
#' If requested the time point of dithionite addition to a sample is determined 
#' using \pkg{wmtsa}-supplied methodology and the acquisition time reset
#' accordingly (\code{0} henceforth corresponds to the time of addition).
#' @param spec_file Path to a \file{*.txt} file as a \code{\link{character}} 
#' object.
#' @param adjust A \code{\link{logical}} indicating of whether (default) or not
#' acquisition time should be reset to have \code{0} (zero) coincide with the 
#' addition of dithionite (see 'Details' section).
#' @return Returns a \code{\link{list}} with the follwoing keys:
#' \describe{
#'  \item{\code{Data}}{A \code{\link{data.frame}} representing the actual 
#'    spectrum with the columns \code{Time.in.sec} and 
#'    \code{Fluorescence.Intensity} (all \code{\link{numeric}}).}
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
#' @seealso \code{scramblase_assay_input_validation},
#' \code{\link[flippant]{parseFelixGxOutput}}, 
#' \code{\link[flippant]{parseFelix32Output}},
#' \code{\link[flippant]{parse_manual_output}}
#' @author Johannes Graumann
#' @keywords manip IO file
#' @examples
#' # stop("Function is missing examples!")
#' @importFrom assertive assert_all_are_readable_files
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_a_string
#' @importFrom RcppRoll roll_mean
#' @importFrom wmtsa wavCWT
#' @importFrom wmtsa wavCWTPeaks
#' @importFrom wmtsa wavCWTTree
#' @export
parse_fluorimeter_output <- function(
  spec_file = NULL,
  adjust = TRUE){
  #######################
  # Check Prerequisites #
  #######################
  assert_is_a_string(spec_file)
  assert_all_are_readable_files(spec_file)
  assert_is_a_bool(adjust)
  ##############
  # Processing #
  ##############
  # Aspirate the file
  ###################
  linesInSpecFile <- readLines(spec_file)
  # Divine the output-producing fluorimeter
  #########################################
  if(
    grepl(pattern="^<Trace>\\s*$",x=linesInSpecFile[1],ignore.case=TRUE) &
      grepl(pattern="^X\\tY\\s*$",x=linesInSpecFile[4],ignore.case=TRUE) &
      grepl(pattern="^</Trace>\\s*$",x=tail(linesInSpecFile,n=1),ignore.case=TRUE)){
    formatOfSpecFile <- "FelixGXv4.1.0.3096"
  } else if(
    grepl(pattern="^1\\s*$",x=linesInSpecFile[1],ignore.case=TRUE) &
      grepl(pattern="^X\\tY\\s*$",x=linesInSpecFile[4],ignore.case=TRUE) &
      grepl(pattern="^\\d+\\.{1,1}\\d+\t\\d+\\s*$",x=tail(linesInSpecFile,n=1),ignore.case=TRUE)
  ){
    formatOfSpecFile <- "Felix32v1.20"
  } else if(
    grepl(pattern="^Time \\(sec\\)\tFluorescense Intensity\\s*$",x=linesInSpecFile[1],ignore.case=TRUE)
  ){
    formatOfSpecFile <- "Manual"
  } else {
    stop("Unsupported data format in ",spec_file)
  }
  # Extract data
  ##############
  if(formatOfSpecFile == "FelixGXv4.1.0.3096"){
    output <- parseFelixGxOutput(linesInSpecFile)
  } else if(formatOfSpecFile == "Felix32v1.20"){
    output <- parseFelix32Output(linesInSpecFile)
  } else if(formatOfSpecFile == "Manual"){
    output <- parse_manual_output(linesInSpecFile)
  }
  # Determine timepoint of dithionite addition and adjust time axis accordingly
  #############################################################################
  # Following http://stackoverflow.com/a/24560355/2103880
  if(adjust){
    # Smooth the intensity using a simple rolling mean
    n <- 10
    smoothedTmpData <- RcppRoll::roll_mean(output$Data$Fluorescence.Intensity, n)
    # Wavelet-based peak detection
    gradient <- -diff(smoothedTmpData)
    cwt <- wavCWT(gradient)
    tree <- wavCWTTree(cwt)
    peaks <- wavCWTPeaks(tree)
    # Functionality appears to have memory mapping issues ...
    gc()
    # Safeguard against artifact peaks (based on incomplete understanding of wmtsa)
    peakMax <- peaks$x[which.max(peaks$y)]
    # Reset acquisition time
    if(length(peakMax) < 1){
      stop("No ditionite addition detected in '",spec_file,"'.")
    }
    zeroTimePoint <- output$Data$Time.in.sec[peakMax[1]]
    output$Data$Time.in.sec <- output$Data$Time.in.sec - zeroTimePoint
  }
  output$Max.Acquisition.Time.in.sec <- max(output$Data$Time.in.sec,na.rm=TRUE)
  output$Min.Acquisition.Time.in.sec <- min(output$Data$Time.in.sec,na.rm=TRUE)
  #################
  # Return result #
  #################
  invisible(output)
}