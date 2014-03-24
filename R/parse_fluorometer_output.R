#' @title parse_fluorometer_output
#' @description Parse fluorometer spectra
#' @details A function to read fluorometer output directly. Intended as a helper
#' function to flippase activity determinations from dithionite assays.
#' 
#' The function is currently capable to deal with input derived from QuantMaster
#' instruments running software versions \code{FelixGX v4.1} and 
#' \code{Felix32 v1.20}. The version is devined from the data structure and
#' appropriate internal parsing functions called.
#' @param spec_file Path to a \file{*.txt} file as a \code{\link{character}} 
#' object.
#' @return Returns a \code{\link{list}} with the follwoing keys:
#' \describe{
#'  \item{\code{Data}}{A \code{\link{data.frame}} representing the actual 
#'    spectrumwith the columns \code{Time.in.sec} and \code{Fluorescense.Intensity}
#'    (all \code{\link{numeric}}).}
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
#' @seealso \code{\link{dithionite_flippase_assay_input_validation}} 
#' \code{\link{parse_felixgx_output}} \code{\link{parse_felix32_output}}
#' \code{\link{parse_manual_output}}
#' @author Johannes Graumann
#' @export
#' @keywords manip IO file
#' @examples
#' stop("Function is missing examples!")
parse_fluorometer_output <- function(spec_file=NA){
  #######################
  # Check Prerequisites #
  #######################
  if(is.na(spec_file[1])){
    stop("'spec_file' must be defined.")
  }
  if(!is.character(spec_file)){
    stop("'spec_file' must be of class 'character'.")
  }
  if(length(spec_file) != 1){
    stop("'spec_file' must be of length 1.")
  }
  if(!file.exists(spec_file)){
    stop("'spec_file' must represent an existing path.")
  }
  if(file.access(names=spec_file,mode=4) != 0){
    stop("'spec_file' must be readable.")
  }
  ##############
  # Processing #
  ##############
  # Aspirate the file
  ###################
  lines_in_spec_file <- readLines(spec_file)
  # Divine the output-producing fluorometer
  #########################################
  if(
    grepl(pattern="^<Trace>$",x=lines_in_spec_file[1],ignore.case=TRUE) &
      grepl(pattern="^X\tY\t$",x=lines_in_spec_file[4],ignore.case=TRUE) &
      grepl(pattern="^</Trace>$",x=tail(lines_in_spec_file,n=1),ignore.case=TRUE)){
    format_of_spec_file <- "FelixGXv4.1.0.3096"
  } else if(
    grepl(pattern="^1$",x=lines_in_spec_file[1],ignore.case=TRUE) &
      grepl(pattern="^X\tY$",x=lines_in_spec_file[4],ignore.case=TRUE) &
      grepl(pattern="^\\d+\\.{1,1}\\d+\t\\d+$",x=tail(lines_in_spec_file,n=1),ignore.case=TRUE)
  ){
    format_of_spec_file <- "Felix32v1.20"
  } else if(
    grepl(pattern="^Time \\(sec\\)\tBlank$",x=lines_in_spec_file[1],ignore.case=TRUE)
  ){
    format_of_spec_file <- "Manual"
  } else {
    stop("Unsupported data format in ",spec_file)
  }
  # Extract data
  ##############
  if(format_of_spec_file == "FelixGXv4.1.0.3096"){
    output <- parse_felixgx_output(lines_in_spec_file)
  } else if(format_of_spec_file == "Felix32v1.20"){
    output <- parse_felix32_output(lines_in_spec_file)
  } else if(format_of_spec_file == "Manual"){
    output <- parse_manual_output(lines_in_spec_file)
  }
  #################
  # Return result #
  #################
  invisible(output)
}