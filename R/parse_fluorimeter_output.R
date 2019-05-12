#' @title parse_fluorimeter_output
#' @description Parse fluorimeter spectra
#' @details A function to read fluorimeter output directly. Intended as a helper
#' function to scramblase activity determinations from dithionite assays.
#' 
#' The function is currently capable to deal with input derived from 
#' QuantaMaster instruments (Photon Technology International, Inc., Edison, 
#' New Jersey)running software versions \code{FelixGX v4.1} 
#' (see \code{\link{parse_felix_gx_output}}), \code{Felix32 v1.20} (see 
#' \code{\link{parse_felix_32_output}}) as well as Horiba fluorimeters 
#' (HORIBA Europe GmbH, Oberursel, Germany) using \code{FluorS Essence v3.8}.
#' The format used in a given file is divined from the data structure and
#' appropriate internal parsing functions are called.
#' 
#' If requested the time point of dithionite addition to a sample is determined 
#' using \pkg{wmtsa}-supplied methodology and the acquisition time reset
#' accordingly (\code{0} henceforth corresponds to the time of addition).
#' @param spec_file Path to a \file{*.txt} file as a \code{\link{character}} 
#' object.
#' @param determine_zero_time A \code{\link{logical}} indicating whether
#' (default) or not the timepoint of dithionite addition should be determined
#' using \pkg{wmtsa}-derived functionality.
#' @param adjust A \code{\link{logical}} indicating of whether (default) or not
#' acquisition time should be reset to have \code{0} (zero) coincide with the 
#' addition of dithionite (see 'Details' section).
#' @param file_type A string specifying whether or the file was created using
#' Felix GX or Felix 32 or FluorS Essence v3.8 or is a "manual" tab delimited
#' file.
#' @return A data frame with two columns:
#' \describe{
#' \item{Time.in.sec}{Numeric. Number of seconds since the start of experiment.}
#' \item{Fluorescence.Intensity}{Numeric. Intensity of fluorescence (relative 
#' scale, no official unit).}
#' }
#' If \code{determine_zero_time} and/or \code{adjust} are set to \code{TRUE},
#' the return value will have an attribute \code{ZeroTimePoint} corresponding to
#' the determined time point of dithionite addition (always \code{0} (zero)
#' where \code{adjust == TRUE}).
#' 
#' For Felix GX, if the file contains the information, the return value will  
#' also have an attribute \code{WavelengthsInNanometres}, which contains the 
#' excitation and emission wavelengths.
#' @seealso \code{scramblase_assay_input_validation},
#' \code{\link[flippant]{parse_felix_gx_output}}, 
#' \code{\link[flippant]{parse_felix_32_output}},
#' \code{\link[flippant]{parse_FluorS_Essence_3.8_output}},
#' \code{\link[flippant]{parse_manual_output}}
#' @author Johannes Graumann
#' @keywords manip IO file
#' @examples
#' library(magrittr)
#' # Extract example data
#' analysis_dir <- file.path(tempdir(), "flippant-case-study")
#' fluor_file <- extract_case_study_data(analysis_dir, "wildtypeEx1_0.txt")
#' # Parse an exemplary file
#' parse_fluorimeter_output(fluor_file) %>%
#'   str()
#' @export
parse_fluorimeter_output <- function(
  spec_file = NULL,
  determine_zero_time = TRUE,
  adjust = TRUE,
  file_type = c("auto", "FelixGXv4.1.0.3096", "Felix32v1.20", "FluorSEssencev3.8", "manual")) {
  # Check prerequisites -----------------------------------------------------
  assertive.types::assert_is_a_string(spec_file)
  assertive.files::assert_all_are_readable_files(spec_file, warn_about_windows = FALSE)
  assertive.types::assert_is_a_bool(determine_zero_time)
  assertive.types::assert_is_a_bool(adjust)

  # Processing --------------------------------------------------------------

  # Divine the output-producing fluorimeter
  #########################################
  file_type <- match.arg(file_type)
  if (file_type == "auto") {
    first_lines <- readLines(spec_file, 4)
    formatOfSpecFile <- if (
      grepl(pattern = "^<Trace>\\s*$",x = first_lines[1],ignore.case = TRUE) &
        grepl(pattern = "^X\\tY\\s*$",x = first_lines[4],ignore.case = TRUE)) {
      "FelixGXv4.1.0.3096"
    } else if (
      grepl(pattern = "^1\\s*$",x = first_lines[1],ignore.case = TRUE) &
        grepl(pattern = "^X\\tY\\s*$",x = first_lines[4],ignore.case = TRUE)) {
      "Felix32v1.20"
    } else if (
      grepl(pattern = "^Time\\tS1\\s*$",x = first_lines[1],ignore.case = TRUE) &
        grepl(pattern = "^s\\tCPS\\s*$",x = first_lines[2],ignore.case = TRUE)) {
      "FluorSEssencev3.8"
    } else if (
      grepl(pattern = "^Time.\\(sec\\)\tFluorescense.Intensity\\s*$",
            x = first_lines[1],ignore.case = TRUE)) {
      "manual"
    } else {
      stop("Unsupported data format in ",spec_file)
    }
  }
  # Extract data
  ##############
  output <- switch(
    formatOfSpecFile,
    FelixGXv4.1.0.3096 = read_felix_gx(spec_file),
    Felix32v1.20 = {
      lines <- readLines(spec_file)
      y <- parse_felix_32_output(lines)
      structure(y$Data, File.Name = y$File.Name) 
    },
    FluorSEssencev3.8 = {
      lines <- readLines(spec_file)
      y <- parse_FluorS_Essence_3.8_output(lines)
      structure(y$Data, File.Name = basename(spec_file))
    },
    manual = {
      lines <- readLines(spec_file)
      y <- parse_manual_output(lines)
      structure(y$Data, File.Name = y$File.Name) 
    }
  )
  # Determine timepoint of dithionite addition and adjust time axis accordingly
  #############################################################################
  if (adjust || determine_zero_time) {
    # Following http://stackoverflow.com/a/24560355/2103880
    # Smooth the intensity using a simple rolling mean
    n <- 10
    smoothedTmpData <- RcppRoll::roll_mean(output$Fluorescence.Intensity, n)
    # Wavelet-based peak detection
    gradient <- -diff(smoothedTmpData)
    cwt <- wmtsa::wavCWT(gradient)
    tree <- wmtsa::wavCWTTree(cwt)
    peaks <- wmtsa::wavCWTPeaks(tree)
    # Functionality appears to have memory mapping issues ...
    gc()
    # Safeguard against artifact peaks (based on incomplete understanding of wmtsa)
    peakMax <- peaks$x[which.max(peaks$y)]
    # Reset acquisition time
    if (length(peakMax) < 1) {
      stop("No dithionite addition detected in '",spec_file,"'.")
    }
    zeroTimePoint <- output$Time.in.sec[peakMax[1]]
    if (adjust) {
      output$Time.in.sec <- output$Time.in.sec - zeroTimePoint
      zeroTimePoint <- 0
    }
    attr(output, "ZeroTimePoint") <- zeroTimePoint
  }

  # Return results ----------------------------------------------------------
  invisible(output)
}