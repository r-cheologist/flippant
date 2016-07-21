#' Read Felix GX data files
#' 
#' Reads spectra from files provided by a QuantaMaster fluorimeter 
#' (Photon Technology International, Inc., Edison, New Jersey) using 
#' \code{FelixGX v4.1}
#' @param data_file A string giving a path to the data file.
#' @return A data frame with two columns:
#' \describe{
#' \item{Time.in.sec}{Numeric. Number of seconds since the start of experiment.}
#' \item{Fluorescence.Intensity}{Numeric. Intensity of fluorescence (relative 
#' scale, no official unit).}
#' }
#' If the file contains the information, the return value will also have an
#' attribute \code{WavelengthsInNanometres}, which contains the excitation and 
#' emission wavelengths.
#' @noRd
read_felix_gx <- function(data_file)
{
  metadata <- readLines(data_file, 3)
  expected_nrows <- as.integer(metadata[2])
  data <- data.table::fread(
    data_file, 
    nrows = expected_nrows, 
    skip = 4, 
    col.names = c("Time.in.sec", "Fluorescence.Intensity"), 
    data.table = FALSE)
  matches <- stringi::stri_match_first_regex(metadata[3], "([0-9]+):([0-9]+)")
  if(!is.na(matches[1, 1]))
  {
    attr(data, "WavelengthsInNanometres") <- c(
      excitation = as.numeric(matches[, 2]), 
      emission   = as.numeric(matches[, 3])
    )
  }
  data
}
