#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
extract_fluorescence_extrema <- function(
  x,
  timepoint_of_measurement,
  n_averaging = 10
) {
  # Check prerequisites -----------------------------------------------------
  assert_has_all_attributes(x, "ZeroTimePoint")
  assert_is_a_number(attr(x, "ZeroTimePoint"))
  if (!is.null(timepoint_of_measurement))
      assert_is_a_number(timepoint_of_measurement)
  assert_all_are_in_closed_range(
    timepoint_of_measurement,
    lower = min(x$Time.in.sec, na.rm = TRUE),
    upper = max(x$Time.in.sec, na.rm = TRUE))
  is_a_number(n_averaging)
  assert_all_are_positive(n_averaging)
  assert_all_are_whole_numbers(n_averaging)

  # Extract baseline fluorescence -------------------------------------------
  min_fluorescence_to_extract <- x %>%
    magrittr::extract(.$Time.in.sec < attr(x, "ZeroTimePoint"), ) %>% 
    utils::tail(n_averaging)
  if (nrow(min_fluorescence_to_extract) < n_averaging) {
    stop("Insufficient baseline data points (n = ",
         nrow(min_fluorescence_to_extract),
         "to extract where 'n_averaging' equals '", n_averaging,
         "' and 'timepoint_of_measurment' is '", timepoint_of_measurement, "'.")
  }

  # Extract endpoint fluorescence -------------------------------------------
  max_fluorescence_to_extract <- x %>%
    magrittr::extract(.$Time.in.sec >= timepoint_of_measurement, ) %>%
    utils::head(n_averaging)
  if (nrow(max_fluorescence_to_extract) < n_averaging) {
    stop("Insufficient endpoint data points (n = ",
         nrow(max_fluorescence_to_extract),
         ") to extract where 'n_averaging' equals '", n_averaging,
         "' and 'timepoint_of_measurment' is '", timepoint_of_measurement, "'.")
  }

  # Averaging ---------------------------------------------------------------
  fluorescence_extrema <- list(
    min_fluorescence_to_extract, max_fluorescence_to_extract) %>%
    lapply(magrittr::extract2, "Fluorescence.Intensity") %>%
    sapply(stats::median, na.rm = TRUE) %>%
    unlist() %>%
    magrittr::set_names(c("Baseline.Fluorescence", "Endpoint.Fluorescence"))

  # Return ------------------------------------------------------------------
  return(fluorescence_extrema)
}
