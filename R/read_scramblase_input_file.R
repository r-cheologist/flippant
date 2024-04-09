#' @importFrom magrittr %>%
read_scramblase_input_file <- function(x){
  x %>%
    assert_is_a_non_empty_string() %>%
    assert_all_are_readable_files(warn_about_windows = FALSE) %>%
    utils::read.delim(
      header=TRUE,
      sep="\t",
      quote="\"'",
      fill=TRUE,
      comment.char="#",
      stringsAsFactors=FALSE,
      check.names=FALSE) %>%
    return()
}