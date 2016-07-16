#' @importFrom magrittr %>%
read_scramblase_input_file <- function(x){
  x %>%
    assertive.strings::assert_is_a_non_empty_string() %>%
    assertive.files::assert_all_are_readable_files() %>%
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