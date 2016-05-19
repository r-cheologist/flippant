#' @importFrom assertive assert_all_are_readable_files
#' @importFrom assertive assert_is_a_non_empty_string
#' @importFrom utils read.delim
read_scramblase_input_file <- function(x){
  assert_is_a_non_empty_string(x)
  assert_all_are_readable_files(x)
  output <- read.delim(
    file=x,
    header=TRUE,
    sep="\t",
    quote="\"'",
    fill=TRUE,
    comment.char="#",
    stringsAsFactors=FALSE,
    check.names=FALSE)
  return(output)
}