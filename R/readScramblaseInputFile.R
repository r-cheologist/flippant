#' @importFrom assertive assert_is_a_string
readScramblaseInputFile <- function(x){
  assert_is_a_string(x)
  if(!file.exists(x)){
    stop("'x' must represent an existing path.")
  }
  if(file.access(names=x,mode=4) != 0){
    stop("'x' must be readable.")
  }
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