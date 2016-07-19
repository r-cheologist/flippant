#' Extract the case study dataset
#' 
#' Extracts the data files used by the case study from the zip archive.
#' @param files A character vector of files to extract, or \code{NULL} to 
#' extract all files.  Passed to \code{\link[utils]{unzip}}.
#' @param exdir A string giving the path to the extraction directory. Passed to
#' \code{\link[utils]{unzip}}.
#' @return A character vector of the extracted file paths is invisibly returned.
#' @seealso \code{\link[utils]{unzip}}
#' @author Richard Cotton
#' @examples 
#' extract_case_study_data(tempfile("flippant"))
#' @importFrom utils unzip
#' @export
extract_case_study_data <- function(exdir = ".", files = NULL)
{
  "extdata/PloierEtAl_Data.zip" %>% 
    system.file(package = "flippant") %>% 
    unzip(files = files, exdir = exdir)
}
