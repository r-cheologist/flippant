#' @rdname scramblase_assay_plot
#' @importFrom assertive.files assert_any_are_existing_files
#' @importFrom assertive.strings assert_is_a_non_empty_string
#' @importFrom magrittr %>%
#' @importFrom magrittr extract
#' @importFrom magrittr extract2
#' @importFrom magrittr inset2
#' @importFrom plyr rbind.fill
#' @importFrom utils write.table
#' @export
scramblase_assay_input_template <- function(path="scramblase_assay_input_template.txt"){
# Check input -------------------------------------------------------------
  path %>%
    assert_is_a_non_empty_string() %>%
    assert_any_are_existing_files()

# Generate template data.frame --------------------------------------------
  # Create data structure in list form
  dataStructure <- list(
    ## Available fileds/columns
    columnName = c(
      "Path",
      "Protein in Reconstitution (mg)",
      "Fluorescence Assay Vol. w/o DT (ul)",
      "Fluorescence Assay Vol. with DT (ul)",
      "Lipid in Reconstitution (mmol)",
      "Timepoint of Measurement (s)",
      "Experiment",
      "Experimental Series"),
    ## Is the column required or facultative?
    columnRequired = c(
      rep(x="Required",times=2),
      rep(x="Facultative",times=6)),
    ## What's the expected data type?
    columnType = c(
      "character",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "character",
      "character"),
    ## Default values used if the column isn't defined
    columnDefault = c(
      rep(x="-",times=2),
      "2000",
      "2040",
      "0.0045",
      "400",
      "NA",
      "NA"))
  # Outcomment the appropriate lines
  commentedDataStructure <- dataStructure %>%
    names() %>%
    lapply(
      function(x){
        if(x == "columnName"){
          dataStructure %>%
            extract2(x) %>%
            return()
        } else {
          dataStructure %>%
            extract2(x) %>%
            inset2(1,paste("#", .[1])) %>%
            return()
        }
      })
  # Assemble the data.frame
  commentedDataStructure %>%
    extract(2:length(.)) %>%
    lapply(
      function(x){
        x %>%
          rbind() %>%
          as.data.frame(stringsAsFactors = FALSE) %>%
          return()
      }
    ) %>%
    rbind.fill() %>%
    set_names(dataStructure %>% extract2(1)) %>%
  # Write the table out
    write.table(
      file=path,
      sep="\t",
      row.names=FALSE)
}