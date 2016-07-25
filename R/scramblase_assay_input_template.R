#' @rdname scramblase_assay_plot
#' @importFrom magrittr %>%
#' @export
scramblase_assay_input_template <- function(path="scramblase_assay_input_template.txt"){
# Check input -------------------------------------------------------------
  path %>%
    assertive.strings::assert_is_a_non_empty_string() %>%
    assertive.files::assert_any_are_existing_files()

# Generate template data.frame --------------------------------------------
  # Create data structure in list form
  dataStructure <- list(
    ## Available fileds/columns
    columnName = c(
      "Path",
      "Protein Reconstituted (mg)",
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
            magrittr::extract2(x) %>%
            return()
        } else {
          dataStructure %>%
            magrittr::extract2(x) %>%
            magrittr::inset2(1,paste("#", .[1])) %>%
            return()
        }
      })
  # Assemble the data.frame
  commentedDataStructure %>%
    magrittr::extract(2:length(.)) %>%
    lapply(
      function(x){
        x %>%
          rbind() %>%
          as.data.frame(stringsAsFactors = FALSE) %>%
          return()
      }
    ) %>%
    plyr::rbind.fill() %>%
    magrittr::set_names(dataStructure %>% magrittr::extract2(1)) %>%
  # Write the table out
    utils::write.table(
      file=path,
      sep="\t",
      row.names=FALSE)
}

# Add global variables to satisfy 'R CMD check' ---------------------------
utils::globalVariables(".")
