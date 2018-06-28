#' @rdname scramblase_assay_plot
#' @importFrom magrittr %>%
#' @export
scramblase_assay_input_template <- function(path="scramblase_assay_input_template.txt", input_directory = NULL){
# Check input -------------------------------------------------------------
  path %>%
    assertive.strings::assert_is_a_non_empty_string() %>%
    assertive.files::is_existing_file() %>%
    assertive.base::assert_all_are_false()
  if(!is.null(input_directory))
  {
    input_directory %>%
      assertive.strings::assert_is_a_non_empty_string() %>%
      assertive.files::assert_all_are_dirs()
  }

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
  commentedDataStructure %<>%
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
    magrittr::set_names(dataStructure %>% magrittr::extract2(1))
  # Get files in input_directory for prepopulation
  if(!is.null(input_directory))
  {
    paths <- input_directory %>%
        list.files()
    tmp_content <- matrix(
      nrow = length(paths),
      ncol = ncol(commentedDataStructure)) %>%
      as.data.frame() %>%
      magrittr::set_colnames(names(commentedDataStructure))
    tmp_content$Path <- paths
    tmp_content[is.na(tmp_content)] <- ''
    commentedDataStructure %<>%
      plyr::rbind.fill(
        tmp_content)
  }
  
  # Write the table out
  commentedDataStructure %>%
    utils::write.table(
      file=path,
      sep="\t",
      row.names=FALSE)
}

# Add global variables to satisfy 'R CMD check' ---------------------------
utils::globalVariables(".")
