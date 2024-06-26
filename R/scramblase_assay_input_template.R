#' @rdname scramblase_assay_plot
#' @importFrom magrittr %>%
#' @export
scramblase_assay_input_template <- function(
  path = "scramblase_assay_input_template.txt", input_directory = NULL,
  overwrite = FALSE){
# Check input -------------------------------------------------------------
  assert_is_a_bool(overwrite)
  assert_is_a_non_empty_string(path)
  if (!overwrite) {
   assert_all_are_false(
      assert_any_are_existing_files(path))
  }
  if (!is.null(input_directory)) {
    assert_is_a_non_empty_string(input_directory)
    assert_all_are_dirs(input_directory)
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
  if (!is.null(input_directory))
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
    utils::write.table(file = path, quote = FALSE, sep = "\t",
                       row.names = FALSE)
  # Silently return path
  invisible(path)
}

# Add global variables to satisfy 'R CMD check' ---------------------------
utils::globalVariables(".")
