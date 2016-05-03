#' @rdname scramblase_assay_plot
#' @importFrom assertive assert_any_are_existing_files
#' @importFrom assertive assert_is_a_non_empty_string
#' @importFrom plyr rbind.fill
#' @export
scramblaseAssayInputTemplate <- function(path="scramblaseAssayInputTemplate.txt"){
# Check input -------------------------------------------------------------
  assert_is_a_non_empty_string(path)
  assert_any_are_existing_files(path)

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
  commentedDataStructure <- lapply(
    names(dataStructure),
    function(x){
      if(x == "columnName"){
        return(dataStructure[[x]])
      } else {
        tmpVec <- dataStructure[[x]]
        tmpVec[1] <- paste("#",tmpVec[1])
        return(tmpVec)
      }
    })
  # Assemble the data.frame
  outputDataFrameAsList <- lapply(
    commentedDataStructure[2:length(commentedDataStructure)],
    function(x){as.data.frame(rbind(x),stringsAsFactors=FALSE)})
  outputDataFrame <- plyr::rbind.fill(outputDataFrameAsList)
  names(outputDataFrame) <- dataStructure[[1]]

# Write the file out ------------------------------------------------------
  write.table(x=outputDataFrame,file=path,sep="\t",row.names=FALSE)

}