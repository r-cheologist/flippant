dithionite_flippase_assay_input_validation <- function(x,scale_to){
  # Check x
  ## General DF characteristics
  if(!is.data.frame(x)){
    stop("'x' must be of class 'data.frame'.")
  }
  if(nrow(x)==0){
    stop("'x' must have rows.")
  }
  if(any(is.na(x))){
    stop("'x' cannot contain 'NA'.")
  }
  ## Required parameters
  required_columns_in_x <- list(
    Name = c(
      "Path",
      "Protein in Reconstitution (mg)"),
    Class = c(
      "character",
      "numeric"))
  if(!all( required_columns_in_x$Name %in% names(x))){
    stop(
      "'x' must hold at least the following columns: '",
      paste(required_columns_in_x,collapse="', '"),
      "'.")
  }
  if(!identical(
    unname(vapply(x[required_columns_in_x$Name],class,c(A="A"))),
    required_columns_in_x$Class)){
    stop(
      "Required columns '",
      paste(required_columns_in_x$Name,collapse="', '"),
      "' must be of classes '",
      paste(required_columns_in_x$Class,collapse="', '"),
      "'.")
  }
  ## Check paths
  if(!all(file.exists(x$Path))){
    stop("All entries in column 'Path' must refer to existing files.")
  }
  if(any(file.access(x$Path,mode=4) == -1)){
    stop("All entries in column 'Path' must refer to existing files.")
  }
  ## Facultative parameters
  facultative_columns_in_x <- list(
    Name = c(
      "Fluorescence Assay Vol. w/o DT (ul)",
      "Fluorescence Assay Vol. with DT (ul)",
      "Egg PC in Reconstitution (mmol)",
      "Timepoint of Measurement (s)",
      "Experiment",
      "Experimental Series"),
    Class = c(
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "character",
      "character"),
    Default = list(
      2000,
      2040,
      0.0045,
      NA,
      NA_character_,
      NA_character_))
  missing_facultative_columns_in_x <- which(!(facultative_columns_in_x$Name %in% names(x)))
  if(length(missing_facultative_columns_in_x) != 0){
    for (y in missing_facultative_columns_in_x){
      if(facultative_columns_in_x$Name[y] == "Timepoint of Measurement (s)"){
        warning(
          "Providing missing column '",
          facultative_columns_in_x$Name[y],
          "' from spectra ('Path').")
        to_be_added_on <- timepoint_of_measurement(x$Path)
      } else if(facultative_columns_in_x$Name[y] %in% c("Experiment","Experimental Series")) {
        to_be_added_on <- facultative_columns_in_x$Default[[y]]
      } else {
        warning(
          "Providing missing column '",
          facultative_columns_in_x$Name[y],
          "' from defaults (",
          facultative_columns_in_x$Default[[y]],
          "). Make sure this is correct.")
        to_be_added_on <- facultative_columns_in_x$Default[[y]]
      }
      output <- cbind(
        x,
        rep(x=to_be_added_on,times=nrow(x)),
        stringsAsFactors=FALSE)
      names(output)[ncol(output)] <- facultative_columns_in_x$Name[y]
      x <- output
    }
  }
  if(!identical(
    unname(vapply(x[facultative_columns_in_x$Name],class,c(A="A"))),
    facultative_columns_in_x$Class)){
    stop(
      "Facultative columns '",
      paste(facultative_columns_in_x$Name,collapse="', '"),
      "' must be of classes '",
      paste(facultative_columns_in_x$Class,collapse="', '"),
      "'.")
  }
  ## Check "Timepoint of Measurement (s)" consistency
  if(length(unique(x$"Timepoint of Measurement (s)")) != 1){
    stop("Column 'Timepoint of Measurement (s)' contains multiple values. Exiting.")
  }
  
  # Check scale_to
  scale_to <- match.arg(
    arg=scale_to,
    choices=c("model","data"),
    several.ok=FALSE)
  if(scale_to == "model"){
    warning("Data will be scaled to the plateau of its monoexponential fit.")
  }
  # Return
  return(list(x=x,scale_to=scale_to))
}