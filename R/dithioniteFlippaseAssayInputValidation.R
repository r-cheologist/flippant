#' @importFrom assertive assert_is_data.frame
#' @importFrom assertive assert_is_a_bool
dithioniteFlippaseAssayInputValidation <- function(x,scaleTo,forceThroughOrigin){
  # Check x
  ## General DF characteristics
  assert_is_data.frame(x)
  if(nrow(x)==0){
    stop("'x' must have rows.")
  }
  if(any(is.na(x))){
    stop("'x' cannot contain 'NA'.")
  }
  ## Enforce use of 'numeric' over 'integer'
  integerColumns <- names(x)[sapply(names(x),function(y){is.integer(x[[y]])})]
  if(length(integerColumns) > 0){
    for(column in integerColumns){
      x[[column]] <- as.numeric(x[[column]])
    }
  }
  ## Required parameters
  requiredColumnsInX <- list(
    Name = c(
      "Path",
      "Protein in Reconstitution (mg)"),
    Class = c(
      "character",
      "numeric"))
  if(!all( requiredColumnsInX$Name %in% names(x))){
    stop(
      "'x' must hold at least the following columns: '",
      paste(requiredColumnsInX,collapse="', '"),
      "'.")
  }
  if(!identical(
    unname(vapply(x[requiredColumnsInX$Name],class,c(A="A"))),
    requiredColumnsInX$Class)){
    stop(
      "Required columns '",
      paste(requiredColumnsInX$Name,collapse="', '"),
      "' must be of classes '",
      paste(requiredColumnsInX$Class,collapse="', '"),
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
  facultativeColumnsInX <- list(
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
      2000.0,
      2040.0,
      0.0045,
      400.0,
      NA_character_,
      NA_character_))
  missingFacultativeColumnsInX <- which(!(facultativeColumnsInX$Name %in% names(x)))
  if(length(missingFacultativeColumnsInX) != 0){
    for (y in missingFacultativeColumnsInX){
      if(facultativeColumnsInX$Name[y] %in% c("Experiment","Experimental Series")) {
        toBeAddedOn <- facultativeColumnsInX$Default[[y]]
      } else {
        message(
          "Providing missing column '",
          facultativeColumnsInX$Name[y],
          "' from defaults (",
          facultativeColumnsInX$Default[[y]],
          "). Make sure this is correct.")
        toBeAddedOn <- facultativeColumnsInX$Default[[y]]
      }
      output <- cbind(
        x,
        rep(x=toBeAddedOn,times=nrow(x)),
        stringsAsFactors=FALSE)
      names(output)[ncol(output)] <- facultativeColumnsInX$Name[y]
      x <- output
    }
  }
  
  if(!identical(
    unname(vapply(x[facultativeColumnsInX$Name],class,c(A="A"))),
    facultativeColumnsInX$Class)){
    stop(
      "Facultative columns '",
      paste(facultativeColumnsInX$Name,collapse="', '"),
      "' must be of classes '",
      paste(facultativeColumnsInX$Class,collapse="', '"),
      "'.")
  }
  ## Check "Timepoint of Measurement (s)" consistency
  if(length(unique(x$"Timepoint of Measurement (s)")) != 1){
    stop("Column 'Timepoint of Measurement (s)' contains multiple values. Exiting.")
  }
  
  # Check scaleTo
  scaleTo <- match.arg(
    arg=scaleTo,
    choices=c("model","data"),
    several.ok=FALSE)
  if(scaleTo == "model"){
    message("Data will be scaled to the plateau of its monoexponential fit.")
  }
  # Check forceThroughOrigin
  assert_is_a_bool(forceThroughOrigin)
  # Return
  return(list(x=x,scaleTo=scaleTo,forceThroughOrigin=forceThroughOrigin))
}