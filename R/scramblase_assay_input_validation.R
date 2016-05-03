#' @importFrom assertive assert_has_rows
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_data.frame
scramblase_assay_input_validation <- function(
  x,
  scale_to,
  generation_of_algorithm,
  force_through_origin,
  split_by_experiment,
  verbose=TRUE){
  # Check verbose
  assert_is_a_bool(verbose)
  # Check x
  ## General DF characteristics
  assert_is_data.frame(x)
  assert_has_rows(x)
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
      "Lipid in Reconstitution (mmol)",
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
        if(verbose){
          message(
            "Providing missing column '",
            facultativeColumnsInX$Name[y],
            "' from defaults (",
            facultativeColumnsInX$Default[[y]],
            "). Make sure this is correct.")
        }
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

  # Check scale_to
  scale_to <- match.arg(
    arg=scale_to,
    choices=c("model","data"),
    several.ok=FALSE)
  if(scale_to == "model" & verbose){
    message("Data will be scaled to the plateau of its monoexponential fit.")
  }
  # Check generation_of_algorithm
  generation_of_algorithm <- as.character(generation_of_algorithm)
  if(identical(generation_of_algorithm, c("2","1"))){
    generation_of_algorithm <- "2"
  }
  generation_of_algorithm <- match.arg(
    arg = generation_of_algorithm,
    choices = c("second","first", "2", "1"),
    several.ok = FALSE)
  if(generation_of_algorithm %in% c("first","1")){
    generation_of_algorithm <- 1
  } else if(generation_of_algorithm %in% c("second", "2")){
    generation_of_algorithm <- 2
  }
  if(verbose){
    message("Data will be fitted to algorithm generation ", generation_of_algorithm, ".")
  }
  # Check force_through_origin
  assert_is_a_bool(force_through_origin)
  # Check split_by_experiment
  assert_is_a_bool(split_by_experiment)
  # Return
  return(
    list(
      x = x,
      scale_to = scale_to,
      generation_of_algorithm = generation_of_algorithm,
      force_through_origin = force_through_origin,
      split_by_experiment = split_by_experiment))
}