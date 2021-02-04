#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
scramblase_assay_input_validation <- function(
  x,
  scale_to,
  ppr_scale_factor,
  generation_of_algorithm,
  force_through_origin,
  split_by_experiment,
  n_averaging,
  r_bar,
  sigma_r_bar,
  protein_is_factor = FALSE,
  verbose = TRUE){
  # Check switches
  assertive.types::assert_is_a_bool(protein_is_factor)
  assertive.types::assert_is_a_bool(verbose)
  # Check x
  ## General DF characteristics
  assertive.types::assert_is_data.frame(x)
  assertive.properties::assert_has_rows(x)
  if (any(is.na(x))) stop("'x' cannot contain 'NA'.")
  ## Enforce use of 'numeric' over 'integer'
  integerColumns <- names(x)[sapply(names(x),function(y){is.integer(x[[y]])})]
  if (length(integerColumns) > 0) {
    for (column in integerColumns)  x[[column]] <- as.numeric(x[[column]])
  }
  ## Required parameters
  requiredColumnsInX <- list(
    Name = c(
      "Path",
      "Protein Reconstituted (mg)"),
    Class = c(
      "character",
      ifelse(protein_is_factor, "character", "numeric")))
  if (!all( requiredColumnsInX$Name %in% names(x))) {
    stop(
      "'x' must hold at least the following columns: '",
      paste(requiredColumnsInX,collapse = "', '"),
      "'.")
  }
  if (!identical(
    unname(vapply(x[requiredColumnsInX$Name],class,c(A = "A"))),
    requiredColumnsInX$Class)) {
    stop(
      "Required columns '",
      paste(requiredColumnsInX$Name,collapse = "', '"),
      "' must be of classes '",
      paste(requiredColumnsInX$Class,collapse = "', '"),
      "'.")
  }
  ## Check paths
  if (!all(file.exists(x$Path))) {
    stop("All entries in column 'Path' must refer to existing files.")
  }
  if (any(file.access(x$Path,mode = 4) == -1)) {
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
      ifelse(
        "Experiment" %in% names(x) && is.factor(x$Experiment),
        "factor", "character"),
      ifelse(
        "Experimental Series" %in% names(x) &&
          is.factor(x[["Experimental Series"]]),
        "factor", "character")),
    Default = list(
      2000.0,
      2040.0,
      0.0045,
      400.0,
      NA_character_,
      NA_character_))
  missingFacultativeColumnsInX <- which(!(facultativeColumnsInX$Name %in% names(x)))
  if (length(missingFacultativeColumnsInX) != 0) {
    for (y in missingFacultativeColumnsInX) {
      if (facultativeColumnsInX$Name[y] %in% c("Experiment","Experimental Series")) {
        toBeAddedOn <- facultativeColumnsInX$Default[[y]]
      } else {
        if (verbose) {
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
        rep(x = toBeAddedOn, times = nrow(x)),
        stringsAsFactors = FALSE)
      names(output)[ncol(output)] <- facultativeColumnsInX$Name[y]
      x <- output
    }
  }
  
  if (!identical(
    unname(vapply(x[facultativeColumnsInX$Name],class,c(A = "A"))),
    facultativeColumnsInX$Class)) {
    stop(
      "Facultative columns '",
      paste(facultativeColumnsInX$Name,collapse = "', '"),
      "' must be of classes '",
      paste(facultativeColumnsInX$Class,collapse = "', '"),
      "'.")
  }

  # Check scale_to
  scale_to %<>%
    match.arg(
      choices = c("model","data"),
      several.ok = FALSE)
  if (scale_to == "model" & verbose) {
    message("Data will be scaled to the plateau of its monoexponential fit.")
  }
  # Check ppr_scale_factor
  if (!is.null(ppr_scale_factor)) {
    ppr_scale_factor %>%
      assertive.types::assert_is_a_number()
    if (verbose) {
      message("PPR will be scaled by a factor of ", ppr_scale_factor, ".")
    }
  }
  # Check generation_of_algorithm
  generation_of_algorithm %<>% 
    as.character()
  if (identical(generation_of_algorithm, c("2","1"))) {
    generation_of_algorithm <- "2"
  }
  generation_of_algorithm %<>%
    match.arg(
      choices = c("second","first", "2", "1"),
      several.ok = FALSE)
  if (generation_of_algorithm %in% c("first","1")) {
    generation_of_algorithm <- 1
  } else if (generation_of_algorithm %in% c("second", "2")) {
    generation_of_algorithm <- 2
  }
  if (verbose) {
    message("Data will be fitted to algorithm generation ",
            generation_of_algorithm, ".")
  }
  # Check force_through_origin
  assertive.types::assert_is_a_bool(force_through_origin)
  # Check split_by_experiment
  assertive.types::assert_is_a_bool(split_by_experiment)
  # check r_bar
  assertive.types::assert_is_a_number(r_bar)
  assertive.numbers::assert_all_are_greater_than(r_bar, 0)
  # check sigmar_bar
  assertive.types::assert_is_a_number(sigma_r_bar)
  assertive.numbers::assert_all_are_greater_than(sigma_r_bar, 0)
  # Return
  list(
    x = x,
    scale_to = scale_to,
    ppr_scale_factor = ppr_scale_factor,
    generation_of_algorithm = generation_of_algorithm,
    force_through_origin = force_through_origin,
    split_by_experiment = split_by_experiment,
    r_bar = r_bar,
    sigma_r_bar = sigma_r_bar) %>%
  return()
}