#' @title version_check
#' @description Function executed when attaching the package to a session
#' @details The function checks the parameters of the current R session and
#' produces a warning character string when the R/package versions used
#' diverge from the expectations.
#' 
#' \code{expected_Packages} must be a named list with the name providing the 
#' package name and the associated object being a \code{\link{package_version}} 
#' object.
#' @param expectedR The R version used to process/plot data for publication as
#' a \code{\link{package_version}} object.
#' @param expectedPackages A named list of \code{\link{package_version}} objects 
#' to be tested for (see "Details" for more).
#' @return In case a version difference is found, a character string is returned, 
#' verbalising the results. Otherwise \code{invisible(NULL) is returned.}
#' @author Johannes Graumann
#' @export
#' @seealso \code{\link[base]{.onAttach}}
#' @keywords misc
#' @examples
#' # Silent - Version matches
#' version_check(expectedR=getRversion())
#' \dontrun{
#' # String produced as nonsense version requested
#' version_check(expectedR=package_version(""88.88.88"))
#' }
#' @importFrom utils packageVersion
version_check <-  function(
  expectedR = getRversion(),
  expectedPackages = NULL
  ){
  # Input Checking
  if(!inherits(expectedR,"package_version")){
    stop("\"expectedR\" must be a \"package_version\" objects.")
  }
  if(length(expectedR) != 1){
    stop("\"expectedR\" must be of lenght 1.")
  }
  if(!is.null(expectedPackages)){
    if(!is.list(expectedPackages)){
      stop("\"expectedPackages\" must be a list.")
    }
    if(is.null(names(expectedPackages))){
      stop("\"expectedPackages\" must be a named list.")
    }
    if(sum(sapply(expectedPackages,function(x){inherits(x,"package_version")})) != length(expectedPackages)){
      stop("All entries in \"expectedPackages\" must be of class \"package_version\".")
    }
  }
  # Outputvariable
  output <- NULL
  # Check R version
  currR <- getRversion()
  if(currR != expectedR){
    locString <-paste(
      "WARNING: The works referring to this package are based on the usage of R version \'",
      expectedR,
      "\'. The current version \'",
      currR,
      "\' deviates from that and results/output may thus differ.\n",sep="")
    output <- paste(output,locString,collapse="\n")
  }
  # Check package versions
  for(name in names(expectedPackages)){
    currPkg <- packageVersion(name)
    expPkg <- expectedPackages[[name]]
    if(currPkg != expPkg){
      locString <- paste(
        "WARNING: The works referring to this package are based on the usage of \'",
        name,
        "\' version \'",
        expPkg,
        "\'. The current version \'",
        currPkg,
        "\' deviates from that and results/output may thus differ.\n",sep="")
      output <- paste(output,locString,collapse="\n")
    }
  }
  #Return Output
  if(is.null(output)){
    return(invisible(output))
  } else {
    return(output)
  }
}
