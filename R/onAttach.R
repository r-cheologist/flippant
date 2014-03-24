#' @title onAttach
#' @name onAttach
#' @description Function executed when attaching the package to a session
#' @details The function is executed upon attaching the package (e.g. via
#' \code{\link{library}}). It checks the parameters of the current R session and
#' produces a warning when the packages used diverge from the ones used for 
#' publication.
#' @author Johannes Graumann
#' @seealso \code{\link[base]{.onAttach}} \code{\link{version_check}}
#' @keywords internal misc
#.onLoad <- function(...){
#  packageStartupMessage(versionCheck())
#}
.onAttach <- function(...){
  vC <- version_check(
    expectedR=package_version("3.0.3"),
    expectedPackages=list(
      scales=package_version("0.2.3"),
      plyr=package_version("1.8.1"),
      ggplot2=package_version("0.9.3.1"),
      robustbase=package_version("0.90-2")))
  if(!is.null(vC)){
    packageStartupMessage(vC)
  }
}
