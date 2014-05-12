#' @title onAttach
#' @name onAttach
#' @description Function executed when attaching the package to a session
#' @details The function is executed upon attaching the package (e.g. via
#' \code{\link{library}}). It checks the parameters of the current R session and
#' produces a warning when the packages used diverge from the ones used for 
#' publication.
#' @author Johannes Graumann
#' @seealso \code{\link[base]{.onAttach}} \code{\link[pdapbase]{checkVersion}}
#' @import pdapbase 
#.onLoad <- function(...){
#  packageStartupMessage(versionCheck())
#}
.onAttach <- function(...){
  vC <- checkVersion(
    expectedR=package_version("3.1.0"),
    expectedPackages=list(
      assertive=package_version("0.1-8"),
      ggplot2=package_version("0.9.3.1"),
      pdapbase=package_version("0.0-1"),
      plyr=package_version("1.8.1"),
      robustbase=package_version("0.91-1"),
      scales=package_version("0.2.4")))
  if(!is.null(vC)){
    packageStartupMessage(vC)
  }
}
