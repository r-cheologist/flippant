#' @title onAttach
#' @name onAttach
#' @description Function executed when attaching the package to a session
#' @details The function is executed upon attaching the package (e.g. via
#' \code{\link{library}}). It checks the parameters of the current R session and
#' produces a warning when the packages used diverge from the ones used for 
#' publication.
#' @author Johannes Graumann
#' @seealso \code{\link[base]{.onAttach}} \code{\link{VersionCheck}}
#' @keywords internal misc
#.onLoad <- function(...){
#  packageStartupMessage(versionCheck())
#}
.onAttach <- function(...){
  vC <- VersionCheck(
    expectedR=package_version("3.0.2"),
    expectedPackages=list(
      gridExtra=package_version("0.9.1"),
      reshape2=package_version("1.2.2"),
      RCFPD=package_version("1.2.6.2"),
      XML=package_version("3.98-1.1"),
      scales=package_version("0.2.3"),
      seqinr=package_version("3.0.7"),
      R.utils=package_version("1.29.8"),
      plyr=package_version("1.8"),
      plotrix=package_version("3.5.3"),
      gsubfn=package_version("0.6-5"),
      png=package_version("0.1-7"),
      gtable=package_version("0.1.2"),
      ggplot2=package_version("0.9.3.1"),
      digest=package_version("0.6.4"),
      biomaRt=package_version("2.18.0"),
      robustbase=package_version("0.90-1")))
  if(!is.null(vC)){
    packageStartupMessage(vC)
  }
}
