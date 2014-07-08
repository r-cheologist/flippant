#' @importFrom pdapbase checkVersion
.onAttach <- function(...){
  vC <- checkVersion(
    expectedR=package_version("3.1.0"),
    expectedPackages=list(
      assertive=package_version("0.1-8"),
      ggplot2=package_version("1.0.0"),
      pdapbase=package_version("0.0-1"),
      plyr=package_version("1.8.1"),
      RcppRoll=package_version("0.1.0"),
      robustbase=package_version("0.91-1"),
      scales=package_version("0.2.4"),
      wmtsa=package_version("2.0-0")))
  if(!is.null(vC)){
    packageStartupMessage(vC)
  }
}
