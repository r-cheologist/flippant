#' @importFrom pdapbase checkVersion
.onAttach <- function(...){
  vC <- checkVersion(
    expectedR = package_version("3.1.1"),
    expectedPackages = list(
      assertive = package_version("0.1-8"),
      ggplot2 = package_version("1.0.0"),
      nlmrt = package_version("2013-9.24"),
      pdapbase = package_version("0.0-4"),
      plyr = package_version("1.8.1"),
      RcppRoll = package_version("0.1.0"),
      splus2R = package_version("1.2-0"),
      wmtsa = package_version("2.0-0")))
  if(!is.null(vC)){
    packageStartupMessage(vC)
  }
}
