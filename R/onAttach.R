#' @importFrom pdapbase checkVersion
.onAttach <- function(...){
  vC <- checkVersion(
    expectedR = package_version("3.2.2"),
    expectedPackages = list(
      assertive = package_version("0.3-0"),
      ggplot2 = package_version("1.0.1"),
      nlmrt = package_version("2013-9.25"),
      pdapbase = package_version("0.5.2"),
      plyr = package_version("1.8.3"),
      RcppRoll = package_version("0.2.2"),
      splus2R = package_version("1.2-0"),
      wmtsa = package_version("2.0-0")))
  if(!is.null(vC)){
    packageStartupMessage(vC)
  }
}
