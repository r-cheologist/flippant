#' @title Relativate
#' @description A function returning each value in a vector as relative to the 
#' vector's maximum value.
#' @details Divides each value in a numerical vector by the maximum value found 
#' in the vector, using \code{max(x,na.rm=TRUE)} to find the latter.
#' @param x Vector of numerical values.
#' @return Returns a vector of numeric values.
#' @author Johannes Graumann
#' @export
#' @keywords methods manip
#' @examples
#' require(MDAA)
#' Relativate(c(0.5,10,1,NA))
Relativate <- function(x){
  if(!is.vector(x)){
    stop("\"x\" must be of class vector.")
  }
  if(!is.numeric(x)){
    stop("\"x\" must be of class numeric.")
  }
  return(x/max(x,na.rm=TRUE))
}