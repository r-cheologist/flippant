#' @title RatioInversion
#' @description A function to invert ratios
#' @details As a ratio of 2 indecates a two-fold enrichment and a ratio of 0.5 
#' the equivalent two-fold depletion, transformation of one into the other has 
#' to make use of the fact that in logarithmic space 0.5 and 2 are equivalent to
#' each other (-1 and 1, respectively for logarithm to the base of 2).
#' 
#' Accordingly this inversion function logarithmizes its numerical input, 
#' proceeds with inversion (multiplication with -1) and delogarithmizes the 
#' result.
#' @param x Object of numerical ratio values to be inverted.
#' @param base Numeric base of the logarithm to be used (defaulting to 2).
#' @return Returns an object equivalent to \code{x}, but with ratios inverted.
#' @author Johannes Graumann
#' @export
#' @keywords distribution methods manip
#' @examples
#' require(MDAA)
#' RatioInversion(0.5)
#' Ratioinversion(0.1,base=10)
#' RatioInversion(c(0.5,1,2))
#' RatioInversion(data.frame(A=c(0.5,1,2),B=c(0.5,1,2))
RatioInversion <- function(x,base=2){
  if(length(base) != 1){
    stop("\"base\" must be a single numeric value.")
  }
  returnX <- base^(log(x,base=base)*(-1))
  inputClass <- class(x)
  if(inputClass == "data.frame"){
    return(as.data.frame(returnX))
  } else {
    return(returnX)
  }
}