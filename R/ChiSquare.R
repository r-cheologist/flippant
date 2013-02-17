#' @title ChiSquare
#' @description Goodness-of-fit calculation as used in the original Protein 
#' Correlation Profiling Publication
#' @details The function iteratively calculates a \eqn{\chi^2} value for either 
#' all columns or all rows in a \code{\link{data.frame}} against asingle 
#' \code{\link{vector}}. See "Reference" for origin.
#' 
#' Note that \eqn{\chi^2} here is not directly related to 
#' \code{\link{chisq.test}} etc.
#' 
#' Further investigation revelas that \eqn{\chi^2} is actually related to the 
#' Mahalanobis distance under particular variance/covarinace matrix conditions
#' (see "Examples" section).
#' @param x Data frame of numerical values.
#' @param y Vector of numerical values rows/columns in x are correlated with.
#' Must be of length \code{nrow(x)} for \code{use.rows == TRUE} (default)
#' and length \code{ncol(x)} for \code{use.rows == FALSE}.
#' @param use.rows Whether to correlate all rows (default) or all columns 
#' iteratively with \code{y}.
#' @return Returns a numerical vector of \eqn{\chi^2} values with an additional
#' \code{attribute} representing the number of data points used (see "Examples").
#' @author Johannes Graumann
#' @references Andersen, J.S., Wilkinson, C.J., Mayor, T., Mortensen, P., 
#' Nigg, E.A., and Mann, M. (2003). Proteomic characterization of the human 
#' centrosome by protein correlation profiling. Nature 426, 570–574.
#' @export
#' @keywords methods manip
#' @examples
#' require(MDAA)
#' myData <- DataPrep(pGroups,na.replacement=NULL)
#' cs <- ChiSquare(myData$RelativeInverted,Relativate(flippaseActivity[["Av. Spec. Activity"]]))
#' cs
#' 
#' # Comparison with Mahalanobis distance
#' testA <- data.frame(rbind(Relativate(rnorm(100))))
#' testB <- Relativate(rnorm(100))
#' ChiSquare(testA,testB)
#' mahalanobis(x=testA,center=testB,cov=diag(length(testA)))
#' attr(cs,"DataPoints")
ChiSquare <- function(x,y,use.rows=TRUE){
  # Input checks
  ##############
  if(!is.data.frame(x)){
    stop("\"x\" must be a data.frame,")
  }
  rowsplit <- lapply(split(x,seq(nrow(x))),unlist)
  if(sum(sapply(rowsplit,is.numeric)) != length(rowsplit)){
    stop("\"x\" must be numeric.")
  }
  if(sum(x > 1,na.rm=TRUE) != 0){
    stop("Values in \"x\" don't seem to be relativated.")
  }
  if(!is.vector(y,mode="numeric")){
    stop("\"y\" must be a numeric vector.")
  }
  if((sum(y > 1,na.rm=TRUE) + sum(y < 0,na.rm=TRUE)) != 0){
    stop("Values in \"y\" don't seem to be relativated.")
  }
  if(!is.logical(use.rows)){
    stop("\"use.rows\" must be logical.")
  }
  if(length(use.rows) != 1){
    stop("\"use.rows\" must be of length 1.")
  }
  if(use.rows){
    if(length(y) != ncol(x)){
      stop("With \"use.rows=TRUE\" the length of \"y\" must equal \"nrow(x)\".")
    }
  } else {
    if(length(y) != nrow(x)){
      stop("With \"use.rows=FALSE\" the length of \"y\" must equal \"ncol(x)\".")
    }
    rm(rowsplit)
  }
  # Calculate Goodness-of-Fit
  ###########################
  # Rows or columns?
  if(use.rows){
    chiSqData <- rowsplit
    rm(rowsplit)
  } else {
    chiSqData <- lapply(seq(ncol(x)),function(z){return(x[[z]])})
  }
  # Calculation
  tmpChiSquare <- lapply(
    chiSqData,
    function(x){
      ## A consensus fractionation profile was obtained by averaging the 
      ## normalized peptide abundance profile for all peptides identifying the 
      ## known centrosomal proteins centrin 2, γ-tubulin and pericentrin.
      ## ‘χ2 values’ were calculated as the squared deviation of the normalized 
      ## profile for all peptides divided by the number of data points, ...
      pepProfile <- unlist(x)
      dps <- sum(!is.na(pepProfile))
      # Unsure whether to square sum or elements ... latter option correlates 
      # nicely with Spearman CC
      chiSquare <- sum((pepProfile - y)^2,na.rm=TRUE)/dps
      return(
        list(ChiSquare=chiSquare,DataPoints=dps))
    })
  # Assemble output
  chiSquare <- sapply(
    X=tmpChiSquare,
    FUN=function(x){x[["ChiSquare"]]},
    USE.NAMES=FALSE)
  attr(x=chiSquare,which="DataPoints") <- sapply(
    X=tmpChiSquare,
    FUN=function(x){x[["DataPoints"]]},
    USE.NAMES=FALSE)
  return(chiSquare)
}
  