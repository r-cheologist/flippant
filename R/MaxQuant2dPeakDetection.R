#' @title MaxQuant2dPeakDetection
#' @aliases MaxQuant2dPeakDetection 
#' @description 2D peak detection for mass spectra as implement in MaxQuant.
#' @details Citing from Cox and Mann (2008; see \cite{References} below):
#' \sQuote{In each MS scan, peaks are detected in a conventional two-dimensional
#' (2D) way by first searching for local maxima of the intensity as a function of 
#' m/z. The lower and upper limits of the m/z interval for a 2D peak ... are then 
#' determined by moving from the maximum position to smaller and larger m/z values,
#' until either the intensity has dropped to zero ..., or a local intensity minimum 
#' has been reached ... . This straightforward approach of peak detection without 
#' any de-convolution, smoothing or de-noising is sufficient for MS data generated 
#' by modern high precision mass spectrometers such as LTQ FT or Orbitrap.
#' }
#'
#' From the figures in Cox and Mann (2008) it seems apparent that the peak starts and
#' stops determined as above are trimmed and the default \code{trim} replicates that
#' behavior.
#' @param x A \code{\link{numeric}} vector of intensities in the order of their
#' profile contribution (e.g. the \code{y} values in the order of their \code{x}
#' partners).
#' @param trim Whether to return peaks with the terminal points removed or not 
#' (defaulting to \code{TRUE}).
#' @param return.half Whether or not to return not fully qualified peaks at the
#' sequence extremes (defaulting to \code{TRUE}). 
#' @return Returns a \code{\link{list}} of \code{\link{numeric}} vectors of 
#' length 2, indicating the starting and stopping index of detected peaks. 
#' Returns \code{\link{NULL}} if no peak is detected.
#' @author Johannes Graumann
#' @references
#' \cite{Cox, J., and Mann, M. (2008). MaxQuant enables high peptide 
#' identification rates, individualized p.p.b.-range mass accuracies 
#' and proteome-wide protein quantification. Nat. Biotechnol 26, 
#' 1367-1372.}
#' @export
#' @examples
#' # Some sample data
#' y <- c(9,8,7,5,4,1,1,2,1,1,3,4,5,6,7,5,4,3,2,1,1,3,4,5,6,7,8)
#' x <- seq(length(y))
#' # Accumulate the output of a couple of argument-variations
#' peakList <- list(
#'  MaxQuant2dPeakDetection(y,trim=FALSE,return.half=TRUE),
#'  MaxQuant2dPeakDetection(y,trim=FALSE,return.half=FALSE),
#'  MaxQuant2dPeakDetection(y,trim=TRUE,return.half=FALSE))
#' # Showcase the output
#' for(peaks in peakList){
#'  plot(x,y)
#'  rect(
#'    xleft=sapply(peaks,function(z){x[z[1]]}),
#'    xright=sapply(peaks,function(z){x[z[2]]}),
#'    ybottom=0,
#'    ytop=10,
#'    col=rainbow(length(peaks)),
#'    density=10)
#'  readline("Hit <ENTER> to proceed ...")
#' }
MaxQuant2dPeakDetection <- function(x,trim=TRUE,return.half=TRUE){
  #######################
  # Check Prerequisites #
  #######################
  if(!is.numeric(x)){
    stop("\"x\" must be numeric.")
  }
  if(!is.logical(trim) | length(trim) != 1){
    stop("\"trim\" must be a single logical.")
  }
  if(!is.logical(return.half) | length(return.half) != 1){
    stop("\"return.half\" must be a single logical.")
  }
  #################
  # Function body #
  #################
  # Determine peak boundaries
  ###########################
  # ..., peaks are detected in a conventional two-dimensional
  # (2D) way by first searching for local maxima of the intensity as a function of 
  # m/z. The lower and upper limits of the m/z interval for a 2D peak ... are then 
  # determined by moving from the maximum position to smaller and larger m/z values,
  # until either the intensity has dropped to zero ..., or a local intensity minimum 
  # has been reached ... . This straightforward approach of peak detection without 
  # any de-convolution, smoothing or de-noising is sufficient for MS data generated 
  # by modern high precision mass spectrometers such as LTQ FT or Orbitrap.
  y1 <- c(NA,x,NA)
  y2 <- c(x,NA,NA)
  y3 <- c(NA,NA,x)
  # Where does a peak start?
  lThanNext <- y1 < y2
  leThanLast <- y1 <= y3
  peakStart <- which((lThanNext + leThanLast) == 2) - 1
  # Is there a maximum?
  hasMax <- sum(na.omit((lThanNext + leThanLast) == 0)) > 0
  # Where does a peak stop?
  leThanNext <- y1 <= y2
  lThanLast <- y1 < y3
  peakStop <- which((leThanNext + lThanLast) == 2) - 1
  # Deal with not fully qualified peaks at the sequence extremes:
  if(suppressWarnings(min(peakStop,na.rm=TRUE)) < suppressWarnings(min(peakStart,na.rm=TRUE))){
    if(return.half){
      peakStart <- c(1,peakStart)
    } else {
      peakStop <- peakStop[seq(from=2,to=length(peakStop))]
    }}
  if(suppressWarnings(max(peakStart,na.rm=TRUE)) > suppressWarnings(max(peakStop,na.rm=TRUE))){
    if(return.half){
      peakStop <- c(peakStop,length(x))
    } else {
      peakStart <- peakStart[seq(length(peakStart)-1)]
    }
  }
  # Deal with zero return
  if(length(peakStart) == 0 & length(peakStop) == 0){
    # Take care of single big peak case
    if(hasMax){
      peakStart <- c(1,peakStart)
      peakStop <- c(peakStop,length(x))
    } else {
      return(NULL)
    }
  }
  # Split the data into the peaks
  ###############################
  peaks <- lapply(
    X = seq(length(peakStart)),
    FUN = function(x){
      return(c(peakStart[x],peakStop[x]))
    })
  # Peak trimming
  ###############
  if(trim){
    peaks <- lapply(
      peaks,
      function(x){
        # Points In Peak
        pip <- x[2]-x[1]+1
        # Only trim peaks with pip > 2
        if(pip < 3){return(x)}
        return(x + c(1,-1))
      })
  }
  # Assemble output #
  ###################
  return(peaks)
}