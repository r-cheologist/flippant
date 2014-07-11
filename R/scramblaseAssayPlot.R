#' @title scramblaseAssayPlot
#' @aliases scrambalseAssayCalculations scramblaseAssayInputTemplate 
#' scramblaseAssayPlot scramblaseAssayStats scramblaseAssayTraces
#' @description Functions for the presentation and evaluaton of dithionite 
#' scramblase assays
#' @details The \code{\link{data.frame}} accepted by the majority of the 
#' functions (\code{x}) must have the following \bold{mandatory} columns:
#' \describe{
#'  \item{\code{Path}:}{Paths to existing and readable \code{ASCII} output files 
#'    of a fluorimeter. See \code{\link{parseFluorimeterOutput}} for details and
#'    supported formats.}
#'  \item{\code{Protein in Reconstitution (mg)}:}{Self-explanatory.}}
#' 
#' Further (\bold{facultative}) columns are:
#' \describe{
#'  \item{\code{Fluorescence Assay Vol. w/o DT (ul)}:}{Volume of the 
#'    fluorescence assay prior to addition of fluorescence-quenching ditihionite
#'    (defaulting to \code{2000}).}
#'  \item{\code{Fluorescence Assay Vol. with DT (ul)}:}{Volume of the 
#'    fluorescence assay after the addition of fluorescence-quenching 
#'    ditihionite (defaulting to \code{2040}).}
#'  \item{\code{Egg PC in Reconstitution (mmol)}:}{Self-explanatory. Defaulting to 
#'    \code{0.0045} (1 ml of a 1 mM solution.}
#'  \item{\code{Timepoint of Measurement (s)}:}{Timepoint (in seconds) used as 
#'    an anchor for the extraction of terminal fluorescence (defaulting to 
#'    \code{400}).}
#'  \item{\code{Experiment}:}{Identifier for any given experiment. Used for 
#'    \code{\link{facet_wrap}} during generation of \code{\link{ggplot}} output.}
#'  \item{\code{Experimental Series}:}{Identifier for a given series/graph (e.g.
#'    \code{Extract} and \code{Depleted Extract}). Used by \code{color} during 
#'    generation of \code{\link{ggplot}} output.}}
#'    
#' Based on MIKE PAPER data is processed as follows (the majority of the 
#' processing is split off into the internal function 
#' \code{flippant:::scramblaseAssayCalculations}):
#' \itemize{
#'  \item{Input is format checked and defaults are injected for facultative 
#'    parameters/columns as appropriate (see input \code{\link{data.frame}} 
#'    format above). The internal function 
#'    \code{flippant:::scramblaseAssayInputValidation} supplies this 
#'    functionality.}
#'  \item{Fluorescence spectra are parsed using 
#'    \code{\link{parseFluorimeterOutput}}. This includes automated 
#'    determination of when dithionite was added to the sample using 
#'    \pkg{wmtsa}-supplied methodology and resetting the acquisition time 
#'    accordingly (\code{0} henceforth corresponds to the time of addition).}
#'  \item{Pre-dithionite-addition \code{Baseline Fluorescence} is determined for
#'    each spectrum by averaging (\code{\link{median}}) over the 10 
#'    values preceding dithionite addition.}
#'  \item{Post-dithinonite-addition \code{Minimum Fluorescence} is determined 
#'    for each spectrum by averaging (\code{\link{median}}) over the last 10 
#'    smaller or equal to \code{400 s} (or \code{Timepoint of Measurement (s)},
#'    see above).}
#'  \item{The \code{Minimum Fluorescence} is volume-corrected based on 
#'    \code{Reaction Volume w/o DT (ul)} and \code{Reaction Volume with DT (ul)}
#'    (see above).}
#'  \item{For each spectrum/datapoint a measured \code{Fluorescence Reduction} 
#'    is calculated as \code{1 - Minimum Fluorescence/Baseline Fluorescence}.}
#'  \item{Data are \code{\link{split}} for parallel treatment using a combined 
#'    \code{Experimental Series}/\code{Experiment} identifier (see above).}
#'  \item{A \code{Relative Fluorescence Reduction} is calculated in comparison
#'    to the liposomes-only/no-protein control).}
#'  \item{A \code{Protein per Phospholipid (mg/mmol)} ratio (\code{PPR}) is 
#'    calculated.}
#'  \item{A probability for a liposome holding >= 1 scramblase molecule are 
#'    calculated using \code{(y - y0)/(ymax - y0)}, where \code{y} is the 
#'    \code{Relative Fluorescence Reduction} and \code{y0} is the 
#'    \code{Relative Fluorescence Reduction} in an experiment without addition 
#'    of protein extract. Depending on the \code{scaleTo} parameter, 
#'    \code{ymax} is either the maximal \code{Relative Fluorescence Reduction} 
#'    in the series (\code{scaleTo = "data"}) or derived from a 
#'    mono-exponential fit to the data (\code{scaleTo = "model"}). The latter 
#'    (default) is a precaution for the case where the protein/phospholipid
#'    titration did not reach the plateau of the saturation curve (yet).}
#'  \item{A monoexponential curve is fitted unsig \code{\link{nlxb}} to either
#'    \code{p(>=1) = b - c*(exp(-PPR/a))} (if \code{forceThroughOrigin = FALSE})
#'    or  \code{p(>=1) = b * (1 - exp(-PPR/a))}.}
#'  \item{Data \code{\link{split}} apart above are recombined and a 
#'    \code{\link{ggplot}} object is assembled with the following layers:
#'    \itemize{
#'      \item{Lines (\code{\link{geom_line}}) representing the monoexponential
#'        fit(s). \code{color} is used to differentiate 
#'        \code{Experimental Series}.}
#'      \item{Segments (\code{\link{geom_segment}}) representing the \code{PPR}
#'        at which the fit constant a is equal to \code{PPR}. This tau value 
#'        has the implication that at this PPR all vesicles on average have 1 
#'        scramblase and 63\% have 1 or more (i.e. are active). \code{color} is 
#'        used to differentiate \code{Experimental Series}.}
#'      \item{Points (\code{\link{geom_point}}) representing the corresponding 
#'        datapoints. \code{color} is used to differentiate 
#'        \code{Experimental Series}.}
#'      \item{Plots are finally \code{\link{facet_wrap}}ed by \code{Experiment} 
#'        and labels adjusted cosmetically.}}
#'  }}
#' @param path \code{\link{character}} object giving the path of an \bold{empty}
#' template for a spreadsheet that can provide \code{x}.
#' @param x \code{\link{data.frame}} as described in "Details".
#' @param scaleTo Defines the source of \code{ymax}, defaulting to 
#' \code{model}. See "Details".
#' @param forceThroughOrigin \code{\link{logical}} indicating whether to force 
#' the fitted curve(s) to penetrate the origin.
#' @param timeMin A single \code{\link{numeric}}. If given, 
#' \code{\link{scramblaseAssayTraces}} produces a time/x axis trimmed to
#' this value.
#' @param timeMax A single \code{\link{numeric}}. If given, 
#' \code{\link{scramblaseAssayTraces}} produces a time/x axis trimmed to
#' this value.
#' @return Returns a \code{\link{ggplot}} object.
#' @author Johannes Graumann
#' @references Menon, I., Huber, T., Sanyal, S., Banerjee, S., Barre, P., Canis, 
#' S., Warren, J.D., Hwa, J., Sakmar, T.P., and Menon, A.K. (2011). Opsin Is a 
#' Phospholipid Flippase. Current Biology 21, 149-153.
#' MIKE PAPER
#' @export
#' @seealso \code{\link{parseFluorimeterOutput}}
#' @import ggplot2
#' @importFrom plyr rbind.fill
#' @examples
#' #stop("Add citation to Mike's manuscript!")
#' #stop("Add example using actually published data.")
scramblaseAssayPlot <- function(x,scaleTo=c("model","data"),forceThroughOrigin=TRUE){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblaseAssayInputValidation(
    x =x ,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin)
  x <- validatedParams[["x"]]
  scaleTo <- validatedParams[["scaleTo"]]
  forceThroughOrigin <- validatedParams[["forceThroughOrigin"]]
# Processing --------------------------------------------------------------
  processedListFromX <- scramblaseAssayCalculations(
    x = x,
    scaleTo = scaleTo,
    forceThroughOrigin = forceThroughOrigin)

# Recombine the processed data --------------------------------------------
  x <- plyr::rbind.fill(lapply(processedListFromX,function(z){z$Raw}))
  fitResultsFromX <- plyr::rbind.fill(lapply(processedListFromX,function(z){z$Fit}))
  annotationsForX <- x[c("Experiment","Experimental Series","Fit Constant (a)")]
  #   annotationsForX <- x[c("Experiment","Experimental Series","PPR at P = 0.5")]
  names(annotationsForX) <- c("Experiment", "Experimental Series", "x1")
  annotationsForX$x2 <- annotationsForX$x1
  annotationsForX$y1 <- 1-exp(-1)
  annotationsForX$LineType <- 2
  #   annotationsForX$y1 <- 0.5
  annotationsForX$y2 <- -Inf
  annotationsForX <- annotationsForX[!duplicated(paste(annotationsForX$Experiment,annotationsForX$"Experimental Series")),]
  if(any(!is.na(annotationsForX$Experiment))){
    annotationListForX <- split(annotationsForX,annotationsForX$Experiment)
  } else {
    annotationListForX <- list(annotationsForX)
  }
  processedAnnotationListForX <- lapply(
    annotationListForX,
    function(x){
      data.frame(
        Experiment=unique(x$Experiment),
        "Experimental Series"=x$"Experimental Series"[1],
        x1=-Inf,
        x2=max(x$"x1",na.rm=TRUE),
        y1=unique(x$y1),
        y2=unique(x$y1),
        LineType=3,
        stringsAsFactors=FALSE)
    }
  )
  annotationsForX <- plyr::rbind.fill(annotationsForX,plyr::rbind.fill(processedAnnotationListForX))

# Assemble the (graphical) output -----------------------------------------
  # Groundwork
  plotOutput <- ggplot(
    data=x,
    aes_string(
      x="`Protein per Phospholipid (mg/mmol)`",
      y="`Probability >= 1 Scramblase in Vesicle`"))
  # Layering
  ## First layer: lines/curves representing the monoexponential fit
  if(any(!is.na(fitResultsFromX$"Experimental Series"))){
    plotOutput <- plotOutput +
      geom_line(data=fitResultsFromX,aes_string(color="`Experimental Series`"),fullRange=TRUE)
  } else {
    plotOutput <- plotOutput +
      geom_line(data=fitResultsFromX,fullRange=TRUE)
  }
  ## Second layer: annotations indicating PPR at tau
  if(any(!is.na(annotationsForX$"Experimental Series"))){
    plotOutput <- plotOutput +
      geom_segment(data=annotationsForX,aes_string(x="x1",xend="x2",y="y1",yend="y2",color="`Experimental Series`"),linetype=2)
  } else {
    plotOutput <- plotOutput +
      geom_segment(data=annotationsForX,aes_string(x="x1",xend="x2",y="y1",yend="y2"),linetype=2)
  }
  ## Third Layer: data points
  if(any(!is.na(x$"Experimental Series"))){
    plotOutput <- plotOutput +
      geom_point(aes_string(color="`Experimental Series`"))
  } else {
    plotOutput <- plotOutput +
      geom_point()
  }
  # Faceting by "Experiment"
  if(any(!is.na(x$"Experiment"))){
    plotOutput <- plotOutput + 
      facet_wrap(~Experiment)
  }
  # Prettifications
  plotOutput <- plotOutput +
    labs(
      x=expression(frac("Protein","Phospholipid")~~bgroup("(",frac("mg","mmol"),")")),
      y=expression(P~bgroup("(",frac("Scramblase","Liposome")>=1,")")),
      color="Experiment")
  # Return
  return(plotOutput)
}
