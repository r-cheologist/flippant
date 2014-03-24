#' @title dithionite_flippase_assay_plot
#' @description A function that automates calculations necessary to interprete
#' dithionite flippase assays
#' @details The function accepts input in form of a \code{\link{data.frame}} 
#' with the following \bold{mandatory} columns:
#' \describe{
#'  \item{\code{Path}:}{Paths to existing and readable \code{ASCII} output files 
#'    of a fluorometer.}
#'  \item{\code{Protein in Reconstitution (mg)}:}{Self-explanatory.}}
#' 
#' Further (facultative) columns are:
#' \describe{
#'  \item{\code{Fluorescence Assay Vol. w/o DT (ul)}:}{Volume of the 
#'    fluorescence assay prior to addition of fluorescense-quenching ditihionite
#'    (defaulting to \code{2000}).}
#'  \item{\code{Fluorescence Assay Vol. with DT (ul)}:}{Volume of the 
#'    fluorescence assay after the addition of fluorescense-quenching 
#'    ditihionite (defaulting to \code{2040}).}
#'  \item{\code{Egg PC in Reconstitution (mmol)}:}{Self-explanatory. Defaulting to 
#'    \code{0.0045} (1 ml of a 1 mM solution.}
#'  \item{\code{Timepoint of Measurement (s)}:}{Timepoint used as an anchor for 
#'    the extraction of terminal fluorescense. 
#'    \code{\link{timepoint_of_measurement}} is used on all \code{Path}s if none 
#'    given.}
#'  \item{\code{Experiment}:}{Identifier for any given experiment. Used for 
#'    \code{\link{facet_wrap}} during generation of \code{\link{ggplot}} output.}
#'  \item{\code{Experimental Series}:}{Identifier for a given series/graph (e.g.
#'    \code{Extract} and \code{Depleted Extract}). Used by \code{color} during 
#'    generation of \code{\link{ggplot}} output.}}
#'    
#' Based on MIKE PAPER the function proceeds as follows (the majority of the 
#' data processing is split off into the internal function 
#' \code{\link{dithionite_flippase_assay_calculations}}):
#' \itemize{
#'  \item{Input is format checked and defaults are injected for facultative 
#'    parameters/columns as appropriate (see input \code{\link{data.frame}} 
#'    format above). The internal function 
#'    \code{\link{dithionite_flippase_assay_input_validation}} supplies this 
#'    functionality.}
#'  \item{Fluorescense spectra are parsed using \code{\link{parse_fluorometer_output}}.}
#'  \item{Pre-dithionite-addition \code{Baseline Fluorescense} is determined for
#'    each spectrum by averaging (\code{\link{median}}) over the first 10 
#'    values.}
#'  \item{Post-ditihinonite-addition \code{Minimum Fluorescense} is determined 
#'    for each spectrum by averaging (\code{\link{median}}) over the last 10 
#'    values common to all spectra.}
#'  \item{The \code{Minimum Fluorescense} is volume-corrected based on 
#'    \code{Reaction Volume w/o DT (ul)} and \code{Reaction Volume with DT (ul)}
#'    (see above).}
#'  \item{For each spectrum/datapoint a measured \code{Fluorescense Reduction} 
#'    is calculated as \code{1-Minimum Fluorescense/Baseline Fluorescense}.}
#'  \item{Data are \code{\link{split}} for parallel treatment using a combined 
#'    \code{Experimental Series}/\code{Experiment} identifier (see above).}
#'  \item{A \code{Relative Fluorescense Reduction} is calculated in comparison
#'    to the liposomes-only/no-protein control).}
#'  \item{A \code{Protein per Phospholipid (mg/mmol)} ratio (\code{PPR}) is 
#'    calculated.}
#'  \item{A probability for a liposome holding >= 1 flippase molecule are 
#'    calculated using \code{(y - y0)/(ymax - y0)}, where \code{y} is the 
#'    \code{Relative Fluorescense Reduction} and \code{y0} is the 
#'    \code{Relative Fluorescense Reduction} in an experiment without addition 
#'    of protein extract. Depending on the \code{scale_to} parameter, 
#'    \code{ymax} is either the maximal \code{Relative Fluorescense Reduction} 
#'    in the series (\code{scale_to = "data"}) or derived from a 
#'    mono-exponential fit to the data (\code{scale_to = "model"}). The latter 
#'    (default) is a precaution for the case where the protein/phospholipid
#'    titration did not reach the plateau of the saturation curve (yet).}
#'  \item{A monoexponential curve is fitted to \code{p(>=1) = 1 - exp(-PPR/a)} 
#'    using \code{\link{nlrob}}.}
#'  \item{Data \code{\link{split}} apart above are recombined and a 
#'    \code{\link{ggplot}} object is assembled with the following layers:
#'    \itemize{
#'      \item{Lines (\code{\link{geom_line}}) representing the monoexponential
#'        fit(s). \code{color} is used to differentiate \code{Experimental Series}.}
#'      \item{Segments (\code{\link{geom_segment}}) representing the \code{PPR}
#'        at which the fit constant a is equal to \code{PPR} and thus
#'        \code{p(>=1) = 1 - exp(-PPR/a) = 1 - exp(-1) ~ 0.63}. This tau value 
#'        has the implication that at this PPR all vesicles on average have 1 
#'        flippase and 63\% have 1 or more (i.e. are active). \code{color} is 
#'        used to differentiate \code{Experimental Series}.}
#'      \item{Points (\code{\link{geom_point}}) representing the corresponding 
#'        datapoints. \code{color} is used to differentiate \code{Experimental Series}.}
#'      \item{Plots are finally \code{\link{facet_wrap}}ed by \code{Experiment} and
#'        lables adjusted cosmetically.}}
#'  }}
#' @param x \code{\link{data.frame}} as described in "Details".
#' @param scale_to Defines the source of \code{ymax}, defaulting to 
#' \code{model}. See "Details".
#' @return Returns a \code{\link{ggplot}} object.
#' @author Johannes Graumann
#' @references Menon, I., Huber, T., Sanyal, S., Banerjee, S., Barre, P., Canis, 
#' S., Warren, J.D., Hwa, J., Sakmar, T.P., and Menon, A.K. (2011). Opsin Is a 
#' Phospholipid Flippase. Current Biology 21, 149-153.
#' MIKE PAPER
#' @export
#' @seealso \code{\link{parse_fluorometer_output}}, 
#' \code{\link{dithionite_flippase_assay_input_validation}},
#' \code{\link{dithionite_flippase_assay_calculations}}
#' @keywords methods manip
#' @import ggplot2
#' @import plyr
#' @import robustbase
#' @examples
#' stop("Add citation to Mike's manuscript!")
#' stop("Revisit workflow description.")
#' stop("Add example using actually published data.")
#' # Build input
#' x <- data.frame(
#'  Path = c(
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-plus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-plus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt"),
#'  "Extract Volume (ul)" = c(0,15,40,75,150,0,15,40,75,150,0,15,40,75,150,0,15,40,75,150),
#'  #     "Reaction Volume w/o DT (ul)" = rep(2000,4),
#'  "Reaction Volume with DT (ul)" = rep(2040,20),
#'  "Concentration Egg PC (mM)" = rep(4.5,20),
#'  "Extract Protein Concentration (mg/ml)" = c(rep(0.67,5),rep(1.26,5),rep(0.67,5),rep(1.26,5)),
#'  #     "Timepoint of Measurement (s)",
#'  "Experimental Series"=c(rep("Extract",5),rep("Depleted Extract",5),rep("Extract",5),rep("Depleted Extract",5)),
#'  Experiment=c(rep("Erg1, Replicate 1",10),rep("Erg1, Replicate 2",10)),
#'  check.names=FALSE,
#'  stringsAsFactors=FALSE)
#'  # Run function
#'  DithioniteFlippaseAssay(x)
dithionite_flippase_assay_plot <- function(x,scale_to=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validated_params <- dithionite_flippase_assay_input_validation(x=x,scale_to=scale_to)
  x <- validated_params[["x"]]
  scale_to <- validated_params[["scale_to"]]

# Processing --------------------------------------------------------------
  processed_list_from_x <- dithionite_flippase_assay_calculations(x=x,scale_to=scale_to)

# Recombine the processed data --------------------------------------------
  x <- rbind.fill(lapply(processed_list_from_x,function(z){z$Raw}))
  fit_results_from_x <- rbind.fill(lapply(processed_list_from_x,function(z){z$Fit}))
  annotations_for_x <- x[c("Experiment","Experimental Series","Fit Constant (a)")]
  #   annotations_for_x <- x[c("Experiment","Experimental Series","PPR at P = 0.5")]
  names(annotations_for_x) <- c("Experiment", "Experimental Series", "x1")
  annotations_for_x$x2 <- annotations_for_x$x1
  annotations_for_x$y1 <- 1-exp(-1)
  annotations_for_x$LineType <- 2
  #   annotations_for_x$y1 <- 0.5
  annotations_for_x$y2 <- -Inf
  annotations_for_x <- annotations_for_x[!duplicated(paste(annotations_for_x$Experiment,annotations_for_x$"Experimental Series")),]
  if(any(!is.na(annotations_for_x$Experiment))){
    annotation_list_for_x <- split(annotations_for_x,annotations_for_x$Experiment)
  } else {
    annotation_list_for_x <- list(annotations_for_x)
  }
  processed_annotation_list_for_x <- lapply(
    annotation_list_for_x,
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
  annotations_for_x <- rbind.fill(annotations_for_x,rbind.fill(processed_annotation_list_for_x))

# Assemble the (graphical) output -----------------------------------------
  # Groundwork
  plot_output <- ggplot(
    data=x,
    aes_string(
      x="`Protein per Phospholipid (mg/mmol)`",
      y="`Probability >= 1 Flippase in Vesicle`"))
  # Layering
  ## First layer: lines/curves representing the monoexponential fit
  if(any(!is.na(fit_results_from_x$"Experimental Series"))){
    plot_output <- plot_output +
      geom_line(data=fit_results_from_x,aes_string(color="`Experimental Series`"))
  } else {
    plot_output <- plot_output +
      geom_line(data=fit_results_from_x)
  }
  ## Second layer: annotations indicating PPR at tau
  if(any(!is.na(annotations_for_x$"Experimental Series"))){
    plot_output <- plot_output +
      geom_segment(data=annotations_for_x,aes_string(x="x1",xend="x2",y="y1",yend="y2",color="`Experimental Series`"),linetype=2)
  } else {
    plot_output <- plot_output +
      geom_segment(data=annotations_for_x,aes_string(x="x1",xend="x2",y="y1",yend="y2"),linetype=2)
  }
  ## Third Layer: data points
  if(any(!is.na(x$"Experimental Series"))){
    plot_output <- plot_output +
      geom_point(aes_string(color="`Experimental Series`"))
  } else {
    plot_output <- plot_output +
      geom_point()
  }
  # Faceting by "Experiment"
  if(any(!is.na(x$"Experiment"))){
    plot_output <- plot_output + 
      facet_wrap(~Experiment)
  }
  # Prettifications
  plot_output <- plot_output +
    labs(
      x=expression(frac("Protein","Phospholipid")~~bgroup("(",frac("mg","mmol"),")")),
      y=expression(P~bgroup("(",frac("Flippase","Liposome")>=1,")")),
      color="Experiment")
  # Return
  return(plot_output)
}
