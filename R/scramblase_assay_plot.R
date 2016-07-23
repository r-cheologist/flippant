#' @title scramblase_assay_plot
#' @aliases scramblaseAssayCalculations scramblase_assay_input_template 
#' scramblase_assay_plot scramblase_assay_stats scramblase_assay_traces
#' @description Functions for the presentation and evaluaton of dithionite 
#' scramblase assays
#' @details The \code{\link{data.frame}} accepted by the majority of the 
#' functions a an \code{R} object or path to a corresponding file (\code{x}) 
#' must have the following \bold{mandatory} columns:
#' \describe{
#'  \item{\code{Path}:}{Paths to existing and readable \code{ASCII} output files 
#'    of a fluorimeter. See \code{\link{parse_fluorimeter_output}} for details and
#'    supported formats.}
#'  \item{\code{Protein Reconstituted (mg)}:}{Self-explanatory.}}
#' 
#' Further (\bold{facultative}) columns are:
#' \describe{
#'  \item{\code{Fluorescence Assay Vol. w/o DT (ul)}:}{Volume of the 
#'    fluorescence assay prior to addition of fluorescence-quenching ditihionite
#'    (defaulting to \code{2000}).}
#'  \item{\code{Fluorescence Assay Vol. with DT (ul)}:}{Volume of the 
#'    fluorescence assay after the addition of fluorescence-quenching 
#'    ditihionite (defaulting to \code{2040}).}
#'  \item{\code{Lipid in Reconstitution (mmol)}:}{Self-explanatory. For the 
#'    standard phospholipid experiment defaulting to \code{0.0045} (1 ml of a 
#'    4.5 mM solution).}
#'  \item{\code{Timepoint of Measurement (s)}:}{The time to determine terminal 
#'    fluorescence, calculated from the point when dithionite is added, in
#'    seconds, defaulting to \code{400}).}
#'  \item{\code{Experiment}:}{Identifier for any given experiment. Used for 
#'    \code{\link{facet_wrap}} during generation of \code{\link{ggplot}} output.
#'    All data with one \code{Experiment} identifier ends up on one plot/facet.}
#'  \item{\code{Experimental Series}:}{Identifier for a given series/graph (e.g.
#'    \code{Extract} and \code{Depleted Extract}). Used by \code{color} during 
#'    generation of \code{\link{ggplot}} output to differentiate lines in the
#'    same plot/facet.}}
#'    
#' Based on Goren et al. (2014) and Ploier et al. (2016) data is processed as
#' follows (the majority of the processing is split off into the internal
#' function \code{scramblase_assay_calculations}):
#' \itemize{
#'  \item{Input is format checked and defaults are injected for facultative 
#'    parameters/columns as appropriate (see input \code{\link{data.frame}} 
#'    format above). The internal function 
#'    \code{scramblase_assay_input_validation} supplies this 
#'    functionality.}
#'  \item{Fluorescence spectra are parsed using 
#'    \code{\link{parse_fluorimeter_output}}. This includes automated 
#'    determination of when dithionite was added to the sample using 
#'    \pkg{wmtsa}-supplied methodology and resetting the acquisition time 
#'    accordingly (\code{0} henceforth corresponds to the time of addition).}
#'  \item{Pre-dithionite-addition \code{Baseline Fluorescence} is determined for
#'    each spectrum by averaging (\code{\link{median}}) over the 10 
#'    values preceding dithionite addition.}
#'  \item{Post-dithinonite-addition \code{Minimum Fluorescence} is determined 
#'    for each spectrum by averaging (\code{\link{median}}) over the last ten
#'    datapoints \eqn{\leq 400\,\mbox{s}}{\ge 400 s} (or 
#'    \code{Timepoint of Measurement (s)}, see above).}
#'  \item{The \code{Minimum Fluorescence} is volume-corrected based on 
#'    \code{Reaction Volume w/o DT (ul)} and \code{Reaction Volume with DT (ul)}
#'    (see above).}
#'  \item{For each spectrum/datapoint a measured \code{Fluorescence Reduction} 
#'    is calculated as 
#'    \deqn{1 - \left(\frac{\mbox{\small Minimum Fluorescence}}{\mbox{\small Baseline Fluorescence}}\right)}{1 - (Minimum Fluorescence/Baseline Fluorescence)}}
#'  \item{A \code{Relative Fluorescence Reduction} is calculated in comparison
#'    to the liposomes-only/no-protein control).}
#'  \item{A \code{Protein per Phospholipid (mg/mmol)} ratio (\code{PPR}) is 
#'    calculated. If \code{ppr_scale_factor} is not \code{NULL}, the value is
#'    scaled (divided) by that value to account for liposomes that remain
#'    inaccessible to reconstitution with scramblase molecules.}
#'  \item{Depending on \code{split_by_experiment}, data are \code{\link{split}} 
#'    for parallel treatment using either \code{Experimental Series}
#'    (\code{split_by_experiment = TRUE}) or a combined
#'     \code{Experimental Series}/\code{Experiment}
#'     (\code{split_by_experiment = FALSE}) identifier (see above).}
#'  \item{A probability for a liposome holding \eqn{\geq 1}{\ge 1} scramblase 
#'    molecules is calculated using 
#'    \deqn{\frac{y-y_0}{y_{\mbox{\scriptsize max}}-y_0}}{(y - y0)/(ymax - y0)}
#'    where \eqn{y} is the \code{Relative Fluorescence Reduction} and \eqn{y_0}{y0}
#'    is the \code{Relative Fluorescence Reduction} in an experiment without
#'    addition of protein extract. Depending on the \code{scale_to} parameter, 
#'    \eqn{y_{\mbox{\scriptsize max}}}{ymax} is either the maximal \code{Relative Fluorescence Reduction} 
#'    in the series (\code{scale_to = "data"}) or derived from a 
#'    mono-exponential fit to the data (\code{scale_to = "model"}). The latter 
#'    (default) is a precaution for the case where the protein/phospholipid
#'    titration did not reach the plateau of the saturation curve.}
#'  \item{A monoexponential curve is fitted using \code{\link{nlsLM}}.
#'  
#'    If \code{generation_of_algorithm} is \code{1}, the underlying formula is
#'    derived from Goren et al. (2014) and data is fitted to either
#'    \deqn{p(\geq 1)=b\cdot(1-e^{-\frac{\mbox{\tiny PPR}}{a}})}{p(\ge 1) = b * (1 - exp(-PPR/a))} 
#'    (if \code{force_through_origin = TRUE}; default) or
#'    \deqn{p(\geq 1)=b-c\cdot e^{-\frac{\mbox{\tiny PPR}}{a}}}{p(\ge 1) = b - c*exp(-PPR/a)}
#'    (if \code{force_through_origin = FALSE}). The latter implies more degrees
#'    of freedom and occasionaly results in better fits to experimental data.
#'    Mechanistic implication, however, are unclear.
#'    
#'    If \code{generation_of_algorithm} is \code{2} (default), the more
#'    elaborate model put forth in Ploier et al. (2016) is employed, using
#'    either
#'    \deqn{p(\geq 1)=b\cdot(\frac{1}{\sqrt{1+\sigma^2\cdot a \cdot x}})\cdot exp(\frac{-\bar{r}^2\cdot a \cdot x}{1+\sigma^2\cdot a\cdot x})}{p(\ge 1) = b * (1 - (1/\sqrt{1 + \sigma^2 * a * x}) * exp((-\bar{r}^2 * a * x)/(1 + \sigma^2 * a * x))}
#'    (if \code{force_through_origin = TRUE}; default) or
#'    \deqn{p(\geq 1)=b-c\cdot(\frac{1}{\sqrt{1+\sigma^2\cdot a \cdot x}})\cdot exp(\frac{-\bar{r}^2\cdot a \cdot x}{1+\sigma^2\cdot a\cdot x})}{p(\ge 1) = b - c (1 - (1/\sqrt{1 + \sigma^2 * a * x}) * exp((-\bar{r}^2 * a * x)/(1 + \sigma^2 * a * x))}
#'    (if \code{force_through_origin = FALSE}).}
#'  \item{Data \code{\link{split}} apart above are recombined and a 
#'    \code{\link{ggplot}} object is assembled with the following layers:
#'    \itemize{
#'      \item{Lines (\code{\link{geom_line}}) representing the monoexponential
#'        fit(s). \code{color} is used to differentiate 
#'        \code{Experimental Series}.}
#'      \item{If \code{generation_of_algorithm} is \code{1}, segments
#'       (\code{\link{geom_segment}}) representing the \code{PPR}
#'        at which the fit constant \eqn{a} is equal to \code{PPR}. This 
#'        \eqn{\tau}{tau} value has the implication that at this \code{PPR} all 
#'        vesicles on average have one scramblase and 63\% have one or more 
#'        (i.e. are active). \code{color} is used to differentiate 
#'        \code{Experimental Series}. Where \code{gerneration_of_algorithm} is
#'        \code{2}, interpretation of \eqn{a} is less obvious and this layer is
#'        thus ommited in the plot.}
#'      \item{Points (\code{\link{geom_point}}) representing the corresponding 
#'        datapoints. \code{color} is used to differentiate 
#'        \code{Experimental Series}.}
#'      \item{Plots are finally \code{\link{facet_wrap}}ed by \code{Experiment} 
#'        (if \code{split_by_experiment = TRUE}) and labels adjusted
#'        cosmetically.}}
#'  }}
#' @param path \code{\link{character}} object giving the path of an \bold{empty}
#' template for a spreadsheet that can provide \code{x}.
#' @param x \code{\link{data.frame}} or path to a tab delimited file 
#' representing it (see "Details").
#' @param ppr_scale_factor \code{\link{numeric}} object providing a scale factor
#' to adjust internally calculated \code{Protein per Phospholipid (mg/mmol)}
#' ratios (\code{PPR}; see "Details").
#' @param scale_to Defines the source of \code{ymax}, defaulting to 
#' \code{model}. See "Details".
#' @param force_through_origin \code{\link{logical}} indicating whether to force 
#' the fitted curve(s) to penetrate the origin (defaulting to \code{TRUE}). See
#' "Details".
#' @param time_min_sec A single \code{\link{numeric}}. If given, 
#' \code{\link{scramblase_assay_traces}} produces a time/x axis trimmed to
#' this value (in seconds).
#' @param time_max_sec A single \code{\link{numeric}}. If given, 
#' \code{\link{scramblase_assay_traces}} produces a time/x axis trimmed to
#' this value (in seconds).
#' @param adjust A single \code{\link{logical}}, indicating whether (default) or 
#' not spectral traces to be plotted are algorithmically aligned at the time
#' point of dithionite addition.
#' @param generation_of_algorithm Either \code{2} or \code{1}
#' (\code{\link{numeric}}; defaulting to \code{2}). See "Details".
#' @param split_by_experiment A single \code{\link{logical}}, indicating whether or
#' not calculations and plots will treat experimental series from different
#' experiments separately (\code{TRUE}, default) or whether data from all
#' experiments included is used for a single calculation/plot per experimental
#' series (\code{FALSE}). While the former emphasizes reproducibility, the
#' latter likely produces a more reliable fit.
#' @param r_bar A \code{\link{numeric}}, representing the average radius of the
#' liposomes used in the assay. Only used in \code{generatio_of_algorithm = 2}
#' and defaulting to \code{88} (see Ploier et al. 2016 for details).
#' @param sigma_r_bar A \code{\link{numeric}}, representing the standard
#' deviationaverage of the radius distribution  of the liposomes used in the
#' assay. Only used in \code{generatio_of_algorithm = 2} and defaulting to
#' \code{28} (see Ploier et al. 2016 for details).
#' @return \code{scrambalse_assay_traces} and \code{scramblase_assay_plot} return 
#' \code{\link{ggplot}} objects representing the raw fluorescence traces and a
#' complete PPR plot, respectively. \code{scrambalseAssayInputTemplate} 
#' generates a tab-delimited \code{ASCII} file in the file system and does not
#' provide further output. \code{scrambalseAssayStats} assembles (and prints) 
#' assay statistics as a \code{\link{data.frame}}.
#' @author Johannes Graumann
#' @references Menon, I., Huber, T., Sanyal, S., Banerjee, S., Barre, P., Canis, 
#' S., Warren, J.D., Hwa, J., Sakmar, T.P., and Menon, A.K. (2011). Opsin Is a 
#' Phospholipid Flippase. Current Biology 21, 149-153.
#' 
#' Goren, M.A., Morizumi, T., Menon, I., Joseph, J.S., Dittman, J.S., 
#' Cherezov, V., Stevens, R.C., Ernst, O.P., and Menon, A.K. (2014). 
#' Constitutive phospholipid scramblase activity of a G Protein-coupled 
#' receptor. Nat Commun 5, 5115.
#' 
#' Ploier, B., Caro, L.N., Morizumi, T., Pandey, K., Pearring, J.N.,
#' Goren, M.A., Finnemann, S.C., Graumann, J., Arshavsky, V.Y., Dittman, J.S.,
#' Ernst, O.P., Menon, A.K. (2016). Dimerization deficiency of enigmatic
#' retinitis pigmentosa-linked rhodopsin mutants. Nat Commun. In Press.
#' @export
#' @seealso \code{\link{parse_fluorimeter_output}} \code{\link{nlsLM}}
#' @import ggplot2
#' @examples
#' library(magrittr)
#' library(ggplot2)
#' # Extract example data
#' analysis_dir <- file.path(tempdir(), "flippant-case-study")
#' extract_case_study_data(analysis_dir)
#' template_file <- file.path(analysis_dir, "inputTable.txt")
#' # Plot the spectral traces
#' scramblase_assay_traces(
#'   template_file,
#'   time_max_sec = 350)
#' # Plot the PPR plot(s) faceting by experiment
#' scramblase_assay_plot(template_file)
#' # Generate tabular results
#' scramblase_assay_stats(template_file)
#' # Plot the PPR plot(s) forgoing faceting by experiment
#' scramblase_assay_plot(template_file, split_by_experiment = FALSE)
#' # Generate tabular results
#' scramblase_assay_stats(template_file, split_by_experiment = FALSE)
scramblase_assay_plot <- function(
  x,
  scale_to = c("model","data"),
  ppr_scale_factor = 0.65,
  force_through_origin = TRUE,
  generation_of_algorithm = c(2, 1),
  split_by_experiment = TRUE,
  r_bar = 88,
  sigma_r_bar = 28){
  UseMethod("scramblase_assay_plot",x)
}
#' @export
scramblase_assay_plot.data.frame <- function(x, ...){
  base_function_scramblase_assay_plot(x, ...)
}

#' @export
scramblase_assay_plot.character <- function(x, ...){
  parsedInputFile <- read_scramblase_input_file(x)
  withr::with_dir(
    dirname(x),
    base_function_scramblase_assay_plot(x=parsedInputFile, ...)
  )
}
base_function_scramblase_assay_plot <- function(
  x,
  scale_to = c("model","data"),
  ppr_scale_factor = 0.65,
  force_through_origin = TRUE,
  generation_of_algorithm = c(2, 1),
  split_by_experiment = TRUE,
  r_bar = 88,
  sigma_r_bar = 28){
# Check Prerequisites -----------------------------------------------------
  validatedParams <- scramblase_assay_input_validation(
    x = x ,
    scale_to = scale_to,
    ppr_scale_factor = ppr_scale_factor,
    force_through_origin = force_through_origin,
    generation_of_algorithm = generation_of_algorithm,
    split_by_experiment = split_by_experiment,
    r_bar = r_bar,
    sigma_r_bar = sigma_r_bar)
  x <- validatedParams[["x"]]
  scale_to <- validatedParams[["scale_to"]]
  ppr_scale_factor <- validatedParams[["ppr_scale_factor"]]
  force_through_origin <- validatedParams[["force_through_origin"]]
  generation_of_algorithm <- validatedParams[["generation_of_algorithm"]]
  split_by_experiment <- validatedParams[["split_by_experiment"]]
  
# Processing --------------------------------------------------------------
  processedListFromX <- scramblase_assay_calculations(
    x = x,
    scale_to = scale_to,
    ppr_scale_factor = ppr_scale_factor,
    force_through_origin = force_through_origin,
    generation_of_algorithm = generation_of_algorithm,
    split_by_experiment = split_by_experiment,
    r_bar = r_bar,
    sigma_r_bar = sigma_r_bar)

# Recombine the processed data --------------------------------------------
  x <- plyr::rbind.fill(
    lapply(
      processedListFromX,
      function(z){z[["Raw"]]}))
  fitResultsFromX <- plyr::rbind.fill(
    lapply(
      processedListFromX,
      function(z){z[["Fit"]]}))
  annotationsForX <- x[c("Experiment","Experimental Series","Fit Constant (a)")]
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
  annotationsForX <- plyr::rbind.fill(
    annotationsForX,
    plyr::rbind.fill(processedAnnotationListForX))

# Assemble the (graphical) output -----------------------------------------
  # Groundwork
  plotOutput <- ggplot(
    data=x,
    aes_string(
      x="`Protein per Phospholipid (mg/mmol)`",
      y="`Probability >= 1 Scramblase in Vesicle`"))
  # Layering
  ## First layer: lines/curves representing the monoexponential fit
  plotOutput <- plotOutput +
    geom_line(
      data=fitResultsFromX,
      aes_string(color = get_color_var(fitResultsFromX)))
  ## Second layer: annotations indicating PPR at tau
  if(generation_of_algorithm == 1){
      plotOutput <- plotOutput +
        geom_segment(
          data = annotationsForX,
          aes_string(
            x = "x1",
            xend = "x2",
            y = "y1",
            yend = "y2",
            color = get_color_var(annotationsForX)),
          linetype = 2)
  }
  ## Third Layer: data points
  plotOutput <- plotOutput +
    geom_point(aes_string(color=get_color_var(x)))
    
  # Faceting by "Experiment"
  if(any(!is.na(x$"Experiment"))){
    if(split_by_experiment){
      plotOutput <- plotOutput + 
        facet_wrap(~Experiment)
    }
  }
  # Prettifications
  x_label <- if(is.null(ppr_scale_factor)){
    expression(frac("Protein","Phospholipid")~~bgroup("(",frac("mg","mmol"),")"))
  } else {
    expression(frac("Protein","Phospholipid")^"adj."~~bgroup("(",frac("mg","mmol"),")"))
  }
  plotOutput <- plotOutput +
    labs(
      x=x_label,
      y=expression(P~bgroup("(",frac("Scramblase","Liposome")>=1,")")),
      color="Experiment")
  
  # Return
  return(plotOutput)
}

get_color_var <- function(data) {
  if(any(!is.na(data$"Experimental Series"))) 
  {
    "`Experimental Series`"
  }
}
