#' @title DithioniteFlippaseAssayAnalysis
#' @description A function that automates calculations necessary to interprete
#' dithionite flippase assays
#' @details The function accepts input in form of a \code{\link{data.frame}} 
#' with the following \bold{mandatory} columns:
#' \describe{
#'  \item{\code{Path}:}{Paths to existing and readable \code{ASCII} output files 
#'    of a Photon QuantMaster fluorometer.}
#'  \item{\code{Extract Volume (ul)}:}{Volume of Triton X-100 extract used.}
#'  \item{\code{Extract Protein Concentration (mg/ml)}:}{Self-explanatory.}
#'  \item{\code{Experimental Series}:}{Identifier for a given series/graph (e.g.
#'    \code{Extract} and \code{Depleted Extract}).}}
#' 
#' Further (facultative) columns are:
#' \describe{
#'  \item{\code{Reaction Volume w/o DT (ul)}:}{Volume of the reaction prior to 
#'    addition of fluorescense-quenching ditihionite (defaulting to 
#'    \code{2000}).}
#'  \item{\code{Reaction Volume with DT (ul)}:}{Volume of the reaction after the
#'    addition of fluorescense-quenching ditihionite (defaulting to 
#'    \code{2040}).}
#'  \item{\code{Concentration Egg PC (mM)}:}{Self-explanatory. Defaulting to 
#'    \code{4.5}.}
#'  \item{\code{Timepoint of Measurement (s)}:}{Timepoint used as an anchor for 
#'    the extraction of terminal fluorescense. 
#'    \code{\link{TimepointOfMeasurement}} is used on all \code{Path}s if none 
#'    given.}
#'  \item{\code{Panel}:}{Used for \code{\link{facet_wrap}} during generation of
#'    \code{\link{ggplot}} output}.}
#'    
#' Based on MIKE PAPER the function proceeds as follows:
#' \itemize{
#'  \item{Input is format checked and defaults are injected for facultative 
#'    parameters/ columns as appropriate (see input \code{\link{data.frame}} 
#'    format above).}
#'  \item{Fluorescense spectra are parsed using 
#'    \code{\link{ParseFluorometerData2}}.}
#'  \item{Pre-dithionite-addition \code{Baseline Fluorescense} is determined for
#'    each spectrum by averaging (\code{\link{median}}) over the first 10 
#'    values.}
#'  \item{Post-ditihinonite-addition \code{Minimum Fluorescense} is determined 
#'    for each spectrum by averaging (\code{\link{median}}) over the last 10 
#'    values common to all spectra.}
#'  \item{The \code{Minimum Fluorescense} is volume-corrected based on 
#'    \code{Reaction Volume w/o DT (ul)} and \code{Reaction Volume with DT (ul)}
#'    (see above).}
#'  \item{For each spectrum/datapoint a \code{Relative Fluorescense Reduction} 
#'    is calculated as \code{1-Minimum Fluorescense/Baseline Fluorescense}.}
#'  \item{Data are \code{\link{split}} for parallel treatment using the 
#'    \code{Experimental Series} identifier (see above).}
#'  \item{p-values for a liposome holding >= 1 flippase molecule are calculated
#'    using \code{(y – y0)/(ymax – y0)}, where \code{y} is the 
#'    \code{Relative Fluorescense Reduction}, \code{y0} is the 
#'    \code{Relative Fluorescense Reduction} in an experiment without addition 
#'    of protein extract and \code{ymax} is the maximal
#'    \code{Relative Fluorescense Reduction} in the series.}
#'  \item{A \code{Protein per Phospholipid (mg/mmol)} ratio (\code{PPR}) is 
#'    calculated.}
#'  \item{A monoexponential curve is fitted to \code{p(≥1) = 1 – exp(-PPR/α)} 
#'    using \code{\link{nlrob}}.}
#'  \item{Data \code{\link{split}} apart above are recombined and a 
#'    \code{\link{ggplot}} object is assembled with the following layers:
#'    \itemize{
#'      \item{Lines (\code{\link{geom_line}}) representing the monoexponential
#'        fit(s). \code{color} is used to differentiate \code{Experimental Series}.}
#'      \item{Points (\code{\link{geom_point}}) representing the corresponding 
#'        datapoints. \code{color} is used to differentiate \code{Experimental Series}.}
#'      \item{Plots are finally \code{\link{facet_wrap}}ed by \code{Panel} and
#'        lables adjusted cosmetically.}}
#'  }}
#' @param x \code{\link{data.frame}} as described in "Details".
#' @return Returns a \code{\link{ggplot}} object.
#' @author Johannes Graumann
#' @references MIKE PAPER
#' @export
#' @seealso \code{\link{ParseQuantMasterData}}, \code{\link{TimepointOfMeasurement}}
#' @keywords methods manip
#' @import ggplot2
#' @import plyr
#' @import robustbase
#' @examples
#' stop("Add citation to Mike's manuscript!")
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
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt"),
#'  "Extract Volume (ul)" = c(0,15,40,75,150,0,15,40,75,150),
#'  #     "Reaction Volume w/o DT (ul)" = rep(2000,4),
#'  "Reaction Volume with DT (ul)" = rep(2040,10),
#'  "Concentration Egg PC (mM)" = rep(4.5,10),
#'  "Extract Protein Concentration (mg/ml)" = c(rep(0.67,5),rep(1.26,5)),
#'  #     "Timepoint of Measurement (s)",
#'  "Experimental Series"=c(rep("Extract",5),rep("Depleted Extract",5)),
#'  Panel=rep("Erg1, Replicate 1",10),
#'  check.names=FALSE,
#'  stringsAsFactors=FALSE)
#'  # Run function
#'  DithioniteFlippaseAssayAnalysis(x)
DithioniteFlippaseAssayAnalysis <- function(x)
{
  #######################
  # Check prerequisites #
  #######################
  # General DF characteristics
  ############################
  if(!is.data.frame(x)){
    stop("'x' must be of class 'data.frame'.")
  }
  if(nrow(x)==0){
    stop("'x' must have rows.")
  }
  if(any(is.na(x))){
    stop("'x' cannot contain 'NA'.")
  }
  # Required parameters
  #####################
  requirements <- list(
    Name = c(
      "Path",
      "Extract Volume (ul)",
      "Extract Protein Concentration (mg/ml)",
      "Experimental Series"),
    Class = c(
      "character",
      "numeric",
      "numeric",
      "character"))
  if(!all( requirements$Name %in% names(x))){
    stop(
      "'x' must hold at least the following columns: '",
      paste(requirements,collapse="', '"),
      "'.")
  }
  if(!identical(
    unname(vapply(x[requirements$Name],class,c(A="A"))),
    requirements$Class)){
    stop(
      "Required columns '",
      paste(requirements$Name,collapse="', '"),
      "' must be of classes '",
      paste(requirements$Class,collapse="', '"),
      "'.")
  }
  # Check paths
  if(!all(file.exists(x$Path))){
    stop("All entries in column 'Path' must refer to existing files.")
  }
  if(any(file.access(x$Path,mode=4) == -1)){
    stop("All entries in column 'Path' must refer to existing files.")
  }
  # Facultative parameters
  ########################
  facultatives <- list(
    Name = c(
      "Reaction Volume w/o DT (ul)",
      "Reaction Volume with DT (ul)",
      "Concentration Egg PC (mM)",
      "Timepoint of Measurement (s)",
      "Panel"),
    Class = c(
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "character"),
    Default = c(
      2000,
      2040,
      4.5,
      NA,
      NA))
  missing <- which(!(facultatives$Name %in% names(x)))
  if(length(missing) != 0){
    for (y in missing){
      if(facultatives$Name[y] == "Timepoint of Measurement (s)"){
        warning(
          "Providing missing column '",
          facultatives$Name[y],
          "' from spectra ('Path').")
        addOn <- TimepointOfMeasurement(x$Path)
      } else {
        warning(
          "Providing missing column '",
          facultatives$Name[y],
          "' from defaults (",
          facultatives$Default[y],
          "). Make sure this is correct.")
        addOn <- facultatives$Default[y]
      }
      tmpX <- cbind(
        x,
        rep(x=addOn,times=nrow(x)),
        stringsAsFactors=FALSE)
      names(tmpX)[ncol(tmpX)] <- facultatives$Name[y]
      x <- tmpX
    }
  }
  if(!identical(
    unname(vapply(x[facultatives$Name],class,c(A="A"))),
    facultatives$Class)){
    stop(
      "Facultative columns '",
      paste(facultatives$Name,collapse="', '"),
      "' must be of classes '",
      paste(facultatives$Class,collapse="', '"),
      "'.")
  }
  # Check "Timepoint of Measurement (s)" consistency
  if(length(unique(x$"Timepoint of Measurement (s)")) != 1){
    stop("Column 'Timepoint of Measurement (s)' contains multiple values. Exiting.")
  }
  ##############
  # Processing #
  ##############
  # Parsing spectra
  #################
  tmpData <- lapply(
    x$Path,
    function(y){ParseQuantMasterData(SpecFile=y)})
  # What spectral time windows to extract?
  minAT <- unique(vapply(tmpData,function(y){y$"Minimal Acquisition Time (s)"},1))
  if(length(minAT) != 1){
    stop("Minimum acquisition times are not identical - aborting.")
  }
  maxAT <- vapply(tmpData,function(y){y$"Maximal Acquisition Time (s)"},1)
  if(any(maxAT < x$"Timepoint of Measurement (s)")){
    stop("'Timepoint of Measurement (s)' is larger than the shortest spectrum 
         acquisition.")
  } else {
    maxAT <- unique(x$"Timepoint of Measurement (s)")
  }
  # Average over first 10 values for activity baseline
  x$"Baseline Fluorescense" <- vapply(
    tmpData,
    function(z){
      tmpFrom <- min(which(z$Data$"Time (s)" >= minAT))
      baselineSS <- seq(from=tmpFrom,to=tmpFrom+9)
      return(median(z$Data$"Fluorescense Intensity"[baselineSS],na.rm=TRUE))
    },
    1)
  # Average over last 10 values (in common time range) for activity
  x$"Minimum Fluorescense" <- vapply(
    tmpData,
    function(z){
      tmpTo <- max(which(z$Data$"Time (s)" <= maxAT))
      activitySS <- seq(from=tmpTo-9,to=tmpTo)
      return(median(z$Data$"Fluorescense Intensity"[activitySS],na.rm=TRUE))
    },
    1)
  # Apply volume correction factors as needed
  correctionFactor <- x$"Reaction Volume with DT (ul)"/x$"Reaction Volume w/o DT (ul)"
  x$"Minimum Fluorescense, Volume Corrected" <- x$"Minimum Fluorescense" * correctionFactor
  # Calculate relative activity reduction
  x$"Relative Fluorescense Reduction" <- 1-x$"Minimum Fluorescense, Volume Corrected"/x$"Baseline Fluorescense"
  # Split by Experiment
  #####################
  xiList <- split(x,x$"Experimental Series")
  # Generate PPR vs. p>=1Flippase/Liposome Data
  #############################################
  xoList <- lapply(
    xiList,
    function(z){
      # Ensure that there's a data point with liposomes ONLY as a unique 
      # reference point
      liposomesOnlyIndex <- which(z$"Extract Volume (ul)" == 0)
      if(length(liposomesOnlyIndex) == 0){
        stop("Experimental series '",unique(y$"Experimental Series"),"' does not
          have the required liposomes-ONLY ('Extract Volume (ul)' of '0') data 
          point.")
      }
      if(length(liposomesOnlyIndex) > 1){
        stop("Experimental series '",unique(y$"Experimental Series"),"' has more
          than one liposomes-ONLY ('Extract Volume (ul)' of '0') data point.")
      }  
      ## End-point fluorescence reduction data from flippase activity assays were 
      ## obtained for proteoliposomes generated over a range of PPR values.
      ## The data were transformed according to the formula
      ## p(≥1 flippase) = (y – yo)/(ymax – yo),
      ## where yo is the percent reduction obtained with liposomes, ymax is the 
      ## maximum percentage reduction observed and p is the probability that a 
      ## particular vesicle in the ensemble is ‘flippase-active’, i.e it possesses 
      ## at least one flippase.
      # Calculate p>=1Flippase/Liposome
      y <- z$"Relative Fluorescense Reduction"
      y0 <- z$"Relative Fluorescense Reduction"[liposomesOnlyIndex]
      ymax <- max(z$"Relative Fluorescense Reduction",na.rm=TRUE)
      z$"Pvalue >= 1 Flippase in Vesicle" <- (y-y0)/(ymax-y0)
      ## The dependence of p(≥1 flippase) on PPR was analyzed as follows.
      ## Definitions:
      ##   f, number of flippases used for reconstitution
      ##   v, number of vesicles(
      ##   m, number of flippases per vesicle (=f/v)
      ##   PPR, mg protein per mmol phospholipid
      # Calculate PPR
      z$"Protein per Phospholipid (mg/mmol)" <- (z$"Extract Protein Concentration (mg/ml)"*z$"Extract Volume (ul)")/z$"Concentration Egg PC (mM)"
      # DEBUG: plot(z$"Protein per Phospholipid (mg/mmol)",z$"Pvalue >= 1 Flippase in Vesicle")
      ## To calculate p(≥1 flippase) as a function of the PPR, we assume that 
      ## reconstitution of opsin/rhodopsin molecules into vesicles occurs 
      ## independently and that the vesicles are identical and may have more than 
      ## one flippase. The probability that a flippase will be reconstituted into a 
      ## particular vesicle is therefore 1/v. On reconstituting f flippases into an 
      ## ensemble of v vesicles, the probability p(k) that a particular vesicle 
      ## contains k flippases is given by the binomial formula:
      ##   
      ##   p(k) = C(f,k)(1/v)k(1-1/v)f-k(
      ## 
      ## Because f and v are both large, it is convenient to use the Poisson 
      ## approximation:
      ##   
      ##   p(k) = (mk/k!)e-m, where m = f/v is the average number of flippases per vesicle
      ## 
      ## The probability of a particular vesicle having no flippases is p(0) = e-m; 
      ## therefore, the probability that a vesicle has one or more flippases, i.e is
      ## active in the flippase assay, is
      ## 
      ## p(≥1) = 1-p(0) = 1 - e-m 
      ## 
      ## The average number of flippases per vesicle, m, is proportional to PPR and 
      ## can be written as m = PPR/α, where α is a constant with units of mg/mmol. 
      ## Thus,
      ## 
      ## p(≥1) = 1 - e-m = 1 – exp(-PPR/α)
      ## 
      ## The mono-exponential fit constant for a graph of p(≥1) vs PPR is α mg/mmol; 
      ## at this PPR value, m = 1 and ~63% of the vesicles in the population possess 
      ## ≥1 flippase.
      ## If we assume that reconstitution of an opsin/rhodopsin monomer (MW 41.7 
      ## kDa) into a 200 nm diameter vesicle confers flippase activity to that 
      ## vesicle, then α = 0.122 mg/mmol (note: 1 mmol of phospholipids yields ~1.75
      ## x 1015 200 nm-diameter vesicles (Mimms et al, 1981); 1 mg of opsin or 
      ## rhodopsin corresponds to 1.44 x 1016 molecules; m = f/v = 1 corresponds to
      ## 1.75 x 1015 opsin/rhodopsin molecules per mmol phospholipid or 0.12 
      ## mg/mmol). If opsin/rhodopsin molecules exist as preformed dimers and 
      ## reconstitution of such dimers into a vesicle confers flippase activity to 
      ## that vesicle, then α = 0.24 mmol/mg.
      # Fit a monoexponential curve to the data
      tmpData <- data.frame(
        x=z$"Protein per Phospholipid (mg/mmol)",
        y=z$"Pvalue >= 1 Flippase in Vesicle")
      library(robustbase)
      Rmod <- nlrob(y ~ 1-exp(-x/a),data=tmpData, start = list(a=1))
      z$"Fit Constant (a)" <- Rmod$coefficients
      output <- list(Raw=z)
      # Generate data to plot the results of the fit
      predictX <- data.frame(
        x=seq(
          from=min(z$"Protein per Phospholipid (mg/mmol)",na.rm=TRUE),
          to=max(z$"Protein per Phospholipid (mg/mmol)",na.rm=TRUE),
          length.out=200))
      predictedY <- predict(
          object=Rmod,
          newdata=predictX)
      output$Fit <- data.frame(
        "Protein per Phospholipid (mg/mmol)"=predictX$x,
        "Pvalue >= 1 Flippase in Vesicle"=predictedY,
        "Experimental Series"=unique(z$"Experimental Series"),
        "Panel"=unique(z$"Panel"),
        check.names=FALSE,
        stringsAsFactors=FALSE)
      # Return
      return(output)
    })
  # Recombine the processed data
  ##############################
  library(plyr)
  x <- rbind.fill(lapply(xoList,function(z){z$Raw}))
  fitX <- rbind.fill(lapply(xoList,function(z){z$Fit}))
  ###############################
  # Assemble (graphical) output #
  ###############################
  library(ggplot2)
  # Groundwork (color-separation of "Experimental Series")
  tmpPlot <- ggplot(
    data=x,
    aes_string(
      x="`Protein per Phospholipid (mg/mmol)`",
      y="`Pvalue >= 1 Flippase in Vesicle`",
      color="`Experimental Series`"))
  # Layering
  tmpPlot <- tmpPlot +
    # First layer: lines/curves representing the monoexponential fit
    geom_line(data=fitX) +
    # Second Layer: data points
    geom_point() + 
    # Faceting by "Panel"
    facet_wrap(~Panel) +
    # Prettifications
    labs(
      x=expression(frac("Protein","Phospholipid")~~bgroup("(",frac("mg","mmol"),")")),
      y=expression(p~bgroup("(",frac("Flippase","Liposome")>=1,")")),
      color="Experiment")
  # Return
  return(tmpPlot)
}