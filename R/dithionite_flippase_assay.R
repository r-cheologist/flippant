#' @title dithionite_flippase_assay
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
#'    \code{\link{TimepointOfMeasurement}} is used on all \code{Path}s if none 
#'    given.}
#'  \item{\code{Experiment}:}{Identifier for any given experiment. Used for 
#'    \code{\link{facet_wrap}} during generation of \code{\link{ggplot}} output.}
#'  \item{\code{Experimental Series}:}{Identifier for a given series/graph (e.g.
#'    \code{Extract} and \code{Depleted Extract}). Used by \code{color} during 
#'    generation of \code{\link{ggplot}} output.}}
#'    
#' Based on MIKE PAPER the function proceeds as follows:
#' \itemize{
#'  \item{Input is format checked and defaults are injected for facultative 
#'    parameters/columns as appropriate (see input \code{\link{data.frame}} 
#'    format above).}
#'  \item{Fluorescense spectra are parsed
#'   
#'     using 
#'  #'  \item{\code{Fluorometer}:}{
#'    Fluorometer producing the file at \code{Path}.
#'    Currently the function can parse \code{ASCII} output as produced by: 
#'    Photon QuantMaster, Perkin Elmer LS55. The corresponding legal 
#'    \code{\link{character}} values in the column are \code{QuantMaster} and 
#'    \code{LS55}. The default is \code{QuantMaster}.
#'  }

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
#'  \item{Data are \code{\link{split}} for parallel treatment using a combined 
#'    \code{Experimental Series}/\code{Experiment} identifier (see above).}
#'  \item{p-values for a liposome holding >= 1 flippase molecule are calculated
#'    using \code{(y - y0)/(ymax - y0)}, where \code{y} is the 
#'    \code{Relative Fluorescense Reduction}, \code{y0} is the 
#'    \code{Relative Fluorescense Reduction} in an experiment without addition 
#'    of protein extract and \code{ymax} is the maximal
#'    \code{Relative Fluorescense Reduction} in the series.}
#'  \item{A \code{Protein per Phospholipid (mg/mmol)} ratio (\code{PPR}) is 
#'    calculated.}
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
#' @return Returns a \code{\link{ggplot}} object.
#' @author Johannes Graumann
#' @references Menon, I., Huber, T., Sanyal, S., Banerjee, S., Barré, P., Canis, 
#' S., Warren, J.D., Hwa, J., Sakmar, T.P., and Menon, A.K. (2011). Opsin Is a 
#' Phospholipid Flippase. Current Biology 21, 149–153.
#' MIKE PAPER
#' @export
#' @seealso \code{\link{ParseQuantMasterData}}, \code{\link{ParseLS55Data}}, 
#' \code{\link{TimepointOfMeasurement}}
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
dithionite_flippase_assay <- function(x){
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
  required_columns_in_x <- list(
    Name = c(
      "Path",
      "Protein in Reconstitution (mg)"),
    Class = c(
      "character",
      "numeric"))
  if(!all( required_columns_in_x$Name %in% names(x))){
    stop(
      "'x' must hold at least the following columns: '",
      paste(required_columns_in_x,collapse="', '"),
      "'.")
  }
  if(!identical(
    unname(vapply(x[required_columns_in_x$Name],class,c(A="A"))),
    required_columns_in_x$Class)){
    stop(
      "Required columns '",
      paste(required_columns_in_x$Name,collapse="', '"),
      "' must be of classes '",
      paste(required_columns_in_x$Class,collapse="', '"),
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
  facultative_columns_in_x <- list(
    Name = c(
      "Fluorescence Assay Vol. w/o DT (ul)",
      "Fluorescence Assay Vol. with DT (ul)",
      "Egg PC in Reconstitution (mmol)",
      "Timepoint of Measurement (s)",
      "Experiment",
      "Experimental Series"),
    Class = c(
      "numeric",
      "numeric",
      "numeric",
      "numeric",
      "character",
      "character"),
    Default = list(
      2000,
      2040,
      0.0045,
      NA,
      NA_character_,
      NA_character_))
  missing_facultative_columns_in_x <- which(!(facultative_columns_in_x$Name %in% names(x)))
  if(length(missing_facultative_columns_in_x) != 0){
    for (y in missing_facultative_columns_in_x){
      if(facultative_columns_in_x$Name[y] == "Timepoint of Measurement (s)"){
        warning(
          "Providing missing column '",
          facultative_columns_in_x$Name[y],
          "' from spectra ('Path').")
        to_be_added_on <- TimepointOfMeasurement(x$Path)
      } else if(facultative_columns_in_x$Name[y] %in% c("Experiment","Experimental Series")) {
        to_be_added_on <- facultative_columns_in_x$Default[[y]]
      } else {
        warning(
          "Providing missing column '",
          facultative_columns_in_x$Name[y],
          "' from defaults (",
          facultative_columns_in_x$Default[[y]],
          "). Make sure this is correct.")
        to_be_added_on <- facultative_columns_in_x$Default[[y]]
      }
      output <- cbind(
        x,
        rep(x=to_be_added_on,times=nrow(x)),
        stringsAsFactors=FALSE)
      names(output)[ncol(output)] <- facultative_columns_in_x$Name[y]
      x <- output
    }
  }
  if(!identical(
    unname(vapply(x[facultative_columns_in_x$Name],class,c(A="A"))),
    facultative_columns_in_x$Class)){
    stop(
      "Facultative columns '",
      paste(facultative_columns_in_x$Name,collapse="', '"),
      "' must be of classes '",
      paste(facultative_columns_in_x$Class,collapse="', '"),
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
  spectral_data <- lapply(
    x$Path,
    parse_fluorometer_output)
  # What spectral time windows to extract?
  minAT <- unique(vapply(spectral_data,function(y){y$Min.Acquisition.Time.in.sec},1))
  if(length(minAT) != 1){
    stop("Minimum acquisition times are not identical - aborting.")
  }
  maxAT <- vapply(spectral_data,function(y){y$Max.Acquisition.Time.in.sec},1)
  if(any(maxAT < x$"Timepoint of Measurement (s)")){
    stop("'Timepoint of Measurement (s)' is larger than the shortest spectrum 
         acquisition.")
  } else {
    maxAT <- unique(x$"Timepoint of Measurement (s)")
  }
  # Average over first 10 values for activity baseline
  x$"Baseline Fluorescense" <- vapply(
    spectral_data,
    function(z){
      start_index_for_averaging <- min(which(z$Data$Time.in.sec >= minAT))
      indexes_for_averaging <- seq(
        from=start_index_for_averaging,
        to=start_index_for_averaging+9)
      return(median(z$Data$Fluorescense.Intensity[indexes_for_averaging],na.rm=TRUE))
    },
    1)
  # Average over last 10 values (in common time range) for activity
  x$"Minimum Fluorescense" <- vapply(
    spectral_data,
    function(z){
      stop_index_for_averaging <- max(which(z$Data$Time.in.sec <= maxAT))
      indexes_for_averaging <- seq(
        from=stop_index_for_averaging-9,
        to=stop_index_for_averaging)
      return(median(z$Data$Fluorescense.Intensity[indexes_for_averaging],na.rm=TRUE))
    },
    1)
  # Apply volume correction factors as needed
  volume_correction_factor <- x$"Fluorescence Assay Vol. with DT (ul)"/x$"Fluorescence Assay Vol. w/o DT (ul)"
  x$"Minimum Fluorescense, Volume Corrected" <- x$"Minimum Fluorescense" * volume_correction_factor
  # Calculate relative activity reduction
  x$"Relative Fluorescense Reduction" <- 1-x$"Minimum Fluorescense, Volume Corrected"/x$"Baseline Fluorescense"
  # Split by Experiment
  #####################
  x$CombinedId <- paste(x$"Experimental Series",x$"Experiment",sep="_")
  input_list_from_x <- split(x,x$CombinedId)
  # Generate PPR vs. p>=1Flippase/Liposome Data
  #############################################
  processed_list_from_x <- lapply(
    input_list_from_x,
    function(z){
      # Ensure that there's a data point with liposomes ONLY as a unique 
      # reference point
      index_of_liposomes_only_data <- which(z$"Protein in Reconstitution (mg)" == 0)
      if(length(index_of_liposomes_only_data) == 0){
        stop("Experimental series '",unique(y$"Experimental Series"),"' does not
          have the required liposomes-ONLY ('Extract Volume (ul)' of '0') data 
          point.")
      }
      if(length(index_of_liposomes_only_data) > 1){
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
      y0 <- z$"Relative Fluorescense Reduction"[index_of_liposomes_only_data]
      ymax <- max(z$"Relative Fluorescense Reduction",na.rm=TRUE)
      z$"Pvalue >= 1 Flippase in Vesicle" <- (y-y0)/(ymax-y0)
      ## The dependence of p(≥1 flippase) on PPR was analyzed as follows.
      ## Definitions:
      ##   f, number of flippases used for reconstitution
      ##   v, number of vesicles(
      ##   m, number of flippases per vesicle (=f/v)
      ##   PPR, mg protein per mmol phospholipid
      # Calculate PPR
      z$"Protein per Phospholipid (mg/mmol)" <- (z$"Protein in Reconstitution (mg)"/z$"Egg PC in Reconstitution (mmol)")
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
      subset_for_fit <- data.frame(
        x=z$"Protein per Phospholipid (mg/mmol)",
        y=z$"Pvalue >= 1 Flippase in Vesicle")
      Rmod <- nlrob(y ~ 1-exp(-x/a),data=subset_for_fit, start = list(a=1), maxit=40)
      z$"Fit Constant (a)" <- Rmod$coefficients
      z$"PPR at P = 0.5" <- -Rmod$coefficients * log(1-0.5)
      output <- list(Raw=z)
      # Generate data to plot the results of the fit
      x_predicted_from_fit <- data.frame(
        x=seq(
          from=min(z$"Protein per Phospholipid (mg/mmol)",na.rm=TRUE),
          to=max(z$"Protein per Phospholipid (mg/mmol)",na.rm=TRUE),
          length.out=200))
      y_predicted_from_fit <- predict(
          object=Rmod,
          newdata=x_predicted_from_fit)
      output$Fit <- data.frame(
        "Protein per Phospholipid (mg/mmol)"=x_predicted_from_fit$x,
        "Pvalue >= 1 Flippase in Vesicle"=y_predicted_from_fit,
        "Experimental Series"=unique(z$"Experimental Series"),
        "Experiment"=unique(z$"Experiment"),
        check.names=FALSE,
        stringsAsFactors=FALSE)
      # Return
      return(output)
    })
  # Recombine the processed data
  ##############################
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
  ###############################
  # Assemble (graphical) output #
  ###############################
  # Groundwork
  plot_output <- ggplot(
    data=x,
    aes_string(
      x="`Protein per Phospholipid (mg/mmol)`",
      y="`Pvalue >= 1 Flippase in Vesicle`"))
  # Layering
  ## First layer: lines/curves representing the monoexponential fit
  if(any(!is.na(fit_results_from_x$"Experimental Series"))){
    plot_output <- plot_output +
      geom_line(data=fit_results_from_x,aes_string(color="`Experimental Series`"))
  } else {
    plot_output <- plot_output +
      geom_line(data=fit_results_from_x)
  }
  ## Second layer: annotations indicating PPR at p=0.5
  if(any(!is.na(annotations_for_x$"Experimental Series"))){
    plot_output <- plot_output +
      geom_segment(data=annotations_for_x,aes_string(x="x1",xend="x2",y="y1",yend="y2",color="`Experimental Series`"),linetype=2)
  } else {
    plot_output <- plot_output +
      geom_segment(data=annotations_for_x,aes(x=x1,xend=x2,y=y1,yend=y2),linetype=2)
  }
  ## Third Layer: data points
  if(any(!is.na(x$"Experimental Series"))){
    plot_output <- plot_output +
      geom_point(aes_string(color="`Experimental Series`"))
  } else {
    plot_output <- plot_output +
      geom_point()
  }
  ## Faceting by "Experiment"
  if(any(!is.na(x$"Experiment"))){
    plot_output <- plot_output + 
      facet_wrap(~Experiment)
  }
  # Prettifications
  plot_output <- plot_output +
    labs(
      x=expression(frac("Protein","Phospholipid")~~bgroup("(",frac("mg","mmol"),")")),
      y=expression(p~bgroup("(",frac("Flippase","Liposome")>=1,")")),
      color="Experiment")
  # Return
  return(plot_output)
}