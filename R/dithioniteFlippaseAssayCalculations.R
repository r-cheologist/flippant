#' @rdname dithioniteFlippaseAssayPlot
dithioniteFlippaseAssayCalculations <- function(x,scale_to){
# Parsing spectra ---------------------------------------------------------
  spectral_data <- lapply(
    x$Path,
    parseFluorometerOutput)

# Read out data -----------------------------------------------------------  
  # What spectral time windows to extract?
  min_acquisition_time <- unique(vapply(spectral_data,function(y){y$Min.Acquisition.Time.in.sec},1))
  if(length(min_acquisition_time) != 1){
    stop("Minimum acquisition times are not identical - aborting.")
  }
  max_acquisition_time <- vapply(spectral_data,function(y){y$Max.Acquisition.Time.in.sec},1)
  if(any(max_acquisition_time < x$"Timepoint of Measurement (s)")){
    stop("'Timepoint of Measurement (s)' is larger than the shortest spectrum 
         acquisition. Aborting.")
  } else {
    max_acquisition_time <- unique(x$"Timepoint of Measurement (s)")
  }
  # Average over first 10 values for activity baseline
  x$"Baseline Fluorescense" <- vapply(
    spectral_data,
    function(z){
      start_index_for_averaging <- min(which(z$Data$Time.in.sec >= min_acquisition_time))
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
      stop_index_for_averaging <- max(which(z$Data$Time.in.sec <= max_acquisition_time))
      indexes_for_averaging <- seq(
        from=stop_index_for_averaging-9,
        to=stop_index_for_averaging)
      return(median(z$Data$Fluorescense.Intensity[indexes_for_averaging],na.rm=TRUE))
    },
    1)
  # Apply volume correction factors as needed
  volume_correction_factor <- x$"Fluorescence Assay Vol. with DT (ul)"/x$"Fluorescence Assay Vol. w/o DT (ul)"
  x$"Minimum Fluorescense, Volume Corrected" <- x$"Minimum Fluorescense" * volume_correction_factor
  ## Calculate activity reduction
  x$"Fluorescense Reduction" <- 1-x$"Minimum Fluorescense, Volume Corrected"/x$"Baseline Fluorescense"

# Generate PPR v.s P(>=1 Flippase/Liposome) -------------------------------
  # Split by Experiment
  x$CombinedId <- paste(x$"Experimental Series",x$"Experiment",sep="_")
  input_list_from_x <- split(x,x$CombinedId)

  # Calculations
  processed_list_from_x <- lapply(
    input_list_from_x,
    function(z){
      ## Ensure that there's a data point with liposomes ONLY as a unique 
      ## reference point
      index_of_liposomes_only_data <- which(z$"Protein in Reconstitution (mg)" == 0)
      if(length(index_of_liposomes_only_data) == 0){
        stop("Experimental series '",unique(y["Experimental Series"]),"' does 
          not have the required liposomes-ONLY ('Extract Volume (ul)' of '0') 
          data point. Aborting.")
      }
      if(length(index_of_liposomes_only_data) > 1){
        stop("Experimental series '",unique(y["Experimental Series"]),"' has 
          more than one liposomes-ONLY ('Extract Volume (ul)' of '0') data 
          point.")
      }  
      ##> End-point fluorescence reduction data from flippase activity assays 
      ##> were obtained for proteoliposomes generated over a range of PPR 
      ##> values.
      ##> The data were transformed according to the formula
      ##> p(≥1 flippase) = (y – yo)/(ymax – yo),
      ##> where yo is the percent reduction obtained with liposomes, ymax is the 
      ##> maximum percentage reduction observed and p is the probability that a 
      ##> particular vesicle in the ensemble is ‘flippase-active’, i.e it 
      ##> possesses at least one flippase.
      ## Calcualte the relative fluorescence reduction
      z$"Relative Fluorescense Reduction" <- z$"Fluorescense Reduction" - z$"Fluorescense Reduction"[index_of_liposomes_only_data]
      ## Calculate PPR
      z$"Protein per Phospholipid (mg/mmol)" <- (z$"Protein in Reconstitution (mg)"/z$"Egg PC in Reconstitution (mmol)")
      ## Calculate p>=1Flippase/Liposome
      y <- z$"Relative Fluorescense Reduction"
      y0 <- z$"Relative Fluorescense Reduction"[index_of_liposomes_only_data]
      if(scale_to == "model"){
        subset_for_fit <- data.frame(
          x=z$"Protein per Phospholipid (mg/mmol)",
          y=z$"Relative Fluorescense Reduction")
        Rmod <- nlrob(y ~ b-exp(-x/a),data=subset_for_fit, start = list(a=1,b=max(z$"Relative Fluorescense Reduction",na.rm=TRUE)), maxit=40)
        ymax <- max(z$"Relative Fluorescense Reduction",na.rm=TRUE) * Rmod$coefficient["b"]
      } else {
        ymax <- max(z$"Relative Fluorescense Reduction",na.rm=TRUE)
      }
      z$"Probability >= 1 Flippase in Vesicle" <- (y-y0)/(ymax-y0)
      ##> The dependence of p(≥1 flippase) on PPR was analyzed as follows.
      ##> Definitions:
      ##>   f, number of flippases used for reconstitution
      ##>   v, number of vesicles(
      ##>   m, number of flippases per vesicle (=f/v)
      ##>   PPR, mg protein per mmol phospholipid
      ##> To calculate p(≥1 flippase) as a function of the PPR, we assume that 
      ##> reconstitution of opsin/rhodopsin molecules into vesicles occurs 
      ##> independently and that the vesicles are identical and may have more
      ##> than one flippase. The probability that a flippase will be 
      ##> reconstituted into a particular vesicle is therefore 1/v. On 
      ##> reconstituting f flippases into an ensemble of v vesicles, the 
      ##> probability p(k) that a particular vesicle contains k flippases is 
      ##> given by the binomial formula:
      ##>   
      ##>   p(k) = C(f,k)(1/v)k(1-1/v)f-k(
      ##> 
      ##> Because f and v are both large, it is convenient to use the Poisson 
      ##> approximation:
      ##>   
      ##>   p(k) = (mk/k!)e-m, where m = f/v is the average number of flippases
      ##> per vesicle
      ##>
      ##> The probability of a particular vesicle having no flippases is 
      ##> p(0) = e-m; therefore, the probability that a vesicle has one or more
      ##> flippases, i.e isactive in the flippase assay, is
      ##> 
      ##> p(≥1) = 1-p(0) = 1 - e-m 
      ##> 
      ##> The average number of flippases per vesicle, m, is proportional to PPR
      ##> and can be written as m = PPR/α, where α is a constant with units of
      ##> mg/mmol. 
      ##> Thus,
      ##> 
      ##> p(≥1) = 1 - e-m = 1 – exp(-PPR/α)
      ##> 
      ##> The mono-exponential fit constant for a graph of p(≥1) vs PPR is 
      ##> α mg/mmol; at this PPR value, m = 1 and ~63% of the vesicles in the 
      ##> population possess ≥1 flippase.
      ## Fit a monoexponential curve to the data
      subset_for_fit <- data.frame(
        x=z$"Protein per Phospholipid (mg/mmol)",
        y=z$"Probability >= 1 Flippase in Vesicle")
      Rmod <- nlrob(y ~ 1-exp(-x/a),data=subset_for_fit, start = list(a=1), maxit=40)
      z$"Fit Constant (a)" <- Rmod$coefficients
      z$"PPR at P = 0.5" <- -Rmod$coefficients * log(1-0.5)
      output <- list(Raw=z)
      ## Generate data to plot the results of the fit
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
        "Probability >= 1 Flippase in Vesicle"=y_predicted_from_fit,
        "Experimental Series"=unique(z$"Experimental Series"),
        "Experiment"=unique(z$"Experiment"),
        check.names=FALSE,
        stringsAsFactors=FALSE)
      ## Return
      return(output)
    }
  )

# Output ------------------------------------------------------------------
  return(processed_list_from_x)

}