#' @importFrom minpack.lm nlsLM
scramblase_assay_calculations <- function(
  x,
  scale_to,
  generation_of_algorithm = 2,
  force_through_origin = FALSE,
  split_by_experiment = TRUE){
# Set parameters ----------------------------------------------------------
  nlsControl <- list(
    minFactor=1/20480,
    maxit=100)

# Parsing spectra ---------------------------------------------------------
  spectralData <- lapply(
    x[["Path"]],
    parse_fluorimeter_output)

# Read out data -----------------------------------------------------------
  # What spectral time windows to extract?
  minAcquisitionTime <- vapply(
    spectralData,
    function(y){
      y[["Min.Acquisition.Time.in.sec"]]},
    1)
  if(!all(minAcquisitionTime < 0)){
    stop("Minimum acquisition times are not all negative - aborting.")
  }
  maxAcquisitionTime <- vapply(
    spectralData,
    function(y){
      y[["Max.Acquisition.Time.in.sec"]]
      },
    1)
  if(any(maxAcquisitionTime < x[["Timepoint of Measurement (s)"]])){
    stop("At least one 'Timepoint of Measurement' is larger than the 
      corresponding spectrum acquisition. Aborting.")
  } else {
    maxAcquisitionTime <- max(x[["Timepoint of Measurement (s)"]])
  }
  # Average over 10 values before dithionite addition for activity baseline
  x[["Baseline Fluorescence"]] <- vapply(
    spectralData,
    function(z){
      indexesForAveraging <- tail(
        which(
          z[["Data"]][["Time.in.sec"]] < 0),
        n=10)
      return(
        median(
          z[["Data"]][["Fluorescence.Intensity"]][indexesForAveraging],
          na.rm=TRUE))
    },
    1)
  # Average over last 10 values (in common time range) for activity
  x[["Minimum Fluorescence"]] <- vapply(
    seq_along(spectralData),
    function(z){
      stopIndexForAveraging <- max(
        which(
          spectralData[[z]][["Data"]][["Time.in.sec"]] <= x[z,"Timepoint of Measurement (s)"]))
      indexesForAveraging <- seq(
        from=stopIndexForAveraging-9,
        to=stopIndexForAveraging)
      return(
        median(
          spectralData[[z]][["Data"]][["Fluorescence.Intensity"]][indexesForAveraging],
          na.rm=TRUE))
    },
    1)
  # Apply volume correction factors as needed
  volumeCorrectionFactor <- x[["Fluorescence Assay Vol. with DT (ul)"]]/x[["Fluorescence Assay Vol. w/o DT (ul)"]]
  x[["Minimum Fluorescence, Volume Corrected"]] <- x[["Minimum Fluorescence"]] * volumeCorrectionFactor
  ## Calculate activity reduction
  x[["Fluorescence Reduction"]] <- 1 - x[["Minimum Fluorescence, Volume Corrected"]]/x[["Baseline Fluorescence"]]

# Generate PPR v.s P(>=1 Scramblase/Liposome) -------------------------------
  # Split by Experiment
  x[["CombinedId"]] <- paste(
    x[["Experimental Series"]],
    x[["Experiment"]],
    sep="_")
  inputListFromX <- split(
    x,
    x[["CombinedId"]])

  # Calculations A
  processedListFromX <- lapply(
    inputListFromX,
    function(z){
      ## Ensure that there's a data point with liposomes ONLY as a unique 
      ## reference point
      indexOfLiposomesOnlyData <- which(z[["Protein in Reconstitution (mg)"]] == 0)
      if(length(indexOfLiposomesOnlyData) == 0){
        stop("Experimental series '",unique(z[["Experimental Series"]]),"' does 
          not have the required liposomes-ONLY ('Extract Volume (ul)' of '0') 
          data point. Aborting.")
      }
      if(length(indexOfLiposomesOnlyData) > 1){
        stop("Experimental series '",unique(z[["Experimental Series"]]),"' has 
          more than one liposomes-ONLY ('Extract Volume (ul)' of '0') data 
          point.")
      }  
      ##> End-point fluorescence reduction data from flippase activity assays 
      ##> were obtained for proteoliposomes generated over a range of PPR 
      ##> values.
      ##> The data were transformed according to the formula
      ##> p(≥1 flippase) = (y – yo)/(yMax – yo),
      ##> where yo is the percent reduction obtained with liposomes, yMax is the 
      ##> maximum percentage reduction observed and p is the probability that a 
      ##> particular vesicle in the ensemble is ‘flippase-active’, i.e it 
      ##> possesses at least one flippase.
      ## Calcualte the relative fluorescence reduction
      z[["Relative Fluorescence Reduction"]] <- z[["Fluorescence Reduction"]] - z[["Fluorescence Reduction"]][indexOfLiposomesOnlyData]
      ## Calculate PPR
      z[["Protein per Phospholipid (mg/mmol)"]] <- calculate_ppr(z)
      ## Calculate p>=1Scramblase/Liposome
      y <- z[["Relative Fluorescence Reduction"]]
      y0 <- z[["Relative Fluorescence Reduction"]][indexOfLiposomesOnlyData]
      if(scale_to == "model"){
        subsetForFit <- data.frame(
          x=z[["Protein per Phospholipid (mg/mmol)"]],
          y=z[["Relative Fluorescence Reduction"]])
        ### Determine a sensible start point for 'a'
        pointSixY <- max(
          subsetForFit[["y"]],
          na.rm=TRUE) * 0.6
        estimatedA <- subsetForFit[["x"]][which.min(abs(subsetForFit[["y"]] - pointSixY))]
        estimatedB <- max(
          z[["Relative Fluorescence Reduction"]],
          na.rm=TRUE)
        rMod <- minpack.lm::nlsLM(
          formula = if(force_through_origin){
            y ~ b * (1 - exp(-x/a))
          } else {
            y ~ b-c*exp(-x/a)
          },
          data = subsetForFit,
          start = if(force_through_origin){
            list(a=estimatedA,b=estimatedB)
          } else {        
            start = list(a=estimatedA,b=estimatedB,c=estimatedB)
          },
          control = nlsControl)
        yMax <- coef(rMod)[["b"]]
      } else {
        yMax <- max(
          z[["Relative Fluorescence Reduction"]],
          na.rm=TRUE)
      }
      z[["Probability >= 1 Scramblase in Vesicle"]] <- (y-y0)/(yMax-y0)
      return(z)
    })
  # Is separate handling of experiments required?
  if(!split_by_experiment){
    processedListFromX <- rbind.fill(processedListFromX)
    processedListFromX <- split(processedListFromX,processedListFromX$"Experimental Series")
  }
  # Calculatons B
  processedListFromX <- lapply(
    processedListFromX,
    function(z){
      subsetForFit <- data.frame(
        x=z[["Protein per Phospholipid (mg/mmol)"]],
        y=z[["Probability >= 1 Scramblase in Vesicle"]])
      ## Determine a sensible start point for 'a' in the monoexponential fit
      pointSixY <- max(
        subsetForFit[["y"]],
        na.rm=TRUE) * 0.6
      estimatedA <- subsetForFit[["x"]][which.min(abs(subsetForFit[["y"]] - pointSixY))]
      if(generation_of_algorithm == 1){
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
        rMod <- minpack.lm::nlsLM(
          formula = if(force_through_origin){
            y ~ b * (1 - exp(-x/a))
          } else {
            y ~ b-c*exp(-x/a)
          },
          data = subsetForFit,
          start = if(force_through_origin){
            list(a = estimatedA, b = 1)
          } else {        
            list(a = estimatedA, b = 1, c = 1)
          },
          control = nlsControl)
        gc()
      } else if(generation_of_algorithm == 2){
        ## p(≥1) = 1-(1/sqrt(1+784*α*PPR))*exp(-3872*α*PPR/(1+784*α*PPR))
        ## p(≥1) = 1 - e-m = 1 – exp(-PPR/α)
        # ## Fit a monoexponential curve to the data
        rMod <- minpack.lm::nlsLM(
          formula = if(force_through_origin){
            y ~ b * (1 - (1/sqrt(1 + 784 * a * x)) * exp(-3873 * a * x / (1 + 784 * a * x)))
          } else {
            y ~ b - c * (1/sqrt(1 + 784 * a * x)) * exp(-3873 * a * x / (1 + 784 * a * x))
          },
          data = subsetForFit,
          start = if(force_through_origin){
            list(a=estimatedA,b=1)
          } else {        
            list(a=estimatedA,b=1,c=1)
          },
          control = nlsControl)
        gc()
      }
      z[["Fit Constant (a)"]] <- coef(rMod)[["a"]]
      output <- list(
        Raw = z,
        FitObject = rMod)
      ## Generate data to plot the results of the fit
      xPredictedFromFit <- seq(
        from = min(
          z[["Protein per Phospholipid (mg/mmol)"]],
          na.rm=TRUE),
        to = max(
          z[["Protein per Phospholipid (mg/mmol)"]],
          na.rm=TRUE),
        length.out=200)
      yPredictedFromFit <- predict(
        rMod,
        newdata = list(x = xPredictedFromFit))
      output[["Fit"]] <- data.frame(
        "Protein per Phospholipid (mg/mmol)"=xPredictedFromFit,
        "Probability >= 1 Scramblase in Vesicle"=yPredictedFromFit,
        "Experimental Series"=unique(z[["Experimental Series"]]),
        "Experiment"=unique(z[["Experiment"]]),
        check.names=FALSE,
        stringsAsFactors=FALSE)
      ## Return
      return(output)
    })

# Output ------------------------------------------------------------------
  return(processedListFromX)

}

calculate_ppr <- function(x){
  return(x[["Protein in Reconstitution (mg)"]]/x[["Lipid in Reconstitution (mmol)"]])
}