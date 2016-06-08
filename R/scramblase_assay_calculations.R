#' @importFrom assertive.numbers assert_all_are_greater_than_or_equal_to
#' @importFrom assertive.numbers assert_all_are_less_than
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom magrittr equals
#' @importFrom magrittr extract
#' @importFrom magrittr extract2
#' @importFrom magrittr is_less_than
#' @importFrom magrittr is_weakly_less_than
#' @importFrom magrittr multiply_by
#' @importFrom magrittr set_names
#' @importFrom magrittr subtract
#' @importFrom minpack.lm nlsLM
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats median
#' @importFrom stats predict
scramblase_assay_calculations <- function(
  x,
  scale_to,
  ppr_scale_factor = 0.65,
  generation_of_algorithm = 2,
  force_through_origin = TRUE,
  split_by_experiment = TRUE){
# Set parameters ----------------------------------------------------------
  nlsControl <- list(
    minFactor=1/20480,
    maxit=100)
  formulae <- list(
    pre_fit = list(
      "TRUE" = as.formula("y ~ b * ( 1 - exp( -x / a ))"),
      "FALSE" = as.formula("y ~ b - c * exp( -x / a )")),
    first_generation_algorithm = list(
      "TRUE" = as.formula("y ~ b * ( 1 - exp( -x / a ) )"),
      "FALSE" = as.formula("y ~ b - c * exp( -x / a )")),
    second_generation_algorithm = list(
      "TRUE" = as.formula("y ~ b * ( 1 - ( 1 / sqrt( 1 + 784 * a * x ) ) * exp( ( -3872 * a * x ) / ( 1 + 784 * a * x ) ) )"),
      "FALSE" = as.formula("y ~ b - c * ( ( 1 / sqrt( 1 + 784 * a * x ) ) * exp( ( -3872 * a * x ) / ( 1 + 784 * a * x ) ) )")))

  generation_of_algorithm %<>%
    switch(
      "1" = "first_generation_algorithm",
      "2" = "second_generation_algorithm")
  
# Parsing spectra ---------------------------------------------------------
  spectralData <- x %>%
    extract2("Path") %>%
    lapply(parse_fluorimeter_output)
  
# Read out data -----------------------------------------------------------
  # What spectral time windows to extract?
  minAcquisitionTime <- spectralData %>%
    sapply(
      function(y){
        y %>%
          extract2("Min.Acquisition.Time.in.sec") %>%
          return()
      }
    ) %>%
    assert_all_are_less_than(0)
  
  spectralData %>%
    sapply(
      function(y){
        y %>%
          extract2("Max.Acquisition.Time.in.sec") %>%
          return()
      }
    ) %>%
    assert_all_are_greater_than_or_equal_to(x[["Timepoint of Measurement (s)"]])
  maxAcquisitionTime <- x %>%
      extract2("Timepoint of Measurement (s)") %>%
      max(na.rm = TRUE)
  
  # Average over 10 values before dithionite addition for activity baseline
  x[["Baseline Fluorescence"]] <- spectralData %>%
    vapply(
      function(z){
        indexesForAveraging <- z %>%
          extract2("Data") %>%
          extract2("Time.in.sec") %>%
          is_less_than(0) %>%
          which() %>%
          tail(n = 10)
        z %>%
          extract2("Data") %>%
          extract2("Fluorescence.Intensity") %>%
          extract(indexesForAveraging) %>%
          median(na.rm = TRUE) %>%
          return()
      },
      1)
  # Average over last 10 values (in common time range) for activity
  x[["Minimum Fluorescence"]] <- spectralData %>%
    seq_along() %>%
    vapply(
      function(z){
        stopIndexForAveraging <- spectralData %>%
          extract2(z) %>%
          extract2("Data") %>%
          extract2("Time.in.sec") %>%
          is_weakly_less_than(x[z,"Timepoint of Measurement (s)"]) %>%
          which() %>%
          max(na.rm = TRUE)
        indexesForAveraging <- seq(
          from=stopIndexForAveraging-9,
          to=stopIndexForAveraging)
        spectralData %>%
          extract2(z) %>%
          extract2("Data") %>%
          extract2("Fluorescence.Intensity") %>%
          extract(indexesForAveraging) %>%
          median(na.rm = TRUE) %>%
          return()
      },
      1
    )

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
  processedListFromX <- inputListFromX %>%
    lapply(
      function(z){
        ## Ensure that there's a data point with liposomes ONLY as a unique 
        ## reference point
        indexOfLiposomesOnlyData <- z %>%
          extract2("Protein in Reconstitution (mg)") %>%
          equals(0) %>%
          which()
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
        z[["Relative Fluorescence Reduction"]]<- z[["Fluorescence Reduction"]] -
          z[["Fluorescence Reduction"]][indexOfLiposomesOnlyData]
        ## Calculate PPR
        z[["Protein per Phospholipid (mg/mmol)"]] <- z %>%
          calculate_ppr(ppr_scale_factor = ppr_scale_factor)
        ## Calculate p>=1Scramblase/Liposome
        y <- z %>%
          extract2("Relative Fluorescence Reduction")
        y0 <- y %>%
          extract(indexOfLiposomesOnlyData)
        if(scale_to == "model"){
          fit_prep <- z %>%
            fit_prep()
          fitStart <- if(force_through_origin){
            list(
              a = fit_prep[["estimated_a"]],
              b = fit_prep[["estimated_b"]])
          } else {        
            start = list(
              a = fit_prep[["estimated_a"]],
              b = fit_prep[["estimated_b"]],
              c = fit_prep[["estimated_b"]])
          }
          rMod <- minpack.lm::nlsLM(
            formula = formulae[["pre_fit"]][[as.character(force_through_origin)]],
            data = fit_prep[["fit_set"]],
            start = fitStart,
            control = nlsControl)
          yMax <- rMod %>%
            coef() %>%
            extract2("b")
        } else {
          yMax <- z %>%
            extract2("Relative Fluorescence Reduction") %>%
            max(na.rm = TRUE)
        }
        z[["Probability >= 1 Scramblase in Vesicle"]] <- (y-y0)/(yMax-y0)
        return(z)
      })
  # Is separate handling of experiments required?
  if(!split_by_experiment){
    processedListFromX <- processedListFromX %>%
      rbind.fill()
    processedListFromX <- processedListFromX %>%
      split(processedListFromX$"Experimental Series")
  }
  # Calculatons B
  processedListFromX %<>%
    lapply(
      function(z){
        fit_prep <- z %>%
          fit_prep()
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
          formula = formulae[[generation_of_algorithm]][[as.character(force_through_origin)]],
          data = fit_prep[["fit_set"]],
          start = if(force_through_origin){
            list(a = fit_prep[["estimated_a"]], b = 1)
          } else {        
            list(a = fit_prep[["estimated_a"]], b = 1, c = 1)
          },
          control = nlsControl)
        gc()
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

calculate_ppr <- function(x, ppr_scale_factor = 0.65){
  ppr <- x[["Protein in Reconstitution (mg)"]]/x[["Lipid in Reconstitution (mmol)"]]
  if(!is.null(ppr_scale_factor)){
    ppr <- ppr/ppr_scale_factor
  }
  return(ppr)
}


# Helper Functions --------------------------------------------------------
fit_prep <- function(z){
  fit_set <- z %>%
    extract2(
      c("Protein per Phospholipid (mg/mmol)",
        "Probability >= 1 Scramblase in Vesicle")) %>%
    set_names(c("x", "y"))
  ### Determine a sensible start point for 'a'
  pointSixY <- fit_set %>%
    extract2("y") %>%
    max(na.rm = TRUE) %>%
    multiply_by(0.6)
  estimated_a <- fit_set %>%
    extract2("x") %>%
    extract(
      fit_set %>%
        extract2("y") %>%
        subtract(pointSixY) %>%
        abs() %>%
        which.min())
  estimated_b <- z %>%
    extract2("Relative Fluorescence Reduction") %>%
    max(na.rm = TRUE)
  
  list(
    fit_set = fit_set,
    estimated_a = estimated_a,
    estimated_b = estimated_b) %>%
    return()
}

