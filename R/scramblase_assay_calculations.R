#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
scramblase_assay_calculations <- function(
  x,
  scale_to,
  ppr_scale_factor = 0.65,
  generation_of_algorithm = 2,
  force_through_origin = TRUE,
  split_by_experiment = TRUE,
  r_bar = 88,
  sigma_r_bar = 28){
# Set parameters ----------------------------------------------------------
  nlsControl <- list(
    minFactor=1/20480,
    maxit=100)
  formulae <- list(
    pre_fit = list(
      "TRUE" = y ~ b * ( 1 - exp( -x / a ) ),
      "FALSE" = y ~ b - c * exp( -x / a )),
    first_generation_algorithm = list(
      "TRUE" = y ~ b * ( 1 - exp( -x / a ) ),
      "FALSE" = y ~ b - c * exp( -x / a )),
    second_generation_algorithm = list(
      "TRUE" = y ~ b * ( 1 - ( 1 / sqrt( 1 + sigma_r_bar^2 * a * x ) ) * exp( ( -a * r_bar^2 / 2 * x ) / ( 1 + sigma_r_bar^2 * a * x ) ) ),
      "FALSE" = y ~ b - c * ( ( 1 / sqrt( 1 + sigma_r_bar^2 * a * x ) ) * exp( ( -a * r_bar^2 / 2 * x ) / ( 1 + sigma_r_bar^2 * a * x ) ) )))

  generation_of_algorithm %<>%
    switch(
      "1" = "first_generation_algorithm",
      "2" = "second_generation_algorithm")
  
# Parsing spectra ---------------------------------------------------------
  spectralData <- x %>%
    magrittr::extract2("Path") %>%
    lapply(parse_fluorimeter_output)
  
# Read out data -----------------------------------------------------------
  # What spectral time windows to extract?
  ## Just checking ...
 
  
 
  maxAcquisitionTime <- x %>%
      magrittr::extract2("Timepoint of Measurement (s)") %>%
      max(na.rm = TRUE)
  
  # Average over 10 values before dithionite addition for activity baseline
  x[["Baseline Fluorescence"]] <- spectralData %>%
    vapply(
      function(spectral_data_i){
        spectral_data_i %>% 
          magrittr::extract(.$Time.in.sec < 0, ) %>% 
          utils::tail(10) %>% 
          magrittr::extract2("Fluorescence.Intensity") %>% 
          stats::median(na.rm = TRUE)
      },
      1)
  # Average over last 10 values (in common time range) for activity
  x[["Minimum Fluorescence"]] <- mapply(
    function(spectral_data_i, timepoint_of_measurement_s)
    {
      # TODO: this behaviour matches the previous implementation, but using <=
      # seems slightly more intuitive
      spectral_data_i %>% 
        magrittr::extract(.$Time.in.sec < timepoint_of_measurement_s, ) %>% 
        utils::tail(10) %>% 
        magrittr::extract2("Fluorescence.Intensity") %>% 
        stats::median(na.rm = TRUE)
    },
    spectralData,
    x$"Timepoint of Measurement (s)"
  )

  # Apply volume correction factors as needed
  volumeCorrectionFactor <- x[["Fluorescence Assay Vol. with DT (ul)"]] / 
    x[["Fluorescence Assay Vol. w/o DT (ul)"]]
  x[["Minimum Fluorescence, Volume Corrected"]] <- x[["Minimum Fluorescence"]] *
    volumeCorrectionFactor
  ## Calculate activity reduction
  x[["Fluorescence Reduction"]] <- 1 - 
    x[["Minimum Fluorescence, Volume Corrected"]] / 
    x[["Baseline Fluorescence"]]

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
          magrittr::extract2("Protein Reconstituted (mg)") %>%
          magrittr::equals(0) %>%
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
        ##> P(>=1 flippase) = (y - yo)/(yMax - yo),
        ##> where yo is the percent reduction obtained with liposomes, yMax is the 
        ##> maximum percentage reduction observed and p is the probability that a 
        ##> particular vesicle in the ensemble is 'flippase-active', i.e it 
        ##> possesses at least one flippase.
        ## Calculate the relative fluorescence reduction
        z[["Relative Fluorescence Reduction"]]<- z[["Fluorescence Reduction"]] -
          z[["Fluorescence Reduction"]][indexOfLiposomesOnlyData]
        ## Calculate PPR
        z[["Protein per Phospholipid (mg/mmol)"]] <- z %>%
          calculate_ppr(ppr_scale_factor = ppr_scale_factor)
        ## Calculate p>=1Scramblase/Liposome
        y <- z %>%
          magrittr::extract2("Relative Fluorescence Reduction")
        y0 <- y %>%
          magrittr::extract(indexOfLiposomesOnlyData)
        if(scale_to == "model"){
          fit_prep <- z %>%
            fit_prep(y = "Relative Fluorescence Reduction")
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
            stats::coef() %>%
            magrittr::extract2("b")
        } else {
          yMax <- z %>%
            magrittr::extract2("Relative Fluorescence Reduction") %>%
            max(na.rm = TRUE)
        }
        z[["Probability >= 1 Scramblase in Vesicle"]] <- (y-y0)/(yMax-y0)
        return(z)
      })
  # Is separate handling of experiments required?
  if(!split_by_experiment){
    processedListFromX <- processedListFromX %>%
      plyr::rbind.fill()
    processedListFromX <- processedListFromX %>%
      split(processedListFromX$"Experimental Series")
  }
  # Calculatons B
  processedListFromX %<>%
    lapply(
      function(z){
        fit_prep <- z %>%
          fit_prep(y = "Probability >= 1 Scramblase in Vesicle")
        ##> Reference: Menon 2011, "Opsin is a phospholipid flippase", 
        ##> supplementary material 01, p2-3, section 
        ##> "Protein-dependence of lipid flipping; analysis of Figure 1f."
        ##> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3057128/#SD1
        ##> 
        ##> The dependence of P(>=1 flippase) on PPR was analyzed as follows.
        ##> Definitions:
        ##>   f, number of flippases used for reconstitution
        ##>   v, number of vesicles(
        ##>   m, number of flippases per vesicle (=f/v)
        ##>   PPR, mg protein per mmol phospholipid
        ##> To calculate P(>=1 flippase) as a function of the PPR, we assume that 
        ##> reconstitution of opsin/rhodopsin molecules into vesicles occurs 
        ##> independently and that the vesicles are identical and may have more
        ##> than one flippase. The probability that a flippase will be 
        ##> reconstituted into a particular vesicle is therefore 1/v. On 
        ##> reconstituting f flippases into an ensemble of v vesicles, the 
        ##> probability P(k) that a particular vesicle contains k flippases is 
        ##> given by the binomial formula:
        ##>   
        ##>   P(k) = C(f,k)*(1/v)^k*(1-1/v)^(f-k)
        ##> 
        ##> where C(f,k) is the number of combinations of k elements chosen 
        ##> from a set of size f.
        ##> 
        ##> Because f and v are both large, it is convenient to use the Poisson 
        ##> approximation:
        ##>   
        ##>   P(k) = (m^k/k!)e^(-m), where m = f/v is the average number of 
        ##> flippases per vesicle
        ##>
        ##> The probability of a particular vesicle having no flippases is 
        ##> P(0) = e^(-m); therefore, the probability that a vesicle has one or 
        ##> more flippases, i.e, is active in the flippase assay, is
        ##> 
        ##> P(>=1) = 1-P(0) = 1 - e^(-m)
        ##> 
        ##> The average number of flippases per vesicle, m, is proportional to PPR
        ##> and can be written as m = PPR/alpha, where alpha is a constant with units of
        ##> mg/mmol. 
        ##> Thus,
        ##> 
        ##> P(>=1) = 1 - e^(-m) = 1 - exp(-PPR/alpha)
        ##> 
        ##> The mono-exponential fit constant for a graph of P(>=1) vs PPR is 
        ##> alpha mg/mmol; at this PPR value, m = 1 and ~63% (= 1 - e^(-1)) of the 
        ##> vesicles in the population possess >=1 flippase.
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
        z[["Fit Constant (a)"]] <- stats::coef(rMod)[["a"]]
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
        yPredictedFromFit <- stats::predict(
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
  ppr <- x[["Protein Reconstituted (mg)"]]/x[["Lipid in Reconstitution (mmol)"]]
  if(!is.null(ppr_scale_factor)){
    ppr <- ppr/ppr_scale_factor
  }
  return(ppr)
}


# Helper Functions --------------------------------------------------------
fit_prep <- function(z, y = c("Relative Fluorescence Reduction", "Probability >= 1 Scramblase in Vesicle")){
  y %<>%
    match.arg(
      choices = c(
        "Relative Fluorescence Reduction",
        "Probability >= 1 Scramblase in Vesicle"),
      several.ok =  FALSE)
  fit_set <- z %>%
    magrittr::extract(
      c("Protein per Phospholipid (mg/mmol)",
        y)) %>%
    magrittr::set_names(c("x", "y"))
  ### Determine a sensible start point for 'a'
  pointSixY <- fit_set %>%
    magrittr::extract2("y") %>%
    max(na.rm = TRUE) %>%
    magrittr::multiply_by(0.6)
  estimated_a <- fit_set %>%
    magrittr::extract2("x") %>%
    magrittr::extract(
      fit_set %>%
        magrittr::extract2("y") %>%
        magrittr::subtract(pointSixY) %>%
        abs() %>%
        which.min())
  estimated_b <- z %>%
    magrittr::extract2(y) %>%
    max(na.rm = TRUE)
  
  list(
    fit_set = fit_set,
    estimated_a = estimated_a,
    estimated_b = estimated_b) %>%
    return()
}

