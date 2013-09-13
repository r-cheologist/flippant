#' @seealso \code{\link{ParseFluorometerData}}, 
#' \code{\link{ParseFluorometerData2}}, \code{\link{TimepointOfMeasurement}}
#' @import robustbase
FlippaseDependencySeries <- function(x)
{
  x <- data.frame(
    Path = c(
      "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
      "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
      "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
      "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
      "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt"),
    "Extract Volume (ul)" = c(0,15,40,75,150),
#     "Reaction Volume Blank (ul)" = rep(2000,4),
    "Reaction Volume Reconstituted (ul)" = rep(2040,5),
    "Concentration Egg PC (mM)" = rep(4.5,5),
    "Extract Protein Concentration (mg/ml)" = rep(0.67,5),
#     "Timepoint of Measurement (s)",
    check.names=FALSE,
    stringsAsFactors=FALSE)
      
#   AcquisitionTimePoint <- MinimumAcquisitionTime(sapply(x,function(y){y$Path}))
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
      "Extract Protein Concentration (mg/ml)"),
    Class = c(
      "character",
      "numeric",
      "numeric"))
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
      "Reaction Volume Blank (ul)",
      "Reaction Volume Reconstituted (ul)",
      "Concentration Egg PC (mM)",
      "Timepoint of Measurement (s)"),
    Class = c(
      "numeric",
      "numeric",
      "numeric",
      "numeric"),
    Default = c(
      2000,
      2040,
      4.5,
      NA))
  missing <- which(!(facultatives$Name %in% names(x)))
  if(length(missing) != 0){
    for (y in missing){
      if(facultatives$Name[y] == "Timepoint of Measurement (s)"){
        warning(
          "Providing missing column '",
          facultatives$Name[y],
          "' from spectra ('Path'). Not appropriate for multiple extract ",
          "comparisons. You have been warned. See function ",
          "'TimepointOfMeasurement()'.")
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
  # Parse the data and extract whats needed
  #########################################
  # Parsing
  #########
  tmpData <- lapply(
    x$Path,
    function(y){ParseFluorometerData2(SpecFile=y)})
  # What spectral time windows to extract?
  ########################################
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
  # Calculations
  ##############
  # Average over first 10 values for activity baseline
  baselineIntensity <- vapply(
    tmpData,
    function(y){
      tmpFrom <- min(which(y$Data$"Time (s)" >= minAT))
      baselineSS <- seq(from=tmpFrom,to=tmpFrom+9)
      return(median(y$Data$"Fluorescense Intensity"[baselineSS],na.rm=TRUE))
    },
    1)
  # Average over last 10 values (in common time range) for activity
  activityIntensity <- vapply(
    tmpData,
    function(y){
      tmpTo <- max(which(y$Data$"Time (s)" <= maxAT))
      activitySS <- seq(from=tmpTo-9,to=tmpTo)
      return(median(y$Data$"Fluorescense Intensity"[activitySS],na.rm=TRUE))
    },
    1)
  # Apply volume correction factors as needed
  correctionFactor <- x$"Reaction Volume Reconstituted (ul)"/x$"Reaction Volume Blank (ul)"
  activityIntensity <- activityIntensity * correctionFactor
  # Calculate and relativate intensity differences
  deltaIntensity <- 1-(activityIntensity/baselineIntensity)
  # Calculate Protein to egg phosphatidylcholin ratio
  proteinAmount <- x$"Extract Protein Concentration (mg/ml)" * x$"Extract Volume (ul)"
  proteinToEpcRatio <- proteinAmount/x$"Concentration Egg PC (mM)"
  # DEBUG: plot(proteinToEpcRatio,deltaIntensity)
  # Fit a monoexponential curve to the data
  #########################################
  tmpData <- data.frame(
    x=proteinToEpcRatio,
    y=deltaIntensity)
  library(robustbase)
#  curve(-1+-1/x,0,10)
#   http://people.richland.edu/james/lecture/m116/logs/models.html
#   y = C ( 1 - e-kt ), k > 0
  Rmod <- nlrob(y ~ C*(1-exp(-k*x)),data=tmpData, start = list(C=max(tmpData$x), k=1))
  mod <- nls(y ~ C*(1-exp(-k*x)),data=tmpData, start = list(C=max(tmpData$x), k=1))
  newData <- data.frame(x=seq(from=min(proteinToEpcRatio),to=max(proteinToEpcRatio),length.out=1000))
  plot(tmpData)
  lines(newData$x,predict(mod,newdata=newData),col="red")
  lines(newData$x,predict(Rmod,newdata=newData),col="green")
  
#   mod <- nlrob(y ~ a*exp(x*b), data=tmpData, start = list(a = -1, b = -1))
  ## Test Example:
#   char_fx <- "-1*exp(x*-1)"
#   curve(-1*exp(x*-1),0,10)
#   fx <- expression(-1*exp(x*-1))
#   plot(fx,0,10)
#   tmpEnv <- new.env()
#   assign("x",seq(10),envir=tmpEnv)
#   curve(-1*exp(x*-1),0,10)
#   
#   dfx <- deriv(fx,"x")
#   
#   tmpEnv <- new.env()
#   assign("x",0,envir=tmpEnv)
#   curve(eval(dfx,envir=tmpEnv),0,10,add=TRUE,col="red")
#   eval(expr=dfx,envir=tmpEnv)
  ###################
  # Assemble output #
  ###################
}