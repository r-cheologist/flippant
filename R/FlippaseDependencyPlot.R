#' @seealso \code{\link{ParseFluorometerData}}
FlippaseDependencyPlot <- function(x=NA,ReactionVolumes=c(2000,2040),ePC=4.5){
  x <- list(
    "15ul"=list(
      Path="../FlurometerParsing/Erg9_TE_assay/TE+/RK_30_May-2013_Erg9_TE_plus_15ul.td"),
    "40ul"=list(
      Path="../FlurometerParsing/Erg9_TE_assay/TE+/RK_30_May-2013_Erg9_TE_plus_40ul.td"),
    "80ul"=list(
      Path="../FlurometerParsing/Erg9_TE_assay/TE+/RK_30_May-2013_Erg9_TE_plus_80ul.td"),
    "160ul"=list(
      Path="../FlurometerParsing/Erg9_TE_assay/TE+/RK_30_May-2013_Erg9_TE_plus_160ul.td"))
  ReactionVolumes=c(2000,2040)
  ePC=4.5
  #######################
  # Check prerequisites #
  #######################
  if(is.na(x)[1]){
    stop("'x' must be defined.")
  }
  if(!is.list(x)){
    stop("'x' must be of class 'list'.")
  }
  if(is.null(names(x))){
    stop("All elements on level 1 of 'x' must be named.")
  }
  if(!all(vapply(x,is.list,TRUE))){
    stop("All elements in 'x' must be of class 'list'.")
  }
  if(listDepth(x) != 2){
    stop("'x' must be of depth 2.")
  }
  if(any(vapply(x,function(y){is.null(names(y))},TRUE))){
    stop("All elements on level 2 of 'x' must be named.")
  }
  nesLev2Names <- "Path"
  facLev2Names <- c("ReactionVolumes","ePC")
  legLev2Names <- c(nesLev2Names,facLev2Names)
  if(!all(vapply(x,function(y){identical(intersect(names(y),legLev2Names),names(y))},TRUE))){
    stop("All names on level 2 of 'x' must be from '",paste(legLev2Names,collapse="','"),"'.")
  }
  if(!any(vapply(x,function(y){identical(intersect(names(y),nesLev2Names),names(y))},TRUE))){
    stop("Elements with names '",paste(nesLev2Names,collapse="','"),"' are compulsory on level 2 of 'x'.")
  }
  if(!all(vapply(x,function(y){is.character(y$Path) & length(y$Path) ==1},TRUE))){
    stop("All 'Path' elements in 'x' must be single length 'character' objects.")
  }
  if(!all(vapply(
    x,
    function(y){
      output <- TRUE
      if("ReactionVolumes" %in% names(y)){
        if(!is.numeric(y$ReactionVolumes) | length(y$ReactionVolumes) != 2){
          output <- FALSE
        }
      }
      return(output)
    },
    TRUE))){
    stop("All 'ReactionVolumes' elements in 'x' must be 'numeric' objects of length 2.")
  }
  if(!is.numeric(ePC) | length(ePC) != 1){
    stop("'ePC' must be 'numeric' objects of length 1.")
  }
  if(!all(vapply(
    x,
    function(y){
      output <- TRUE
      if("ePC" %in% names(y)){
        if(!is.numeric(y$ePC) | length(y$ePC) != 1){
          output <- FALSE
        }
      }
      return(output)
    },
    TRUE))){
    stop("All 'ePC' elements in 'x' must be 'numeric' objects of length 1.")
  }
  if(!all(file.exists(vapply(x,function(y){y$Path},"A")))){
    stop("All 'Path' elements in 'x' must refer to existing files.")
  }
  if(any(file.access(vapply(x,function(y){y$Path},"A"),mode=4) == -1)){
    stop("All 'Path' elements in 'x' must refer to readable files.")
  }
  ##############
  # Processing #
  ##############
  # Parse the data and extract whats needed
  #########################################
  # Parsing
  tmpData <- lapply(
    x,
    function(y){ParseFluorometerData(SpecFile=y$Path)})
  # What intensity to extract?
  minAT <- vapply(tmpData,function(y){y$"Minimal Acquisition Time (s)"},1)
  if(!all(minAT[1] == minAT)){
    stop("Minimum acquisition times are not identical - aborting.")
  }
  minAT <- unique(minAT)
  tmpAT <- vapply(tmpData,function(y){y$"Maximal Acquisition Time (s)"},1)
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
  tmpRV <- vapply(
    tmpData,
    function(y){
      if(is.null(y$"Reaction Volumes")){
        return(ReactionVolumes)
      } else {
        return(y$"Reaction Volumes")
      }
    },
    c(1,2))
  tmpCF <- vapply(tmpRV,function(y){y[2]/y[1]},1)
  activityIntensity <- activityIntensity * tmpCF
  # Calculate and relativate intensity differences
  deltaIntensity <- 1-(activityIntensity/baselineIntensity)
}