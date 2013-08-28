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
  if(!all(sapply(x,is.list))){
    stop("All elements in 'x' must be of class 'list'.")
  }
  if(listDepth(x) != 2){
    stop("'x' must be of depth 2.")
  }
  if(any(sapply(x,function(x){is.null(names(x))}))){
    stop("All elements on level 2 of 'x' must be named.")
  }
  nesLev2Names <- "Path"
  facLev2Names <- c("ReactionVolumes","ePC")
  legLev2Names <- c(nesLev2Names,facLev2Names)
  if(!all(sapply(x,function(x){identical(intersect(names(x),legLev2Names),names(x))}))){
    stop("All names on level 2 of 'x' must be from '",paste(legLev2Names,collapse="','"),"'.")
  }
  if(!any(sapply(x,function(x){identical(intersect(names(x),nesLev2Names),names(x))}))){
    stop("Elements with names '",paste(nesLev2Names,collapse="','"),"' are compulsory on level 2 of 'x'.")
  }
  if(!all(sapply(x,function(x){is.character(x$Path) & length(x$Path) ==1}))){
    stop("All 'Path' elements in 'x' must be single length 'character' objects.")
  }
  if(!all(sapply(
    x,
    function(x){
      output <- TRUE
      if("ReactionVolumes" %in% names(x)){
        if(!is.numeric(x$ReactionVolumes) | length(x$ReactionVolumes) != 2){
          output <- FALSE
        }
      }
      return(output)
    }))){
    stop("All 'ReactionVolumes' elements in 'x' must be 'numeric' objects of length 2.")
  }
  if(!all(sapply(
    x,
    function(x){
      output <- TRUE
      if("ePC" %in% names(x)){
        if(!is.numeric(x$ePC) | length(x$ePC) != 1){
          output <- FALSE
        }
      }
      return(output)
    }))){
    stop("All 'ePC' elements in 'x' must be 'numeric' objects of length 1.")
  }
  if(!all(file.exists(sapply(x,function(x){x$Path})))){
    stop("All 'Path' elements in 'x' must refer to existing files.")
  }
  if(any(file.access(sapply(x,function(x){x$Path}),mode=4) == -1)){
    stop("All 'Path' elements in 'x' must refer to readable files.")
  }
}