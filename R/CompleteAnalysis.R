#' @export
CompleteAnalysis <- function(dir=tempdir()){
  # Prep the data for correlation analysis
  fig1Data <- DataPrep(pGroups)
  # Correlate (Spearman)
  fig1Data$CC <- Correlate(
    x=fig1Data,
    y=Relativate(flippaseActivity[["Av. Spec. Activity"]]))
  # Merge with original data
  subsetter <- sapply(rownames(fig1Data),function(x){which(pGroups$id == x)})
  fig1Data <- cbind(pGroups[subsetter,],fig1Data)
  # Save
  write.table(
    x=fig1Data,
    file=file.path(dir,paste(Sys.Date(),"Correlation.txt",sep="_")),
    sep="\t",
    row.names=FALSE,
    col.names=TRUE)
  # Produce first figure
  test <- Figure1(fig1Data,dir=dir)
}

#   # Where do the ratio minima reside 1?
#   #ratioMin <- sapply(split(tmpData,seq(nrow(tmpData))),which.min)
#   # Impute missing values
#   #library(pcaMethods)
#   #imputationMethod <- "bpca"
#   #scaleMethod <- "none"
#   #centering <- FALSE
#   #pcares <- pca(
#   #  object=t(tmpData),
#   #  method=imputationMethod,
#   #  scale=scaleMethod,
#   #  center=centering)
#   #tmpData <- as.data.frame(t(slot(pcares,"completeObs")))
#   #rm(pcares)
#   # Where do the ratio minima reside 2?
#   #imputedRatioMin <- sapply(split(tmpData,seq(nrow(tmpData))),which.min)
#   # 
#   # Correlate
#   correlationMethod <- "spearman"
#   cCoeffs <- sapply(
#     seq(nrow(tmpData)),
#     function(x){
#       pcp <- unlist(tmpData[x,])
#       cCoeff <- cor(
#         pcp,
#         flippaseActivity[["Av. Spec. Activity"]] * -1,
#         method=correlationMethod)
#       return(cCoeff)
#     })
#   # Assembly
#   names(tmpData) <- paste("Imputed",names(tmpData))
#   saveData$RatioMin <- ratioMin
#   saveData <- cbind(saveData,tmpData,stringsAsFactors=FALSE)
#   saveData$cCoeffs <- cCoeffs
#   saveData$ImputedRatioMin <- imputedRatioMin
#   # Order by cCoeff
#   saveData <- saveData[order(saveData$cCoeffs,decreasing=TRUE),]
#   # First x candidates
#   topX <- 50
#   write.table(
#     x = saveData[seq(topX),],
#     file = file.path(
#       "/tmp",
#       paste(
#         Sys.Date(),
#         "_Top",topX,
#         "_Cor-",correlationMethod,
#         "_Imput-",imputationMethod,
#         "_Scale-",scaleMethod,
#         "_Center-",centering,
#         ".txt",
#         sep="")
#     ),
#     sep = "\t",row.names=FALSE,col.names=TRUE)