#' @export
SimpleSort <- function(){
  # Discard reverse/contaminants
  saveData <- pGroups
  saveData <- saveData[saveData$Reverse != "+",]
  saveData <- saveData[saveData$Contaminant != "+",]
  rownames(saveData) <- NULL
  # Extract columns containing NOT-NORMALIZED ratios and assemble
  tmpData <- saveData[grep(pattern="^Ratio H/L Exp. [[:alpha:]]{1}",names(saveData))]
  # Invert the ratios
  tmpData <- RatioInversion(tmpData)
  # Where's the max?'
  ratioMax <- sapply(split(tmpData,seq(nrow(tmpData))),which.max)
  simpleData <- saveData[ratioMax == 4,]
  simpleData <- simpleData[order(simpleData["Ratio H/L Exp. D"]),]
  topX <- 50
  write.table(
    x = simpleData,
    file = file.path(
      "/tmp",
      paste(
        Sys.Date(),
        "MinRatioInFracD.txt",
        sep="")
    ),
    sep = "\t",row.names=FALSE,col.names=TRUE)
}