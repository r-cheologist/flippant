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
  names(tmpData) <- sub(pattern="H/L",replacement="L/H",x=names(tmpData))
  # Where's the max?'
  ratioMax <- sapply(split(tmpData,seq(nrow(tmpData))),which.max)
  # Assemble and filter
  simpleData <- cbind(saveData,tmpData)
  simpleData <- simpleData[ratioMax == 4,]
  simpleData <- simpleData[order(simpleData["Ratio L/H Exp. D"],decreasing=TRUE),]
  return(simpleData)
}