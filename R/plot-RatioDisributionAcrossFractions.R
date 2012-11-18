#' @title RatioDistributionAcrossFractions
#' @aliases RatioDistributionAcrossFractions
#' @description Visuallay exploring the need of ratio scaling to make 
#' fractions of the protein profiling experiment comparable.
#' @details It is possible that different fractions/experiments deliver ratio
#' distributions so widely differing in range/dimension from each other
#' that scaling (e.g. z-scoring) may be advisable.
#' 
#' This plot explores that need visually by using violin and box plots supplied 
#' by \pkg{\link{ggplot2}}.
#' 
#' For the dataset at hand scaling is not judged to be necessary.
#' @return Returns a \code{\link{list}} containing a  
#' \code{\link[ggplot2]{ggplot}} object representing the plot and 
#' \code{\link{character}} sting representing the legend with the handles 
#' \code{Plot} and \code{Legend}, respectively.
#' @author Johannes Graumann
#' @export
#' @seealso \code{\link[ggplot2]{geom_violin}} \code{\link[ggplot2]{geom_boxplot}}
#' @keywords distribution hplot datasets
#' @examples
#' require(ggplot2)
#' # Plot
#' Figure <- RatioDistributionAcrossFractions()
#' Figure$Plot
RatioDistributionAcrossFractions <- function(){
  # Discard reverse/contaminants
  tmpData <- pGroups
  tmpData <- tmpData[tmpData$Reverse != "+",]
  tmpData <- tmpData[tmpData$Contaminant != "+",]
  # Extract columns containing NOT-NORMALIZED ratios and assemble plottable frame
  tmpData <- pGroups[grep(pattern="^Ratio H/L Exp. [[:alpha:]]{1}",names(pGroups))]
  tmpList <- lapply(
    sub(pattern="^Ratio H/L Exp. ",replacement="",names(tmpData)),
    function(x){
      tmpDF <- data.frame(
        Fraction = rep(x,nrow(tmpData)),
        Ratio = tmpData[[paste("Ratio H/L Exp.",x)]])
    })
  tmpData <- do.call("rbind",tmpList)
  # Calculate number of ratios for scaling the box width
  ratiocount <- sapply(split(tmpData,tmpData$Fraction),function(x){sum(!(is.nan(x$Ratio)))})
  relratiocount <- ratiocount/max(ratiocount)
  # Construct Plot
  tmpPlot <- ggplot(data=tmpData,aes(x=Fraction,y=Ratio))
  tmpPlot <- 
    tmpPlot + 
    scale_y_log10() +
    geom_violin() +
    geom_boxplot(width=0.25) + #aes(size=relratiocount*0.25)) +
    labs(y=expression("Ratio "*frac(Input,Fraction)))
  return(
    list(
      Plot=tmpPlot,
      Legend=NA_character_))
}
