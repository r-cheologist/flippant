#' @export
prelimFigure1 <- function(x,directory=tempdir()){
  require(reshape)
  # Prep output
  output <- list()
  # Extract columns 
  tmpData <- x[c("Majority protein IDs","Gene names","Peptides Exp. D","Ratio H/L Exp. D",names(x)[grep(pattern="^Rel. Ratio",names(x))],"Rel.CC","transmembrane_domain")]
  # Fill in missing gene names with Uniprot IDs
  subsetter <- tmpData[["Gene names"]]==""
  tmpData[subsetter,"Gene names"] <- tmpData[subsetter,"Majority protein IDs"]
  tmpData <- tmpData[2:ncol(tmpData)]
  names(tmpData)[1] <- "ID"
  # Tag peptide count/raw ratio onto label
  tmpData$ID <- paste(
    tmpData$ID,
    round(tmpData$Rel.CC,2),
    tmpData[["Peptides Exp. D"]],
    round(RatioInversion(tmpData[["Ratio H/L Exp. D"]]),2),
    tmpData$transmembrane_domain,
    sep="|")
  # Strip down fraction names
  names(tmpData) <- sub(pattern="Rel. Ratio L/H Exp. ",replacement="",names(tmpData))
  # Cut into CC brackets
  tmpData$CorrelationQuality <- cut(x=tmpData$Rel.CC,breaks=c(-1,-0.95,-0.9,-0.85,-0.8,0,0.8,0.85,0.9,0.95,1.0))
  # Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
  tmpData <- melt(data=tmpData,measure.vars=LETTERS[1:13],na.rm=TRUE)
  tmpData$variable <- match(tmpData$variable,LETTERS)
  # Prep the flippase Activity analogously
  tmpFA <- flippaseActivity
  tmpFA[["Av. Spec. Activity"]] <- Relativate(tmpFA[["Av. Spec. Activity"]])
  names(tmpFA) <- c("variable","value")
  tmpFA$variable <- match(tmpFA$variable,LETTERS)
  # Plot A
  tmpPlot <- ggplot(data=tmpData,aes(x=variable,y=value))
  tmpPlot <- tmpPlot +
    geom_line(aes(group=ID)) +
    facet_wrap(~CorrelationQuality) +
    geom_line(data=tmpFA,col="red") +
    labs(
      title="Rel. Ratio Profiles Faceted by Correlation Coefficient Bracket",
      x="Fraction",
      y="Relative Ratio L/H, Relative Flippase Activity")
  output[["A"]] <- tmpPlot
  # Correlation-bracket plots
  tmpData <- split(tmpData,tmpData$CorrelationQuality)
  closeinspection <- list(B="(0.85,0.9]",C="(0.9,0.95]",D="(0.95,1]")
  for(y in names(closeinspection)){
    tmpPlot <- ggplot(data=tmpData[[closeinspection[[y]]]],aes(x=variable,y=value))
    tmpPlot <- tmpPlot +
      geom_line(aes(color=ID)) +
      geom_line(data=tmpFA,col="black") +
      labs(
        x="Fraction",
        y="Relative Ratio L/H, Relative Flippase Activity",
        title=paste("Correlation Coefficient Bracket:",closeinspection[[y]]),
        color="ID|CC|Peptides|Raw Ratio")
    if(y  =="C"){
      tmpPlot <- tmpPlot +
        guides(col=guide_legend(ncol=2))
    } else if(y == "B"){
      tmpPlot <- tmpPlot +
        guides(col=guide_legend(ncol=3))
    }
#     if(y=="C"){
#       tmpPlot <- tmpPlot +
#         guides(col=guide_legend(ncol=2))
#     }
    output[[y]] <- tmpPlot
  }
  png(
    filename=file.path(directory,"Figure1.png"),
    width=210,
    height=297,
    units="mm",
    res=300
  )
  library(grid)
  grid.newpage()
  pushViewport(viewport(x=0,y=1,width=1,height=0.5,,just=c("left","top"),name="A"))
  upViewport()
  pushViewport(viewport(x=0,y=0.5,width=0.5,height=0.5,,just=c("left","top"),name="B"))
  upViewport()
  pushViewport(viewport(x=0.5,y=0.5,width=0.5,height=0.5,,just=c("left","top"),name="C"))
  upViewport()
  print(output$A,vp="A")
  grid.text(label="A",x=0.025,y=0.975,vp="A")
  print(output$B,vp="B")
  grid.text(label="B",x=0.025,y=0.975,vp="B")
  print(output$C,vp="C")
  grid.text(label="C",x=0.025,y=0.975,vp="C")
  dev.off()
  return(invisible(output))
}