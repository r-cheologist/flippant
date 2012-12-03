#' @export
AdditiveFilteredProfilePlots <- function(x,filter){
  require(reshape)
  # Prep output
  output <- list()
  # Filter candidates
  x <- x[na.omit(match(x=tolower(filter),table=tolower(x[["Gene names"]]))),]
  # Fill in missing gene names with Uniprot IDs
  subsetter <- x[["Gene names"]]==""
  x[subsetter,"Gene names"] <- x[subsetter,"Majority protein IDs"]
  names(x)[which(names(x)=="Gene names")] <- "ID"
  # Tag peptide count/raw ratio onto label
  x$ID <- paste(
    x$ID,
    round(x$CC,2),
    x[["Peptides Exp. D"]],
    round(x[["Ratio L/H Exp. D"]],2),
    x$transmembrane_domain,
    sep="|")
  # Extract columns 
  tmpRaw <- x[c("ID",names(x)[grep(pattern="^Ratio L/H",names(x))])]
  # Make column sums and Relativate - assemble into Data Frame
  tmpString <- Relativate(colSums(tmpRaw[grep(pattern="^Ratio L/H",names(tmpRaw))]))
  tmpSum <- data.frame(
    ID=paste(tmpRaw$ID,collapse=" + "),
    A=tmpString[1],
    B=tmpString[2],
    C=tmpString[3],
    D=tmpString[4],
    E=tmpString[5],
    F=tmpString[6],
    G=tmpString[7],
    H=tmpString[8],
    I=tmpString[9],
    J=tmpString[10],
    K=tmpString[11],
    L=tmpString[12],
    M=tmpString[13])
  # Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
  tmpSum <- melt(data=tmpSum,measure.vars=LETTERS[1:13],na.rm=TRUE)
  tmpSum$variable <- match(tmpSum$variable,LETTERS)
  # Prep the flippase Activity analogously
  tmpFA <- flippaseActivity
  tmpFA[["Av. Spec. Activity"]] <- Relativate(tmpFA[["Av. Spec. Activity"]])
  names(tmpFA) <- c("variable","value")
  tmpFA$variable <- match(tmpFA$variable,LETTERS)
  # Plot A
  sumPlot <- ggplot(data=tmpSum,aes(x=variable,y=value))
  sumPlot <- sumPlot +
    geom_line() +
    geom_line(data=tmpFA,col="red") +
    labs(
      title=tmpSum$ID[1],
      x="Fraction",
      y="Relative Summed Ratio L/H/Relative Flippase Activity")
  output[["A"]] <- sumPlot
  return(invisible(output))
}