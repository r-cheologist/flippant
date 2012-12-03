#' @export
FilteredProfilePlots <- function(x,filter){
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
  tmpRel <- x[c("ID",names(x)[grep(pattern="^Rel. Ratio",names(x))])]
  tmpRaw <- x[c("ID",names(x)[grep(pattern="^Ratio L/H",names(x))])]
  # Strip down fraction names
  names(tmpRaw) <- sub(pattern="Ratio L/H Exp. ",replacement="",names(tmpRaw))
  names(tmpRel) <- sub(pattern="Rel. Ratio L/H Exp. ",replacement="",names(tmpRel))
  # Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
  tmpRaw <- melt(data=tmpRaw,measure.vars=LETTERS[1:13],na.rm=TRUE)
  tmpRaw$variable <- match(tmpRaw$variable,LETTERS)
  tmpRel <- melt(data=tmpRel,measure.vars=LETTERS[1:13],na.rm=TRUE)
  tmpRel$variable <- match(tmpRel$variable,LETTERS)
  # Prep the flippase Activity analogously
  tmpFA <- flippaseActivity
  tmpFA[["Av. Spec. Activity"]] <- Relativate(tmpFA[["Av. Spec. Activity"]])
  names(tmpFA) <- c("variable","value")
  tmpFA$variable <- match(tmpFA$variable,LETTERS)
  # Plot A
  rawPlot <- ggplot(data=tmpRaw,aes(x=variable,y=value))
  rawPlot <- rawPlot +
    geom_line() +
    facet_wrap(~ID) +
#     geom_line(data=tmpFA,col="red") +
    labs(
      x="Fraction",
      y="Ratio L/H")
  output[["A"]] <- rawPlot
  # Plot B
  relPlot <- ggplot(data=tmpRel,aes(x=variable,y=value))
  relPlot <- relPlot +
    geom_line() +
    facet_wrap(~ID) +
    geom_line(data=tmpFA,col="red") +
    labs(
      x="Fraction",
      y="Relative Ratio L/H/Relative Flippase Activity")
  output[["B"]] <- relPlot
  return(invisible(output))
}