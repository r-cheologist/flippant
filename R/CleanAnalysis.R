directory = "~/localtmp/Anant/"
# Load Data
protData <- pGroups
# pepData <- peptides

dataStats <- data.frame(
  Parameter = c("Raw Count"),
#   Peptides=nrow(pepData),
  Proteins=nrow(protData),
  stringsAsFactors=FALSE)

# Extract columns containing NOT-NORMALIZED ratios and assemble
protSubsetter <- grep(pattern="^Ratio H/L Exp. [[:alpha:]]{1}$",names(protData))
# pepSubsetter <- grep(pattern="^Ratio H/L Exp. [[:alpha:]]{1}$",names(pepData))

# Reverse
protCols <- sub(pattern="H/L",replacement="L/H",names(protData)[protSubsetter])
protData[protCols] <- RatioInversion(protData[protSubsetter])
rm(protSubsetter)
# pepCols <- sub(pattern="H/L",replacement="L/H",names(pepData)[pepSubsetter])
# pepData[pepCols] <- RatioInversion(pepData[pepSubsetter])
# rm(pepSubsetter)

# Logarithmize
protData[protCols] <- log2(protData[protCols])
names(protData)[which(names(protData) %in% protCols)] <- paste("Log2",protCols)
protCols <- paste("Log2",protCols)

# pepData[pepCols] <- log2(pepData[pepCols])
# names(pepData)[which(names(pepData) %in% pepCols)] <- paste("Log2",pepCols)
# pepCols <- paste("Log2",pepCols)
# 
# # Home grown ratio averaging per protein and fraction
# ## Isolate/split pGroups?
# tmpPGroups <- strsplit(x=pepData[["Protein group IDs"]],split=";")
# ## How many pGroups per entry?
# tmpMultiplicity <- sapply(tmpPGroups,length)
# ## Replicate rows associated with multiple entries
# pepData <- pepData[rep(x=seq(nrow(pepData)),times=tmpMultiplicity),]
# pepData[["Protein group IDs"]] <- unlist(tmpPGroups,use.names=FALSE)
# ## Split by pGroup
# pepDataList <- split(x=pepData,f=pepData[["Protein group IDs"]])
# ## Medianize logarithmized, inverted ratios by protein & fraction
# pepDataList <- lapply(
#   pepDataList,
#   function(x){
#     id <- x[1,"Protein group IDs"]
#     pepratios <- x[pepCols]
#     medians <- apply(pepratios,2,median,na.rm=TRUE)
#     names(medians) <- paste("Median Peptide",names(medians))
#     cbind(id,data.frame(rbind(medians),check.names=FALSE))
#   })
# pepData <- rbind.fill(pepDataList)
# ## Merge with protein data
# protData <- merge(x=protData,y=pepData,by="id",all=TRUE)
# ## Explore MQ protein averaging vs home grown one ...
# protCols2 <- paste("Median Peptide",protCols)
# protData[["CC.MQ.vs.HG.PGroup.Ratios"]] <- apply(
#   X=protData,
#   MARGIN=1,
#   FUN=function(x){
#     cor(x=as.numeric(x[protCols]),y=as.numeric(x[protCols2]),use="na.or.complete",method="spearman")
#   })
# annotation_df <- data.frame(
#   Median = median(protData[["CC.MQ.vs.HG.PGroup.Ratios"]],na.rm=TRUE),
#   MAD = mad(protData[["CC.MQ.vs.HG.PGroup.Ratios"]],na.rm=TRUE))
# tmpPlot <- ggplot(data=protData)
# tmpPlot +
#   geom_rect(data=annotation_df,aes(xmax=Median+MAD,xmin=Median-MAD,ymin=-Inf,ymax=Inf),alpha=0.25) +
#   geom_vline(data=annotation_df,aes(xintercept=Median),col="Red") +
#   geom_density(aes(x=CC.MQ.vs.HG.PGroup.Ratios)) +
#   geom_text(data=annotation_df,aes(x=Median,y=Inf,label=paste("Median:",round(Median,2))),hjust=1.05,vjust=1.5) +
#   labs(
#     x="CC (Spearman) MaxQuant vs. Manual Peptide Ratio Averaging per Fraction and Protein",
#     y="Density")
## --> It does not seem to be worth doing this.

# Relativate
protDataList <- lapply(
  seq(nrow(protData)),
  function(x){
    ratios <- unlist(protData[x,protCols])
    data.frame(rbind(Relativate(ratios)),check.names=FALSE)
  })
protDataList <- rbind.fill(protDataList)
protData[paste("Rel.",protCols)] <- protDataList
# names(protData)[which(names(protData) %in% protCols)] <- paste("Rel.",protCols)
protCols <- paste("Rel.",protCols)

# Discard reverse and contaminant
protData <- protData[protData$Reverse != "+",]
protData <- protData[protData$Contaminant != "+",]
dataStats <- rbind(
  dataStats,
  data.frame(
    Parameter="NoRevCont",
    Proteins=nrow(protData),
    stringsAsFactors=FALSE))

# Filter for minimum number of values
minValueCount <- (length(protCols) - 1)
protData <- protData[rowSums(!is.na(protData[protCols])) >= minValueCount,]
dataStats <- rbind.fill(
  dataStats,
  data.frame(
    Parameter=paste("Min",minValueCount,"Values",sep=""),
    Proteins=nrow(protData),
    stringsAsFactors=FALSE))

# Gradient QC
protData[["MaxFraction"]] <- apply(X=protData[protCols],MARGIN=1,FUN=function(x){which(x=="1")[1]})
tmpPlot <- ggplot(data=data.frame(MaxFraction=factor(protData$MaxFraction),MW=protData[["Mol. weight [kDa]"]]))
tmpPlot +
  geom_boxplot(aes_string(x="MaxFraction",y="MW",group="MaxFraction"),outlier.colour=NA) +
  geom_jitter(aes_string(x="MaxFraction",y="MW",group="MaxFraction")) +
  labs(
    x="Glycerol Gradient Fraction",
    y="Molecular Weight (kDa)")
## --> seems ok to fraction 6, but then?

# Compare protein profiles with Flippase Activity using multiple measures
relLog2Fa <- Relativate(log2(flippaseActivity[["Av. Spec. Activity"]]))
# Graumann's choice
protData[["Spearman.CC"]] <- Correlate(
  x=protData[protCols],
  y=relLog2Fa,
  use.rows=TRUE,
  method="spearman",
  use="na.or.complete")
# # Original Andersen et al. 2003
# protData[["ChiSquare"]] <- ChiSquare(
#   x=protData[protCols],
#   y=relLog2Fa,
#   use.rows=TRUE)# tmpPlot <- ggplot(data=protData)
# tmpPlot <- ggplot(data=protData)
# tmpPlot + 
#   geom_point(aes_string(x="Spearman.CC",y="ChiSquare")) +
#   labs(
#     x="Correlation Coefficient (Spearman)",
#     y=expression(chi^2))
### --> Spearman correlates reasonably with ChiSquare and is easier to evaluate
### --> Sticking with it
## Recommended by Dengjel et al 2010. 
## With this particular variance/covariance matrix the same as ChiSquare ...
# protData[["Mahalanobis"]] <- sapply(
#   seq(nrow(protData)),
#   function(x){
#     protratios <- protData[x,protCols]
#     mahalanobis(x=protratios,center=relFa,cov=diag(length(protratios)))
#   })
# tmpPlot <- ggplot(data=protData)
# tmpPlot + geom_point(aes(x=ChiSquare,y=Mahalanobis))

# # Investigate effect of normalization on CC
# protData[["Unrel.Spearman.CC"]] <- Correlate(
#   x=protData[sub(pattern="^Rel. ",replacement="",x=protCols)],
#   y=log2(flippaseActivity[["Av. Spec. Activity"]]),
#   use.rows=TRUE,
#   method="spearman",
#   use="na.or.complete")
# tmpPlot <- ggplot(data=protData)
# tmpPlot + 
#   geom_point(aes_string(x="Unrel.Spearman.CC",y="Spearman.CC"))

# Add CC bins
protData$Spearman.CC.Bin <- cut(
  x=protData$Spearman.CC,
  breaks=seq(from=-1,to=1,by=0.05),
  include.lowest=TRUE)

# Retrieve Biomart/ENSEMBL/TMHMM annotation regarding potential transmembrane
# nature
ensembl <- useMart(biomart="ensembl",dataset="scerevisiae_gene_ensembl")
leadingMajorUniprotIds <- sapply(
  strsplit(protData[["Majority protein IDs"]],";"),
  function(x){x[1]})
tmTable <- getBM(
  mart=ensembl,
  attributes=c("uniprot_swissprot_accession","transmembrane_domain"),
  filters="uniprot_swissprot_accession",
  values=leadingMajorUniprotIds)
protData <- cbind(
  protData,
  transmembrane_domain=tmTable[match(x=leadingMajorUniprotIds,table=tmTable[[1]]),"transmembrane_domain"])
# # Save
# write.table(
#   x=protData,
#   file=file.path(directory,paste(Sys.Date(),"Correlation.txt",sep="_")),
#   sep="\t",
#   row.names=FALSE,
#   col.names=TRUE)

# # Are tranmembrane proteins per se flippases?
# tmpPlotData <- protData
# tmpPlotData$transmembrane_domain <- factor("All")
# tmpPlotData <- rbind.fill(tmpPlotData,protData[protData$transmembrane_domain != "",])
# tmpPlot <- ggplot(data=rbind.fill(tmpPlotData))
# tmpPlot +
#   stat_density(
#     aes_string(
#       x="Spearman.CC",
#       color="transmembrane_domain",
#       size="..count.."),
#     geom="line",
#     position="identity") +
#   labs(
#     x="Correlation Coefficient (Spearman)",
#     y="Density",
#     color="Protein Groups",
#     size="Group Members")
## --> does not seem to be the case!

# Tag peptide count/raw ratio onto label
tmpLabels  <- sapply(strsplit(protData[["Gene names"]],";"),function(x){x[1]})
tmpLabels[is.na(tmpLabels)] <- protData[is.na(tmpLabels),"Majority protein IDs"]
protData$Label <- paste(
  tmpLabels,
  round(protData$Spearman.CC,2),
  protData[["Peptides Exp. D"]],
  round(RatioInversion(protData[["Ratio H/L Exp. D"]]),2),
  protData$transmembrane_domain,
  sep="|")

# Prep for plotting
plotData <- protData
names(plotData) <- sub(pattern="Rel. Log2 Ratio L/H Exp. ",replacement="",names(plotData))

# Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
plotData <- melt(data=plotData,measure.vars=LETTERS[1:13],na.rm=TRUE)
plotData$variable <- match(plotData$variable,LETTERS)
relLog2Fa <- data.frame(rbind(relLog2Fa))
names(relLog2Fa) <- LETTERS[1:length(relLog2Fa)]
relLog2Fa <- melt(data=relLog2Fa,measure.vars=LETTERS[1:13],na.rm=TRUE)
relLog2Fa$variable <- match(relLog2Fa$variable,LETTERS)

tmpPlot <- ggplot(data=plotData)
tmpPlot + 
  geom_line(aes(x=variable,y=value,group=Label)) +
  geom_line(data=relLog2Fa,aes(x=variable,y=value),color="red") +
  facet_wrap(~Spearman.CC.Bin,nrow=8) +
  labs(
    x="Fraction",
    y="Relative Logarithmized Ratio L/H, Relative Logarithmized Flippase Activity",
    title="Spearman Correlation Coefficient Bins")

tmpPlot <- ggplot(data=plotData[plotData$transmembrane_domain == "TMHMM",])
tmpPlot + 
  geom_line(aes(x=variable,y=value,group=Label)) +
  geom_line(data=relLog2Fa,aes(x=variable,y=value),color="red") +
  facet_wrap(~Spearman.CC.Bin,nrow=8) +
  labs(
    x="Fraction",
    y="Relative Logarithmized Ratio L/H, Relative Logarithmized Flippase Activity",
    title="Spearman Correlation Coefficient Bins - TMD Proteins Only")

for(bin in levels(plotData$Spearman.CC.Bin)){
  tmpPlotData <- plotData[plotData$Spearman.CC.Bin == bin,]
  tmpPlot <- ggplot(data=tmpPlotData)
  tmpPlot <- tmpPlot +
    geom_line(aes(x=variable,y=value,color=Label))+
    geom_line(data=relLog2Fa,aes(x=variable,y=value),color="red") +
    labs(
      x="Fraction",
      y="Relative Logarithmized Ratio L/H, Relative Logarithmized Flippase Activity",
      title=paste("Spearman Correlation Coefficient Bin",bin))
  if(!(bin %in% c("[-1,-0.95]","(-0.95,-0.9]","(-0.9,-0.85]","(0.15,0.2]","(0.35,0.4]","(0.4,0.45]","(0.95,1]"))){
    if(bin %in% c("(-0.65,-0.6]","(-0.6,-0.55]","(-0.55,-0.5]","(-0.4,-0.35]")){
      tmpPlot <- tmpPlot + guides(col=guide_legend(ncol=3))
    } else {
      tmpPlot <- tmpPlot + guides(col=guide_legend(ncol=2))
    }
  }
#   print(tmpPlot)
  tmpBin <- sub("\\[","",bin)
  tmpBin <- sub("\\]","",tmpBin)
  tmpBin <- sub("\\(","",tmpBin)
  tmpBin <- sub("\\)","",tmpBin)
  tmpBin <- sub("\\,","_",tmpBin)
  ggsave(
    filename=file.path(directory,paste("SpearmanCC",tmpBin,"Bin.pdf",sep="_")),
    width=11.69,
    height=8.27,
    unit="in",
    dpi=300)
#   readline("Hit <ENTER> to proceed ...")
}

#######################################

plotData <- protData[protData$transmembrane_domain !="",]
plotData <- plotData[!is.na(plotData$transmembrane_domain),]
pairwisecombinations <- combn(seq(nrow(plotData)),m=2,simplify=FALSE)
pairwiselist <- lapply(
  pairwisecombinations,
  function(x){
    # NA + x currently produces NA!
    ratios <- data.frame(rbind(plotData[x[1],protCols] + plotData[x[2],protCols]))
    names(ratios) <- paste("Comb.",protCols)
    ratios[["Label"]] <- paste(sub(pattern="\\|.*$",replacement="",x=plotData[x,"Label"]),collapse=" + ")
    return(ratios)
  })
pairwise <- rbind.fill(pairwiselist)
# Careful! Based on strange current version of relLog2Fa!
pairwise[["Spearman.CC"]] <- Correlate(
  x=pairwise[paste("Comb.",protCols)],
  y=relLog2Fa$value,
  use.rows=TRUE,
  method="spearman",
  use="na.or.complete")
# Add CC bins
pairwise$Spearman.CC.Bin <- cut(
  x=pairwise$Spearman.CC,
  breaks=seq(from=-1,to=1,by=0.05),
  include.lowest=TRUE)
table(pairwise$Spearman.CC.Bin)
## --> There are approx 100 pairs just in the 0.95-1 bin
## --> likely not much power in this analysis.

# Prep for plotting
plotData <- pairwise
names(plotData) <- sub(pattern="Comb. Rel. Log2 Ratio L/H Exp. ",replacement="",names(plotData))

# Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
plotData <- melt(data=plotData,measure.vars=LETTERS[1:13],na.rm=TRUE)
plotData$variable <- match(plotData$variable,LETTERS)
table()

tmpPlot <- ggplot(data=plotData[plotData$Spearman.CC.Bin == "(0.95,1]",])
tmpPlot +
  geom_line(aes(x=variable,y=value,color=Label)

#######################################

# Reload data - no data point filtration
protData <- pGroups
# Extract columns containing NOT-NORMALIZED ratios and assemble
protSubsetter <- grep(pattern="^Ratio H/L Exp. [[:alpha:]]{1}$",names(protData))
# Reverse
protCols <- sub(pattern="H/L",replacement="L/H",names(protData)[protSubsetter])
protData[protCols] <- RatioInversion(protData[protSubsetter])
rm(protSubsetter)
# Logarithmize
protData[protCols] <- log2(protData[protCols])
names(protData)[which(names(protData) %in% protCols)] <- paste("Log2",protCols)
protCols <- paste("Log2",protCols)
# Relativate
protDataList <- lapply(
  seq(nrow(protData)),
  function(x){
    ratios <- unlist(protData[x,protCols])
    data.frame(rbind(Relativate(ratios)),check.names=FALSE)
  })
protDataList <- rbind.fill(protDataList)
protData[paste("Rel.",protCols)] <- protDataList
# names(protData)[which(names(protData) %in% protCols)] <- paste("Rel.",protCols)
protCols <- paste("Rel.",protCols)
# Discard reverse and contaminant
protData <- protData[protData$Reverse != "+",]
protData <- protData[protData$Contaminant != "+",]
dataStats <- rbind(
  dataStats,
  data.frame(
    Parameter="NoRevCont",
    Proteins=nrow(protData),
    stringsAsFactors=FALSE))
# Compare protein profiles with Flippase Activity using multiple measures
relLog2Fa <- Relativate(log2(flippaseActivity[["Av. Spec. Activity"]]))
# Graumann's choice
protData[["Spearman.CC"]] <- Correlate(
  x=protData[protCols],
  y=relLog2Fa,
  use.rows=TRUE,
  method="spearman",
  use="na.or.complete")
# Add CC bins
protData$Spearman.CC.Bin <- cut(
  x=protData$Spearman.CC,
  breaks=seq(from=-1,to=1,by=0.05),
  include.lowest=TRUE)
# Retrieve Biomart/ENSEMBL/TMHMM annotation regarding potential transmembrane
# nature
ensembl <- useMart(biomart="ensembl",dataset="scerevisiae_gene_ensembl")
leadingMajorUniprotIds <- sapply(
  strsplit(protData[["Majority protein IDs"]],";"),
  function(x){x[1]})
tmTable <- getBM(
  mart=ensembl,
  attributes=c("uniprot_swissprot_accession","transmembrane_domain"),
  filters="uniprot_swissprot_accession",
  values=leadingMajorUniprotIds)
protData <- cbind(
  protData,
  transmembrane_domain=tmTable[match(x=leadingMajorUniprotIds,table=tmTable[[1]]),"transmembrane_domain"])
# Tag peptide count/raw ratio onto label
tmpLabels  <- sapply(strsplit(protData[["Gene names"]],";"),function(x){x[1]})
tmpLabels[is.na(tmpLabels)] <- protData[is.na(tmpLabels),"Majority protein IDs"]
protData$Label <- paste(
  tmpLabels,
  round(protData$Spearman.CC,2),
  protData[["Peptides Exp. D"]],
  round(RatioInversion(protData[["Ratio H/L Exp. D"]]),2),
  protData$transmembrane_domain,
  sep="|")

# Find Peaks
peakData <- protData[grep(pattern="Rel. Log2 Ratio L/H Exp. ",x=names(protData))]
rownames(peakData) <- protData$id
# Identify peaks
peaks <- lapply(
  seq(nrow(peakData)),
  function(x){
    output <- list(
      Peaks = MaxQuant2dPeakDetection(unlist(peakData[x,]),trim=FALSE,return.half=TRUE),
      id=rownames(peakData[x,]))
    return(output)
  })

# Crosscorrelate peaks (>=3 points!)
peaks <- lapply(
  seq(length(peaks)),
  function(x){
    tmpPeaks <- peaks[[x]]$Peaks
    tmpCC <- sapply(
      tmpPeaks,
      function(y){
        subsetter <- seq(from=y$Start,to=y$Stop)
        if(length(subsetter) < 3){return(NA)}
        tmpRatios <- peakData[x,subsetter]
        tmpFA <- relLog2Fa[subsetter]
        return(
          Correlate(
            x=tmpRatios,
            y=tmpFA,
            use.rows=TRUE,
            method="spearman"))
      })
    tmpLength <- sapply(
      tmpPeaks,
      function(x){
        length(seq(from=x$Start,to=x$Stop))
      })
    output <- list(id = peaks[[x]]$id, Peaks = tmpPeaks, CC = tmpCC, PPP = tmpLength)
    return(output)
  })
# Isolate the best peak
bestPeak <- lapply(
  seq(length(peaks)),
  function(x){
    y <- peaks[[x]]
    if(sum(is.na(y$CC))==length(y$Peaks)){
      output <- data.frame(
        PeakStart=NA,
        PeakStop=NA,
        PeakMax=NA,
        PeakLength=NA,
        Peak.CC=NA)          
    } else {
      subsetter <- unname(which.max(y$CC))
      output <- data.frame(
        PeakStart=y$Peaks[[subsetter]][["Start"]],
        PeakStop=y$Peaks[[subsetter]][["Stop"]],
        PeakMax=y$Peaks[[subsetter]][["Max"]],
        PeakLength=length(seq(from=y$Peaks[[subsetter]][["Start"]],to=y$Peaks[[subsetter]][["Stop"]])),
        Peak.CC=y$CC[[subsetter]]
      )
    }
    output$id <- y$id
    return(output)
  })
bestPeak <- rbind.fill(bestPeak)
# Merge with Data
protData <- merge(x=protData,y=bestPeak,by="id",all=TRUE)
# Toss NAs
protData <- protData[!is.na(protData$Peak.CC),]
# Is the peak-derived CC better than the one derived from the total data?
protData <- protData[protData$Spearman.CC < protData$Peak.CC,]
# Is the peak-derived CC >= 0.8
protData <- protData[protData$Peak.CC >= 0.8,]
# Is the PPP > 3?
protData <- protData[protData$PeakLength > 3,]
# Is the peak max where we expect it?
protData <- protData[protData$PeakMax %in% c(3,4,5),]
# Is the total CC <= 0.8?
protData <- protData[protData$Spearman.CC <= 0.8,]
# Tag peptide count/raw ratio onto label
tmpLabels  <- sapply(strsplit(protData[["Gene names"]],";"),function(x){x[1]})
tmpLabels[is.na(tmpLabels)] <- protData[is.na(tmpLabels),"Majority protein IDs"]
protData$Label <- paste(
  tmpLabels,
  round(protData$Peak.CC,2),
  protData[["Peptides Exp. D"]],
  round(RatioInversion(protData[["Ratio H/L Exp. D"]]),2),
  protData$transmembrane_domain,
  sep="|")
# Order by CC
protData <- protData[order(protData$Peak.CC,decreasing=TRUE),]
# Blotting batches
panels <- ceiling(nrow(protData)/20)
protData$Panel <- rep(1:panels,each=20)[1:nrow(protData)]
# Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
tmpRel <- protData
names(tmpRel) <- sub(pattern="Rel. Log2 Ratio L/H Exp. ",replacement="",names(tmpRel))
tmpRel <- melt(data=tmpRel,measure.vars=LETTERS[1:13])
tmpRel$variable <- match(tmpRel$variable,LETTERS)
# Prep the flippase Activity analogously
relLog2Fa <- data.frame(rbind(Relativate(log2(flippaseActivity[["Av. Spec. Activity"]]))))
names(relLog2Fa) <- LETTERS[1:length(relLog2Fa)]
relLog2Fa <- melt(data=relLog2Fa,measure.vars=LETTERS[1:13],na.rm=TRUE)
relLog2Fa$variable <- match(relLog2Fa$variable,LETTERS)
for(panel in seq(panels)){
  tmpPlot <- ggplot(data=tmpRel[tmpRel$Panel == panel,],aes(x=factor(variable),y=value))
  tmpPlot <- tmpPlot +
    geom_rect(aes(xmin=factor(PeakStart),xmax=factor(PeakStop),ymin=-Inf,ymax=Inf),fill="blue",alpha=0.01) +
    geom_line(aes(group=Label)) +
    geom_line(data=relLog2Fa,color="red",aes(group=NA)) +
    facet_wrap(~Label) +
    labs(
      x="Fraction",
      y="Relative Logarithmized Ratio L/H/Relative Flippase Activity")
  ggsave(
    filename=file.path(directory,paste("PeakCorrelation_",panel,".pdf",sep="")),
    plot=tmpPlot,
    width=11.69,
    height=8.27,
    units="in")
}







# Prep for plotting
plotData <- protData
names(plotData) <- sub(pattern="Rel. Log2 Ratio L/H Exp. ",replacement="",names(plotData))

# Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
plotData <- melt(data=plotData,measure.vars=LETTERS[1:13],na.rm=TRUE)
plotData$variable <- match(plotData$variable,LETTERS)
relLog2Fa <- data.frame(rbind(relLog2Fa))
names(relLog2Fa) <- LETTERS[1:length(relLog2Fa)]
relLog2Fa <- melt(data=relLog2Fa,measure.vars=LETTERS[1:13],na.rm=TRUE)
relLog2Fa$variable <- match(relLog2Fa$variable,LETTERS)

tmpPlot <- ggplot(data=plotData)
tmpPlot + 
  geom_line(aes(x=variable,y=value,group=Label)) +
  geom_line(data=relLog2Fa,aes(x=variable,y=value),color="red") +
  facet_wrap(~Spearman.CC.Bin,nrow=8) +
  labs(
    x="Fraction",
    y="Relative Logarithmized Ratio L/H, Relative Logarithmized Flippase Activity",
    title="Spearman Correlation Coefficient Bins")