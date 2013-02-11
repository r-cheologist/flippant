#' @export
CompleteAnalysis <- function(directory=tempdir()){
  ###########################################
  # Profile Correlation Analysis - Figure 1 #
  ###########################################
  # Prep the data for correlation analysis
  fig1Data <- DataPrep(pGroups)
  # Correlate (Spearman)
  fig1Data$Inverted$CC <- Correlate(
    x=fig1Data[["Inverted"]],
    y=flippaseActivity[["Av. Spec. Activity"]])
  fig1Data$RelativeInverted[["Rel.CC"]] <- Correlate(
    x=fig1Data[["RelativeInverted"]],
    y=Relativate(flippaseActivity[["Av. Spec. Activity"]]))
  # Merge with original data
  subsetter <- sapply(rownames(fig1Data$RelativeInverted),function(x){which(pGroups$id == x)})
  fig1Data <- cbind(pGroups[subsetter,],fig1Data$Inverted,fig1Data$RelativeInverted)
  # Retrieve Biomart/ENSEMBL/TMHMM annotation regarding potential transmembrane
  # nature
  ensembl <- useMart(biomart="ensembl",dataset="scerevisiae_gene_ensembl")
  leadingMajorUniprotIds <- sapply(
    strsplit(fig1Data[["Majority protein IDs"]],";"),
    function(x){x[1]})
  tmTable <- getBM(
    mart=ensembl,
    attributes=c("uniprot_swissprot_accession","transmembrane_domain"),
    filters="uniprot_swissprot_accession",
    values=leadingMajorUniprotIds)
  fig1Data <- cbind(fig1Data,transmembrane_domain=tmTable[match(x=leadingMajorUniprotIds,table=tmTable[[1]]),"transmembrane_domain"])
  # Save
  write.table(
    x=fig1Data,
    file=file.path(directory,paste(Sys.Date(),"Correlation.txt",sep="_")),
    sep="\t",
    row.names=FALSE,
    col.names=TRUE)
  # Investigate effect of normalization on CC
#   sFigure1 <- ggplot(data=fig1Data)
#   sFigure1 + 
#     geom_point(aes_string(x="Rel.CC",y="CC"))
  # Produce first figure
  Figure1(fig1Data,directory=directory)
  ############################################################
  # Simple sorting: which ratio profile peaks in fraction 4? #
  ############################################################
  sS <- SimpleSort(Fraction=4)
  write.table(
    x = fig1Data[na.omit(match(sS$id,fig1Data$id)),],
    file = file.path(
      directory,
      paste(
        Sys.Date(),
        "MaxRatioInFrac4.txt",
        sep="_")
    ),
    sep = "\t",row.names=FALSE,col.names=TRUE)
  sS <- SimpleSort(Fraction=3)
  write.table(
    x = fig1Data[na.omit(match(sS$id,fig1Data$id)),],
    file = file.path(
      directory,
      paste(
        Sys.Date(),
        "MaxRatioInFrac3.txt",
        sep="_")
    ),
    sep = "\t",row.names=FALSE,col.names=TRUE)
  ###########################################################################
  # Fish out profiles for flippase candidates tested negatively (Menon lab) #
  ###########################################################################
#   NonFlippases <- c("Cpt1","Ept1","Ist2","Cho1","Pis1","Ale1","Cho2","Slc1","Yop1","Srp102","Pbn1","Sec12","Ost1","Tsc10","Wbp1","Spc3","Pga3","Sop4","Sur2","Lcb2","Pmt1")
  Figure2 <- FilteredProfilePlots(fig1Data,NonFlippases)
  ####################################
  # Fish out profiles for Erg1/Erg11 #
  ####################################
  Figure3 <- FilteredProfilePlots(fig1Data,c("Erg1","Erg11"))
  # p-type ATPases known to flip lipids
  Figure4 <- FilteredProfilePlots(fig1Data,c("Dfn1","Dfn2","Dfn3","Drs2","Neo1"))
  # Organellar markers for plasma membrane (Pma1,Gas1,Can1),Vacuole (Vph1), Golgi (Mnn1)
  Figure5 <- FilteredProfilePlots(fig1Data,c("Pma1","Gas1","Can1","Vph1","Mnn1","Sec61"))
  # Enzymes/transporters involved in protein N-glycosylation
  Figure6 <- FilteredProfilePlots(
    fig1Data,
    c("Alg7","Alg13","Alg14","Alg1","Alg2","Alg11","Rft1","Alg3","Alg9","Alg12",
      "Alg6","Alg8","Alg10","Alg5","Dpm1","Dfg10","Sec59","Ost1","Swp1","Wbp1",
      "Ost4","Stt3","Ost3","Ost6","Ost2"))
  # Enzymes of GPI anchor biosynthesis
  Figure7 <- FilteredProfilePlots(
    fig1Data,
    c("Gpi3","Gpi15","Gpi2","Gpi1","Gpi19","Eri1","Gpi12","Gwt1","Gpi14","Pbn1",
      "Gpi18","Mcd4","Gpi10","Gpi13","Gpi11","Gpi7","Gaa1","Gpi8","Gpi16",
      "Gpi17","Gab1","Bst1"))
  # Religious list of abundant proteins
  Figure8 <- FilteredProfilePlots(
    fig1Data,
    c("ERG1","PHO88","SUR2","LCB2","ERV25","HOM6","YMR122W-A","SAC1","PMT1",
      "IFA38","HXT3","RTN1","ERG3","YNR021W","PGA3","ERP2","SEC61","ERG5",
      "TSC13","VMA13","GSF2","TPO4","LCB1","ERP1","SNA2","STE24","FAT1","SEC62",
      "WBP1","YHL017W","YNL181W","DPL1","VOA1","OST1","SOP4","ARV1"))
  FilteredProfilePlots(fig1Data,c("STE24"))
  # Cooperativity Erg1/11?
  Figure6 <- AdditiveFilteredProfilePlots(fig1Data,c("Erg1","Erg11"))
  ####################################
  # Do correlation analysis on PEAKS #
  ####################################
  # Prep the data for correlation analysis
  fig2Data <- fig1Data[grep(pattern="Rel. Ratio L/H Exp. ",x=names(fig1Data))]
  rownames(fig2Data) <- fig1Data$id
  # Identify peaks
  peaks <- lapply(
    seq(nrow(fig2Data)),
    function(x){
      output <- list(
        Peaks = MaxQuant2dPeakDetection(unlist(fig2Data[x,]),trim=FALSE,return.half=TRUE),
        id=rownames(fig2Data[x,]))
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
          tmpRatios <- fig2Data[x,subsetter]
          tmpFA <- Relativate(flippaseActivity[["Av. Spec. Activity"]])[subsetter]
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
    peaks,
    function(x){
      if(sum(is.na(x$CC))==length(x$Peaks)){
        output <- data.frame(
          PeakStart=NA,
          PeakStop=NA,
          PeakMax=NA,
          PeakLength=NA,
          Peak.CC=NA)          
      } else {
        subsetter <- which.max(x$CC)
        output <- data.frame(
          PeakStart=x$Peaks[[subsetter]][["Start"]],
          PeakStop=x$Peaks[[subsetter]][["Stop"]],
          PeakMax=x$Peaks[[subsetter]][["Max"]],
          PeakLength=length(seq(from=x$Peaks[[subsetter]][["Start"]],to=x$Peaks[[subsetter]][["Stop"]])),
          Peak.CC=x$CC[[subsetter]]
          )
      }
      output$id <- x$id
      return(output)
    })
  bestPeak <- rbind.fill(bestPeak)
  # Merge with Data
  fig2Data <- merge(x=fig1Data,y=bestPeak,by="id",all=FALSE)
  # Toss NAs
  fig2Data <- fig2Data[!is.na(fig2Data$Peak.CC),]
  # Is the peak-derived CC better than the one derived from the total data?
  fig2Data <- fig2Data[fig2Data$CC < fig2Data$Peak.CC,]
  # Is the peak-derived CC >= 0.8
  fig2Data <- fig2Data[fig2Data$Peak.CC >= 0.8,]
  # Is the PPP > 3?
  fig2Data <- fig2Data[fig2Data$PeakLength > 3,]
  # Is the peak max where we expect it?
  fig2Data <- fig2Data[fig2Data$PeakMax %in% c(3,4,5),]
  # Is the total CC <= 0.8?
  fig2Data <- fig2Data[fig2Data$CC <= 0.8,]
  # Fill in missing gene names with Uniprot IDs
  subsetter <- fig2Data[["Gene names"]]==""
  fig2Data[subsetter,"Gene names"] <- fig2Data[subsetter,"Majority protein IDs"]
  names(fig2Data)[which(names(fig2Data)=="Gene names")] <- "ID"
  # Tag peptide count/raw ratio onto label
  fig2Data$ID <- sapply(
      seq(nrow(fig2Data)),
      function(x){
        output <- paste(
          fig2Data[x,"ID"],
          round(fig2Data[x,"Peak.CC"],2),
          fig2Data[x,"Peptides Exp. D"],
          round(fig2Data[x,"Ratio L/H Exp. D"],2),
          fig2Data[x,"transmembrane_domain"],
          sep="|")
        return(output)
      })
  # Order by CC
  fig2Data <- fig2Data[order(fig2Data$Peak.CC,decreasing=TRUE),]
  # Blotting batches
  panels <- ceiling(nrow(fig2Data)/20)
  fig2Data$Panel <- rep(1:panels,each=20)[1:nrow(fig2Data)]
  # Melt the data frame into plottable shape and swtich fraction names from alpha to numeric
  tmpRel <- fig2Data
  names(tmpRel) <- sub(pattern="Rel. Ratio L/H Exp. ",replacement="",names(tmpRel))
  tmpRel <- melt(data=tmpRel,measure.vars=LETTERS[1:13])
  tmpRel$variable <- match(tmpRel$variable,LETTERS)
  # Prep the flippase Activity analogously
  tmpFA <- flippaseActivity
  tmpFA[["Av. Spec. Activity"]] <- Relativate(tmpFA[["Av. Spec. Activity"]])
  names(tmpFA) <- c("variable","value")
  tmpFA$variable <- match(tmpFA$variable,LETTERS)
  for(panel in seq(panels)){
    tmpPlot <- ggplot(data=tmpRel[tmpRel$Panel == panel,],aes(x=variable,y=value))
    tmpPlot <- tmpPlot +
      geom_rect(aes(xmin=PeakStart,xmax=PeakStop,ymin=-Inf,ymax=Inf),fill="blue",alpha=0.01) +
      geom_line() +
      geom_line(data=tmpFA,col="red") +
      facet_wrap(~ID) +
      labs(
        x="Fraction",
        y="Relative Ratio L/H/Relative Flippase Activity")
    ggsave(
      filename=file.path("","tmp",paste("PeakCorrelation_",panel,".pdf",sep="")),
      plot=tmpPlot,
      width=11.69,
      height=8.27,
      units="in")
  }
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