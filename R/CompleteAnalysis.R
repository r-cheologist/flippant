#' @export
CompleteAnalysis <- function(dir=tempdir()){
  ###########################################
  # Profile Correlation Analysis - Figure 1 #
  ###########################################
  # Prep the data for correlation analysis
  fig1Data <- DataPrep(pGroups)
  # Correlate (Spearman)
  fig1Data$RelativeInverted$CC <- Correlate(
    x=fig1Data[["RelativeInverted"]],
    y=Relativate(flippaseActivity[["Av. Spec. Activity"]]))
  # Merge with original data
  subsetter <- sapply(rownames(fig1Data$RelativeInverted),function(x){which(pGroups$id == x)})
  fig1Data <- cbind(pGroups[subsetter,],fig1Data$Inverted,fig1Data$RelativeInverted)
  # Retrieve Biomart/ENSEMBL/TMHMM annotation regarding potential transmembrane
  # nature
  library(biomaRt)
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
    file=file.path(dir,paste(Sys.Date(),"Correlation.txt",sep="_")),
    sep="\t",
    row.names=FALSE,
    col.names=TRUE)
  # Produce first figure
  Figure1(fig1Data,dir=dir)
  ############################################################
  # Simple sorting: which ratio profile peaks in fraction 4? #
  ############################################################
  sS <- SimpleSort(Fraction=4)
  write.table(
    x = fig1Data[na.omit(match(sS$id,fig1Data$id)),],
    file = file.path(
      dir,
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
      dir,
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
  # Cooperativity Erg1/11?
  Figure6 <- AdditiveFilteredProfilePlots(fig1Data,c("Erg1","Erg11"))
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