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