#' @title 20121111_proteinGroups
#' @description "proteinGroups.txt" as produced for a protein correlation 
#' profiling experiment in search of lipid flippses in the endoplasmic 
#' reticulum.
#' @details
#' The data set derives from 156 RAW files corresponding to 13 velocity gradient
#' fractions from a seperation of flippase activity that was subjected to in-gel
#' digest using 12 gel slices each (\eqn{13\times12 = 156}).
#' 
#' Prior to GeLC/MS analysis, each fractions has been mixed 1:1 by protein 
#' content with material equivalent to the input material of the velocity 
#' gradient but metabolically (SILAC, see references below) labeled using 
#' "heavy" Lysine (Lysine 8; \eqn{^13C_6 ^15N_2 L-lysine}).
#' 
#' The raw data was analyze using MaxQuant v.1.3.0.3 against the 
#' \emph{Saccharmocyes cerevisiae} reference proteome from UniProtKB (including 
#' isoform data; downloaded on) using the following parameters:
#' \itemize{
#'  \item Ratio.H.L (SILAC) ratio as calculated by MaxQuant
#'  \item Ratio.H.L.Log2 Logarithmized version of the above
#'  \item Intensity Signal intensity as summed by MaxQuant (binning dimension for 
#'    Significance B calculation)
#'  \item Ratio.H.L.Log2.Significance.A Significance A as calculated by Perseus 1.1.1.17
#'    from the logarithmized ratios
#'  \item Ratio.H.L.Log2.Significance.B Significance B as calculated by Perseus 1.1.1.17
#'    from the logarithmized ratios and corresponding intensity values}
#' @author Johannes Graumann, Roopesh Krisnankutty
#' @format A tab-delimited text file as produced by MaxQuant.
#' @keywords datasets
#' @docType data
#' @name data_20121111_proteinGroups
NULL
