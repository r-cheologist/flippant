#' @title dithionite_flippase_assay_stats
#' @description A function providing tau values for dithionite flippase assays
#' @details Input format for and data processing by the function is documented 
#' at \code{\link{dithionite_flippase_assay_plot}}.
#' @param x \code{\link{data.frame}} as described in \code{\link{dithionite_flippase_assay_plot}}.
#' @param scale_to Defines the source of \code{ymax}, defaulting to 
#' \code{model}. See \code{\link{dithionite_flippase_assay_plot}}.
#' @return Returns a \code{\link{data.frame}}.
#' @author Johannes Graumann
#' @references Menon, I., Huber, T., Sanyal, S., Banerjee, S., Barré, P., Canis, 
#' S., Warren, J.D., Hwa, J., Sakmar, T.P., and Menon, A.K. (2011). Opsin Is a 
#' Phospholipid Flippase. Current Biology 21, 149–153.
#' MIKE PAPER
#' @export
#' @seealso \code{\link{ParseQuantMasterData}}, \code{\link{ParseLS55Data}}, 
#' \code{\link{dithionite_flippase_assay_plot}},
#' \code{\link{dithionite_flippase_assay_input_validation}},
#' \code{\link{dithionite_flippase_assay_calculations}}
#' @keywords methods manip
#' @import plyr
#' @import robustbase
#' @examples
#' stop("Add citation to Mike's manuscript!")
#' stop("Revisit workflow description.")
#' stop("Add example using actually published data.")
#' # Build input
#' x <- data.frame(
#'  Path = c(
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-plus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-minus-150ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/ePC.txt",
#'  "~/localTmp/Fluor Data_Menon Lab/21FEB2013_Erg1 immun_deple/Erg1 TE-plus-15ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-40ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-75ul.txt",
#'  "~/localTmp/Fluor Data_Menon Lab//21FEB2013_Erg1 immun_deple/Erg1 TE-plus-150ul.txt"),
#'  "Extract Volume (ul)" = c(0,15,40,75,150,0,15,40,75,150,0,15,40,75,150,0,15,40,75,150),
#'  #     "Reaction Volume w/o DT (ul)" = rep(2000,4),
#'  "Reaction Volume with DT (ul)" = rep(2040,20),
#'  "Concentration Egg PC (mM)" = rep(4.5,20),
#'  "Extract Protein Concentration (mg/ml)" = c(rep(0.67,5),rep(1.26,5),rep(0.67,5),rep(1.26,5)),
#'  #     "Timepoint of Measurement (s)",
#'  "Experimental Series"=c(rep("Extract",5),rep("Depleted Extract",5),rep("Extract",5),rep("Depleted Extract",5)),
#'  Experiment=c(rep("Erg1, Replicate 1",10),rep("Erg1, Replicate 2",10)),
#'  check.names=FALSE,
#'  stringsAsFactors=FALSE)
#'  # Run function
#'  DithioniteFlippaseAssay(x)
dithionite_flippase_assay_stats <- function(x,scale_to=c("model","data")){
# Check Prerequisites -----------------------------------------------------
  validated_params <- dithionite_flippase_assay_input_validation(x=x,scale_to=scale_to)
  x <- validated_params[["x"]]
  scale_to <- validated_params[["scale_to"]]

# Processing --------------------------------------------------------------
  processed_list_from_x <- dithionite_flippase_assay_calculations(x=x,scale_to=scale_to)
  processed_list_from_x <- lapply(processed_list_from_x,function(y){y[["Raw"]]})

# Assemble the output -----------------------------------------
  output <- rbind.fill(processed_list_from_x)
  output <- output[!duplicated(output$CombinedId),]
  output <- output[which(names(output) %in% c("Fit Constant (a)","PPR at P = 0.5","Experimental Series","Experiment"))]

# Return ------------------------------------------------------------------
  return(output)
}