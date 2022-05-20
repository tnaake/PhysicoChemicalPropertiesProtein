#' @name calculateGravyScore
#' 
#' @title Calculate Gravy scores of an amino acid/protein
#' 
#' @description 
#' The function \code{calculateGravyScore} calculates the GRAVY (grand average 
#' of hydropathy) value of a given amino acid/protein sequence. 
#' The GRAVY value is calculated by adding the hydropathy value for each 
#' residue and dividing by the length of the sequence.
#'
#' @details
#' The hydropathy values are taken from Table 2 in Kyte and Doolittle (1982).
#' 
#' @references
#' Kyte and Doolittle (1982):
#' 'A Simple Method for Displaying the  Hydropathic Character of a Protein',
#' Journal of Molecular Biology, 175, 105-132.
#' doi: 10.1016/0022-2836(82)90515-0
#' 
#' @param aa \code{AAString} containing the amino acid, peptide, or 
#' protein sequence 
#' 
#' @export
#' 
#' @examples 
#' library(Biostrings)
#' aa <- AAString(x = "TEST")
#' calculateGravyScore(aa = aa)
calculateGravyScore <- function(aa) {
    
    aa <- checkAA(aa)
    
    mappings <- rbind(
        c("I", "Ile", "isoleucine", 4.5),
        c("V", "Val", "valine", 4.2),
        c("L", "Leu", "leucine", 3.8),
        c("F", "Phe", "phenylalanine", 2.8),
        c("C", "Cys", "cysteine", 2.5),
        c("M", "Met", "methionine", 1.9),
        c("A", "Ala", "alanine", 1.8),
        c("G", "Gly", "glycine", -0.4),
        c("T", "Thr", "threonine", -0.7),
        c("W", "Trp", "tryptophan", -0.9),
        c("S", "Ser", "serine", -0.8),
        c("Y", "Tyr", "tyrosine", -1.3),
        c("P", "Pro", "proline", -1.6),
        c("H", "His", "histidine", -3.2),
        c("E", "Glu", "glutamic acid", -3.5),
        c("Q", "Gln", "glutamine", -3.5),
        c("D", "Asp", "aspartic acid", -3.5),
        c("N", "Asn", "asparagine", -3.5),
        c("K", "Lys", "lysine", -3.9),
        c("R", "Arg", "arginine", -4.5))
    
    mappings <- as.data.frame(mappings)
    colnames(mappings) <- c("SINGLE_LETTER", "THREE_LETTER", "NAME", "SCORE")
    rownames(mappings) <- mappings$SINGLE_LETTER
    mappings$SCORE <- as.numeric(mappings$SCORE)
    score <- mappings[, "SCORE", drop = FALSE]
    
    ## get all the characters of the sequence
    aa_split <- strsplit(aa, split = "")[[1]]
    
    ## get the corresponding scores for the sequence characters and calculate
    ## the mean of the sequence
    aa_score <- score[aa_split, ]
    mean(aa_score)
}
