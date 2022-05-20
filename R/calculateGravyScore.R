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
    
    ## check the AAString object that it contains the valid residues
    aa <- checkAA(aa)
    
    ## load the object that contains hydropathy values and retrieve scores
    f <- system.file("GRAVY/hydropathy.RDS", 
        package = "PhysicoChemicalPropertiesProtein")
    mappings <- readRDS(f)
    score <- mappings[, "SCORE", drop = FALSE]
    
    ## get all the characters of the sequence
    aa_split <- strsplit(aa, split = "")[[1]]
    
    ## get the corresponding scores for the sequence characters and calculate
    ## the mean of the sequence
    aa_score <- score[aa_split, ]
    mean(aa_score)
}
