#' @name checkAA
#' 
#' @title Check if amino acid sequence contains right characters
#' 
#' @description 
#' The function \code{checkAA} checks if the letters in \code{aa} contain the 
#' correct characters (all letters in \code{aa} have to be in 
#' \code{AA_ALPHABET}). The function will return the sequence as a 
#' character vector or will break otherwise.
#'
#' @details
#' \code{aa} has to be an \code{AAString} object.
#' 
#' @param aa \code{AAString} containing the amino acid, peptide, or protein 
#' sequence 
#' 
#' @return character
#' 
#' @importFrom methods is
#' @importFrom S4Vectors safeExplode
#' @importFrom Biostrings AA_ALPHABET
#' 
#' @examples 
#' library(Biostrings)
#' aa <- AAString(x = "TEST")
#' PhysicoChemicalPropertiesProtein:::checkAA(aa)
checkAA <- function(aa) {
    
    if (!is(aa, "AAString"))
        stop("'aa' has to be an 'AAString' object")
    
    ## convert aa to character
    aa <- as.character(aa)
    
    ## strsplit aa and check if the letters are in AA_ALPHABET
    aa_split <- S4Vectors::safeExplode(aa)
    
    if (!all(aa_split %in% Biostrings::AA_ALPHABET))
        stop("'aa' contains letters not defined in the amino acid alphabet")
    
    return(aa)
}