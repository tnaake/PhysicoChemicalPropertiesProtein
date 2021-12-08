#' @name calculateIsoelectricPoint
#' 
#' @title Calculate isoelectric point of an amino acid/protein
#' 
#' @description 
#' The function `calculateIsoelectricPoint` calculates the isoelectric point 
#' of a given amino acid/protein sequence. 
#' Several methods are implemented that use different pKa values for the 
#' ionizable groups of proteins (N terminal, C terminal, and amino acids
#' C, D, E, H, K, R, and Y). 
#'
#' @details
#' The pKa values are taken from Table 4 in Kozlowski (2016).
#' 
#' For polypeptides, the isoelectric point depends primarily on the dissociation 
#' constants (pKa) for the ionizable groups of seven charged amino acids: 
#' glutamate (\\delta-carboxyl group), aspartate (\\beta-carboxyl group), cysteine 
#' (thiol group), tyrosine (phenol group), histidine (imidazole side chains), 
#' lysine (\\epsilon-ammonium group) and arginine (guanidinium group). Moreover, 
#' the charge of the terminal groups (NH2 and COOH) can greatly affect the pI 
#' of short peptides. Generally, the Glu, Asp, Cys, and Tyr ionizable groups 
#' are uncharged below their pKa and negatively charged above their pKa. 
#' Similarly, the His, Lys, and Arg ionizable groups are positively charged 
#' below their pKa and uncharged above their pKa.
#' 
#' The implemented algorithm (bisection algorithm) is as in Kozlowski (2016).
#' 
#' @references
#' Kozlowski (2016): 'IPC - Isoelectric Point Calculator', Biology Direct, 
#' 11, 55.
#' doi: doi.org/10.1186/s13062-016-0159-9
#' 
#' @param aa `character(1)`, vector containing the amino acid, peptide, or 
#' protein sequence 
#' @param method `character(1)`, vector specifying the method/parametrization
#' for the pKa values
#' 
#' @export
#' 
#' @examples 
#' aa <- "TEST"
#' calculateIsoelectricPoint(aa = aa, method = "EMBOSS")
calculateIsoelectricPoint <- function(aa, 
    method = c("EMBOSS", "DTASelect", "Solomon", "Sillero", "Rodwell", 
        "Lehninger", "Toseland", "Thurlkill", "Nozaki", 
        "IPC_protein", "IPC_peptide")) {
    
    method <- match.arg(method)
    
    if (method == "EMBOSS")
        pKa <- c(NH2 = 8.6, COOH = 3.6, C = 8.5, D = 3.9, E = 4.1, 
                 H = 6.5, K = 10.8, R = 12.5, Y = 10.1)
    
    if (method == "DTASelect") 
        pKa <- c(NH2 = 8, COOH = 3.1, C = 8.5, D = 4.4, E = 4.4, 
                 H = 6.5, K = 10, R = 12, Y = 10)
    
    if (method == "Solomon")
        pKa <- c(NH2 = 9.6, COOH = 2.4, C = 8.3, D = 3.9, E = 4.3, 
                 H = 6, K = 10.5, R = 12.5, Y = 10.1)
    
    if (method == "Sillero")
        pKa <- c(NH2 = 8.2, COOH = 3.2, C = 9, D = 4, E = 4.5, 
                 H = 6.4, K = 10.4, R = 12, Y = 10)
    
    if (method == "Rodwell")
        pKa <- c(NH2 = 8, COOH = 3.1, C = 8.33, D = 3.68, E = 4.25, 
                 H = 6, K = 11.5, R = 11.5, Y = 10.07)
    
    if (method == "Lehninger")
        pKa <- c(NH2 = 9.69, COOH = 2.34, C = 8.33, D = 3.86, E = 4.25, 
                 H = 6, K = 10.5, R = 12.4, Y = 10)
    
    if (method == "Toseland")
        pKa <- c(NH2 = 8.71, COOH = 3.19, C = 6.87, D = 3.6, E = 4.29, 
                 H = 6.33, K = 10.45, R = 12, Y = 9.61)
    
    if (method == "Thurlkill")
        pKa <- c(NH2 = 8, COOH = 3.67, C = 8.55, D = 3.67, E = 4.25, 
                 H = 6.54, K = 10.4, R = 12, Y = 9.84)
    
    if (method == "Nozaki")
        pKa <- c(NH2 = 7.5, COOH = 3.8, C = 9.5, D = 4, E = 4.4, 
                 H = 6.3, K = 10.4, R = 12, Y = 9.6)
    
    if (method == "IPC_protein")
        pKa <- c(NH2 = 9.094, COOH = 2.869, C = 7.555, D = 3.872, E = 4.412, 
                 H = 5.637, K = 9.052, R = 11.84, Y = 10.85)
    
    if (method == "IPC_peptide")
        pKa <- c(NH2 = 9.564, COOH = 2.383, C = 8.297, D = 3.887, E = 4.317, 
                 H = 6.018, K = 10.517, R = 12.503, Y = 10.071)
    
    ## get all the characters of the sequence
    aa_split <- strsplit(aa, split = "")[[1]]
    residues <- c("C", "D", "E", "H", "K", "R", "Y")
    num <- table(aa_split)[residues]
    num <- num[!is.na(names(num))]
    num_missing <- residues[!residues %in% names(num)]
    if (length(num_missing) > 0) {
        num_missing_vec <- numeric(length(num_missing))
        names(num_missing_vec) <- num_missing 
        num <- c(num, num_missing_vec)
    }
    num <- c("NH2" = 1, "COOH" = 1, num)
    
    ## the net charge of the peptide/protein is related to the solution pH,
    ## use the Henderson-Hasselbalch equation to calculate the charge at a 
    ## certain pH
    pH <- 7
    pHprev <- 0
    pHnext <- 14
    pI <- FALSE
    
    ## implement the bisection algorithm, which in each iteration halves
    ## the search space and then moves higher or lower by 3.5 (half of 7)
    ## depending on the charge. In the next ieration, the pH is changed
    ## by 1.75 (half of 3.5), ... The process is repeated until the algorithm
    ## reaches the desired precision
    while (!pI) {
        ## for negatively charged residues
        ## the Glu (E), Asp (D), Cys (C), and Tyr (Y) ionizable groups are 
        ## uncharged below their pKa and negatively charged above their pKa,
        ## charge of COOH group: negative
        
        sum_COOH <- sum_E <- sum_D <- sum_C <- sum_Y <- 0
        
        ## calculate the sum of the Henderson-Hasselbalch equation along the
        ## residues
        sum_COOH <- -num[["COOH"]] / (1 + 10^(pKa[["COOH"]] - pH))
        sum_E <- -num[["E"]] / (1 + 10^(pKa[["E"]] - pH))
        sum_D <- -num[["D"]] / (1 + 10^(pKa[["D"]] - pH))
        sum_C <- -num[["C"]] / (1 + 10^(pKa[["C"]] - pH))
        sum_Y <- -num[["Y"]] / (1 + 10^(pKa[["Y"]] - pH))
        sum_neg <- sum_COOH + sum_E + sum_D + sum_C + sum_Y
        
        ## for positively charged residues
        ## His (H), Lys (K), and Arg (R) are positively charged below their pKa  
        ## and uncharged above their pKa. 
        ## charge of NH2 group: positive
        sum_NH2 <- sum_H <- sum_K <- sum_R <- 0
        
        ## calculate the sum of the Henderson-Hasselbalch equation along the 
        ## residues
        sum_NH2 <- num[["NH2"]] / (1 + 10^(pH - pKa[["NH2"]]))
        sum_H <- num[["H"]] / (1 + 10^(pH - pKa[["H"]]))
        sum_K <- num[["K"]] / (1 + 10^(pH - pKa[["K"]]))
        sum_R <- num[["R"]] / (1 + 10^(pH - pKa[["R"]]))
        sum_pos <- sum_NH2 + sum_H + sum_K + sum_R
        
        ## charge of a macromolecule at a given pH is the sum of the positive
        ## and negative charges of the individual amino acids 
        if (sum_pos + sum_neg < 0) {
            ## the new pH value must be smaller
            pH_tmp <- pH
            pH <- pH - (pH - pHprev) / 2
            pHnext <- pH_tmp
        } else {
            ## the new pH value must be higher
            pH_tmp <- pH
            pH <- pH + (pHnext - pH) / 2
            pHprev <- pH_tmp
        }
        
        if (pH - pHprev < 0.001 & pHnext - pH < 0.001) 
            pI <- TRUE
        
        if (pH > 14) stop("'pH' is greater than 14")
        if (pH < 0) stop("'pH' is smaller than 14")
    }
    
    return(pH)
}
