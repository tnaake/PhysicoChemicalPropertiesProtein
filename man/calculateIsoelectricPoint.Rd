\name{calculateIsoelectricPoint}
\alias{calculateIsoelectricPoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateIsoelectricPoint(aa, method = c("EMBOSS", "DTASelect", "Solomon", "Sillero", "Rodwell", "Lehninger", "Toseland", "Thurlkill", "Nozaki", "IPC_protein", "IPC_peptide"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{
%%     ~~Describe \code{aa} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (aa, method = c("EMBOSS", "DTASelect", "Solomon", "Sillero", 
    "Rodwell", "Lehninger", "Toseland", "Thurlkill", "Nozaki", 
    "IPC_protein", "IPC_peptide")) 
{
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
        pKa <- c(NH2 = 9.69, COOH = 2.34, C = 8.33, D = 3.86, 
            E = 4.25, H = 6, K = 10.5, R = 12.4, Y = 10)
    if (method == "Toseland") 
        pKa <- c(NH2 = 8.71, COOH = 3.19, C = 6.87, D = 3.6, 
            E = 4.29, H = 6.33, K = 10.45, R = 12, Y = 9.61)
    if (method == "Thurlkill") 
        pKa <- c(NH2 = 8, COOH = 3.67, C = 8.55, D = 3.67, E = 4.25, 
            H = 6.54, K = 10.4, R = 12, Y = 9.84)
    if (method == "Nozaki") 
        pKa <- c(NH2 = 7.5, COOH = 3.8, C = 9.5, D = 4, E = 4.4, 
            H = 6.3, K = 10.4, R = 12, Y = 9.6)
    if (method == "IPC_protein") 
        pKa <- c(NH2 = 9.094, COOH = 2.869, C = 7.555, D = 3.872, 
            E = 4.412, H = 5.637, K = 9.052, R = 11.84, Y = 10.85)
    if (method == "IPC_peptide") 
        pKa <- c(NH2 = 9.564, COOH = 2.383, C = 8.297, D = 3.887, 
            E = 4.317, H = 6.018, K = 10.517, R = 12.503, Y = 10.071)
    aa_split <- strsplit(aa, split = "")[[1]]
    residues <- c("C", "D", "E", "H", "K", "R", "Y")
    num <- table(aa_split)[residues]
    num <- num[!is.na(names(num))]
    num_missing <- residues[!residues \%in\% names(num)]
    if (length(num_missing) > 0) {
        num_missing_vec <- numeric(length(num_missing))
        names(num_missing_vec) <- num_missing
        num <- c(num, num_missing_vec)
    }
    num <- c(NH2 = 1, COOH = 1, num)
    pH <- 7
    pHprev <- 0
    pHnext <- 14
    pI <- FALSE
    while (!pI) {
        sum_COOH <- sum_E <- sum_D <- sum_C <- sum_Y <- 0
        sum_COOH <- -num[["COOH"]]/(1 + 10^(pKa[["COOH"]] - pH))
        sum_E <- -num[["E"]]/(1 + 10^(pKa[["E"]] - pH))
        sum_D <- -num[["D"]]/(1 + 10^(pKa[["D"]] - pH))
        sum_C <- -num[["C"]]/(1 + 10^(pKa[["C"]] - pH))
        sum_Y <- -num[["Y"]]/(1 + 10^(pKa[["Y"]] - pH))
        sum_neg <- sum_COOH + sum_E + sum_D + sum_C + sum_Y
        sum_NH2 <- sum_H <- sum_K <- sum_R <- 0
        sum_NH2 <- num[["NH2"]]/(1 + 10^(pH - pKa[["NH2"]]))
        sum_H <- num[["H"]]/(1 + 10^(pH - pKa[["H"]]))
        sum_K <- num[["K"]]/(1 + 10^(pH - pKa[["K"]]))
        sum_R <- num[["R"]]/(1 + 10^(pH - pKa[["R"]]))
        sum_pos <- sum_NH2 + sum_H + sum_K + sum_R
        if (sum_pos + sum_neg < 0) {
            pH_tmp <- pH
            pH <- pH - (pH - pHprev)/2
            pHnext <- pH_tmp
        }
        else {
            pH_tmp <- pH
            pH <- pH + (pHnext - pH)/2
            pHprev <- pH_tmp
        }
        if (pH - pHprev < 0.001 & pHnext - pH < 0.001) 
            pI <- TRUE
        if (pH > 14) 
            stop("'pH' is greater than 14")
        if (pH < 0) 
            stop("'pH' is smaller than 14")
    }
    return(pH)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
