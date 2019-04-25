#' Significance level computation in GSEA
#' 
#'Computes a p-value for a gene set.
#'@param ES Observed enrichment score.
#'@param ESnull Null distribution.
#'@param doplot If \code{TRUE}, a histogram over the normalised null distribution is drawn. Default is \code{FALSE}.
#' @details Function that computes a p-value for a gene set given an enrichment score 
#' and a null distribution. Calculations are done separately for positive and
#' negative enrichment scores. The function also returns normalised observed 
#' enrichment score and normalised null distribution for the gene set. 
#' In addition, a histogram over the normalised null distribution can be plotted.
#' @return \item{p.value}{p-value for the gene set}
#'   \item{NES}{Normalised observed enrichment score}
#'  \item{NESnull}{Normalised estimated null distribution}
#' @author Guro Dorum
#' @importFrom graphics hist
#' @export
significance <-
function(ES, ESnull, doplot=FALSE) { 
    
    #Include observed ES in null distribution to prevent p-values of magnitude 0
    ESnull <- c(ES,ESnull)
    #Calculate p-value separately for positive (>= 0) and negative ES as the 
    #number of ES's in the null distribution as least as extreme as the observed ES
    if(ES >= 0) {
        Nl <- sum(ESnull >= ES)
        pos <- ESnull[ESnull >= 0]
        p.value <- Nl/length(pos)
        NES <- ES/mean(pos)
        NESnull <- ESnull/mean(pos)
    } else {
        Ns <- sum(ESnull <= ES)
        neg <- ESnull[ESnull < 0]
        p.value <- Ns/length(neg)
        NES <- ES/abs(mean(neg))
        NESnull <- ESnull/abs(mean(neg))
    }
    #Remove observed NES from null distribution
    NESnull <- NESnull[-1]
    if(doplot)   
        hist(NESnull,freq=FALSE,xlab="NES",ylab="Density",main="Null distribution for NES")
    
    return(list(p.value=p.value, NES=NES, NESnull=NESnull))
}
