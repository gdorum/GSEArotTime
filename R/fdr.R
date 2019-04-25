#' False discovery rate computation for GSEA.
#' 
#' Computes a false discovery rate (FDR) q-value for a gene
#' set given the observed normalised enrichment scores and a normalised null
#' distribution. Calculations are done separately for positive and negative
#' normalised enrichment scores.
#' @param NES Normalised enrichment score for the given gene set.
#' @param NESobs Vector of normalised enrichment scores for all gene sets to be tested.
#' @param NESnull Matrix of normalised null distribution for all gene sets to be tested. Each column represent a gene set, and each row represent a permutation.
#' @return FDR q-value for the given gene set.
#' @references Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#' Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
#' and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
#' approach for interpreting genome-wide expression pro?les, \emph{PNAS}, \bold{102},
#' 15545-15550.
#' @author Guro Dorum
#' @export
fdr <-
function(NES,NESobs,NESnull) { 
    nrot <- nrow(NESnull)
    
    Nnp <- numeric(nrot)
    Nnl <- numeric(nrot)
    Nnn <- numeric(nrot)
    Nns <- numeric(nrot)

    #For each permutation, counting the number of positive/negative NES and the number of NES more 
    #extreme than the given NES in the null distribution
    for(i in 1:nrot) {
               
        if(NES >= 0) {
            Nnp[i] <- sum(NESnull[i,] >= 0) + 1
            Nnl[i] <- sum(NESnull[i,] >= NES) + 1
        } else {
           Nnn[i] <- sum(NESnull[i,] < 0) + 1
           Nns[i] <- sum(NESnull[i,] <= NES) + 1
        }
    }
    
    #Counting the number of positive/negative NES and NES more extreme than the given NES,
    #among the NES of all gene sets to be tested. Finally, calculating the q-value.
    if(NES >= 0) {
    
        Np <- sum(NESobs >= 0)
        Nl <- sum(NESobs >= NES)
        Nnpl <- Nnl/Nnp 
        q <- mean(Nnpl)/(Nl/Np)
        
    } else {
        
        Nn <- sum(NESobs < 0)
        Ns <- sum(NESobs <= NES)
        Nnns <- Nns/Nnn
        q <- mean(Nnns)/(Ns/Nn)
        
    }
    if( q > 1) q <- 1
    q
}
