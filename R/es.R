#' Enrichment score
#' 
#' Calculates an enrichment score for a gene set.
#' @param L Vector with one expression value per gene.
#' @param S Vector of same length as L where 1 indicates that the gene in L is present in the gene set and 1 indicates that it is not.
#' @param p Weight. Default is 1.
#' @param doplot If TRUE, the running sum is plotted. Default is FALSE.
#' @param index If TRUE, the index at which the ES occurs in the sorted list is returned. Default is FALSE.
#' @details See Subramanian et al. for details.
#' @return  \item{E }{Enrichment score.}
#' \item{ind }{Index of the enrichment score (if \code{index=TRUE})}.
#' @references Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#' Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
#' and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
#' approach for interpreting genome-wide expression profiles, \emph{PNAS}, \bold{102},15545-15550.
#' @author Solve Sabo, Guro Dorum
#' @importFrom graphics abline hist plot
#' @export
es <-
function(L,S,p=1, doplot=FALSE, index=FALSE) {   
    #Sort L and S according to L
    L_sort <- sort(L,decreasing=TRUE,index.return=TRUE)
    L <- L_sort$x
    S <- S[L_sort$ix]

    Sc <- 1-S
    Ns <- sum(S)
    N <- length(L)
    
    #If L and S have none or all genes in common
    if( Ns == 0 )
        stop("No genes are member of the gene set")
    if( Ns == length(S) )
        stop("All genes are members of the gene set")
        
    #Weighting factor (N_R in Subramanian et al., 2005)    
    Ws <- sum(S*abs(L)^p)
    
    pmiss <- cumsum(Sc)/(N-Ns)
    phit <- cumsum(S*abs(L)^p)/Ws
    
    #Running sum
    ph_pm <- phit-pmiss
    
    #The enrichment score is the maximum deviation from 0 of the running sum
    ind <- which.max(abs(ph_pm))    
    E <- ph_pm[ind]
    names(E) <- NULL 
    
    #Plot running sum?
    if(doplot) 
    {
        plot(1:N,ph_pm,"l",col=2,lwd=2,xlab="L",ylab="Phit-Pmiss",main="Running sum")
        abline(h=0)
        abline(v=ind,lty=3)
    }
    
    if(index) return(list(E=E, ind=ind)) #Also return index of occurence of ES
    else return(E) 

}
