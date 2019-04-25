#' QR decomposition
#' 
#' Function that performs a QR decomposition on a square matrix.
#' @param X Square matrix (nxn)
#' @return  \item{Q}{Orientation matrix}
#' \item{R}{Configuration matrix}
#' @author Solve Sabo
#' @export
QR <-
function(X){ 

    n <- dim(X)[1]
    if(dim(X)[2]!=n)stop("Applies only to quadratic matrices \n")

    Ident <- diag(n)
    Q <- Ident
    R <- X
    Emark <- X
    for(j in 1:(n-1)){
        Emark <- R[j:n,j:n]
        Qmark <- householder(Emark[,1])
        Q1 <- Ident
        Q1[j:n,j:n] <- Qmark
        R <- Q1%*%R    
        Q <- Q%*%Q1
    }

    return(list(Q=Q,R=R))
}
