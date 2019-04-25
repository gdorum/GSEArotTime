#' Restricted maximum likelihood to estimate the covariance matrix of longitudinal expression data with gene dependencies within gene sets and correlation between samples
#' 
#' Restricted maximum likelihood for longitudinal expression data where genes are assumed correlated within gene set and uncorrelated between different gene sets. 
#' All genes within the same gene set are assumed to have identical covariance matrix. The gene sets may differ in size and hence have different covariance matrix.
#' @param alphaStart A numeric vector of start values for the parameters to be estimated. The first parameter must be gene variance, the second parameter random error 
#' variance, the third must be a parameter controlling correlation between time points and the fourth must be a parameter controlling correlation between genes. 
#' The remaining parameters can be variances for other random design factors.
#' @param Xlist Design matrix with fixed design variables for a gene set (assuming same design for all gene sets)
#' @param Ylist vector of observed values for all gene sets
#' @param Zlist Design matrix for random design factors for a gene set (assuming same design for all gene sets)
#' @param TsList Structure of time dependencies between samples in a gene set (assuming same design for all gene sets)
#' @param GsList Structure of gene dependencies between samples in a gene set (assuming same design for all gene sets)
#' @param BsList (Optional) structure of batch dependencies between samples in a gene set (assuming same design for all gene sets)
#' @param A Restrictions for parameters. See \code{\link{constrOptim}}
#' @param b Restrictions for parameters. See \code{\link{constrOptim}}
#' @return \item{alpha.hat }{Vector of estimates for the parameters in alphaStart}
#' \item{beta.hat }{Vector of estimated parameters for fixed design factors}
#' \item{sigma }{Estimated common variance}
#' \item{maxlik }{The maximum likelihood value}
#' @details Arguments \code{A} and \code{b} corresponds to arguments \code{ui} and \code{ci} in \code{\link{constrOptim}}.
#' @seealso \code{\link{reml}} for restricted maximum likelihood without gene correlations and \code{\link{reml_genedep}} that does not assume identical gene set size.
#' @importFrom stats constrOptim 
reml_genedep <- function(alphaStart, Xlist, Ylist, Zlist, TsList, GsList, BsList, A, b) {
  
  p <- ncol(Xlist[[1]])             #Number of fixed factors
  nset <- length(Xlist)
  
  res <- constrOptim( theta=alphaStart, f=loglik_genedep, grad=NULL, ui=A, ci=b, X=Xlist, Z=Zlist, Y=Ylist, Ts=TsList, Gs=GsList, Bs=BsList)
  alpha.hat <- res$par
  maxlik <- res$value
  
  phi.hat <- alpha.hat[3]
  ga.hat <- alpha.hat[4]
  #Time variance, gene variance, random error variance, remaining variances
  theta.hat <- c(1,alpha.hat[-c(3,4)])
  
  RSS <- numeric(nset) 
  N <- 0
  for(k in 1:nset) {
    
    Y <- Ylist[[k]]
    X <- Xlist[[k]]
    Z <- Zlist[[k]]
    Bs <- BsList[[k]]
    Gs <- GsList[[k]]
    Ts <- TsList[[k]]
    
    n <- length(Y)
    dim(Y) <- c(n,1)
    N <- N + n #Counting the total number of samples
    
    #Time and gene correlations
    Rt <- exp(-phi.hat*Ts)
    Rg <- exp(-ga.hat*Gs)
    
    #Make covariance matrix
    if( ncol(Z) > 2 ) {
      #only works when we have a batch effect. Else use makecov()
      V.hat <- theta.hat[4]*Bs + theta.hat[1]*Rt + theta.hat[2]*Rg + diag(n)*theta.hat[3]
    }else V.hat <- theta.hat[1]*Rt + theta.hat[2]*Rg + diag(n)*theta.hat[3]#V.hat <- makecov_genedep(theta.hat,Rt=Rt,Rg=Rg,design=Z)  
    
    Vinv <- solve(V.hat)
    XVinv <- t(X)%*%Vinv
    XVinvX <- XVinv%*%X
    XV <- solve(XVinvX)%*%XVinv
    
    beta.hat <- XV%*%Y
    r <- Y - X%*%beta.hat
    
    RSS[k] <- t(r)%*%Vinv%*%r        
  }    
  sigma.hat <- sum(RSS) / (N-nset*p)
  
  return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat, maxlik=maxlik))
  
}