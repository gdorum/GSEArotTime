#' Restricted maximum likelihood to estimate the covariance matrix of longitudinal expression data with correlation between samples
#' 
#' Restricted maximum likelihood for longitudinal expression data where genes are assumed independent. All genes are assumed to have identical covariance matrix.
#' @param alphaStart A numeric vector of start values for the parameters to be estimated. The first parameter must be random error variance, the last must be a 
#' parameter controlling correlation between time points. The remaining parameters can be variances for other random design factors 
#' @param X Design matrix with fixed design variables
#' @param Y Matrix of observed values. Rows are samples and columns are variables.
#' @param Z Design matrix for random design factors
#' @param subject A numeric vector indicating which subject each sample belongs to
#' @param A Restrictions for parameters. See \code{\link{constrOptim}}
#' @param b Restrictions for parameters. See \code{\link{constrOptim}}
#' @return \item{alpha.hat }{Vector of estimates for the parameters in alphaStart}
#' \item{beta.hat }{Vector of estimated parameters for fixed design factors}
#' \item{sigma }{Estimated common variance}
#' \item{maxlik }{The maximum likelihood value}
#' @details Arguments \code{A} and \code{b} corresponds to arguments \code{ui} and \code{ci} in \code{\link{constrOptim}}.
#' @seealso \code{\link{reml2}} and \code{\link{reml_genedep}} for restricted maximum likelihood where gene correlations are assumed within gene set.
#' @importFrom stats constrOptim 
#' @export
reml <-
function(alphaStart, X, Y, Z, subject, A, b) {
  
  G <- ncol(Y) #Number of genes
  N <- nrow(Y) #Number of samples
  p <- ncol(X) #Number of fixed factors
  
  res <- constrOptim( theta=alphaStart, f=loglik, grad=NULL, ui=A, ci=b, X=X, Z=Z, Y=Y, subject=subject)
  alpha.hat <- res$par
  maxlik <- res$value
  
  phi <- alpha.hat[length(alpha.hat)]
  theta.hat <- c(1,alpha.hat[-length(alpha.hat)])
  
  #Make covariance matrix with estimated parameters
  #Find correlation between time points
  R <- matrix(0,N,N)
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      if(subject[i]==subject[j]) R[i,j] <- exp( -phi*( Z[j,1] - Z[i,1] ) )
    }
  }
  R <- R + t(R)
  diag(R) <- 1
  #Find covariance matrix
  V.hat <- makecov(theta.hat,correlation=R,design=Z,subject=subject)
  
  #Estimate beta and sigma for final values of alpha
  Vinv <- solve(V.hat)
  XVinv <- t(X)%*%Vinv
  XVinvX <- XVinv%*%X
  beta.hat <- solve(XVinvX)%*%XVinv%*%Y
  r <- Y - X%*%beta.hat
  RSS <- sum( diag( t(r)%*%Vinv%*%(r) ) )
  sigma.hat <- RSS / ( G*(N-p) )
  
  return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat, maxlik=maxlik))
  
}
