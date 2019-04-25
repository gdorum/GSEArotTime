loglik <-
  function(theta, X, Z, Y, subject) {
    #log-likelihood function to be optimised with reml().
    #
    #Input:
    #    - theta: a vector of paramterers to be optimised over. The first parameter must be random error variance, the last
    #             must be a parameter controlling correlation between time points. The remaining parameters can be variances
    #             for other random design factors. 
    #    - X: design matrix with fixed design variables.
    #    - Z: design matrix with random design factors.
    #    - Y: matrix of observed expression data. Rows are samples and columns are variables.
    #    - subject: vector of indices indicating subject for each sample.
    #Output:
    #    - maxLoglik: maximum likelihood
    
    G <- ncol(Y)
    N <- nrow(Y)
    p <- ncol(X)
    
    phi <- theta[length(theta)]
    alpha <- c(1,theta[-length(theta)])
    
    R <- matrix(0,N,N)
    for(i in 1:(N-1)) {
      for(j in (i+1):N) {
        if(subject[i]==subject[j]) R[i,j] <- exp( -phi*( Z[j,1] - Z[i,1] ) )
      }
    }
    R <- R + t(R)
    diag(R) <- 1
    
    V.hat <- makecov(alpha,correlation=R,design=Z,subject=subject)   
    
    Vinv <- solve(V.hat)
    XVinv <- t(X)%*%Vinv
    XVinvX <- XVinv%*%X
    #beta.hat are parameters for fixed effects and different for all genes
    beta.hat <- solve(XVinvX)%*%XVinv%*%Y
    r <- Y - X%*%beta.hat
    #Gives a GxG symmetric matrix in which we want the numbers on the diagonal
    RSS <- sum( diag( t(r)%*%Vinv%*%(r) ) )
    sigmat.hat <- RSS / ( G*(N-p) )
    
    #Add extra minus before the log-likelihood because R optimisation functions find minimum on default. 
    maxLoglik <- - ( - 0.5*( G*(N-p)*log(sigmat.hat) + G*log( det(V.hat) ) + RSS/sigmat.hat + G*log( det( XVinvX ) ) ) )
    
    maxLoglik     
    
  }


loglik2 <-
  function(theta, X, Z, Y, Ts, Gs, Bs, set) {
    #log-likelihood function to be optimised with reml().
    #
    #Input:
    #    - theta: a vector of parameters to be optimised over. The first parameter must be gene variance, 
    #      the second parameter random error variance, the third must be a parameter controlling correlation between time points and the fourth must
    #      be a parameter controlling correlation between genes. The remaining parameters can be variances for other random design factors. 
    #    - X: design matrix with fixed design variables for a gene set (assuming same design for all gene sets).
    #    - Z: design matrix for random design factors for a gene set (assuming same design for all gene sets).
    #    - Y: vector of observed values for all gene sets.
    #    - Ts: structure of time dependencies between samples in a gene set (assuming same design for all gene sets).
    #    - Gs: structure of gene dependencies between samples in a gene set (assuming same design for all gene sets).
    #    - Bs: structure of batch dependencies between samples in a gene set (assuming same design for all gene sets).
    #    - set: vector of indices indicating gene set for each sample.
    #Output:
    #    - maxLoglik: maximum likelihood
    
    N <- length(Y)
    p <- ncol(X)
    nset <- length(unique(set))
    n <- N/nset
    
    phi <- theta[3]
    ga <- theta[4]    
    #Time variance, gene variance, random error variance, and variances for other random factors.
    alpha <- c(1,theta[-c(3,4)])
    
    #This is identical for all gene sets (assuming identical design)
    #Time and gene correlations
    Rt <- exp(-phi*Ts)
    Rg <- exp(-ga*Gs)
    
    #Make covariance matrix
    if( ncol(Z) > 2 ) {
      #only works when we have a batch effect. Else use makecov()
      V.hat <- alpha[4]*Bs + alpha[1]*Rt + alpha[2]*Rg + diag(n)*alpha[3]
    }else V.hat <- makecov2(alpha,Rt=Rt,Rg=Rg,design=Z)  
    
    Vinv <- solve(V.hat)
    XVinv <- t(X)%*%Vinv
    XVinvX <- XVinv%*%X
    XV <- solve(XVinvX)%*%XVinv
    
    #This is different for each gene set (includes Y)
    RSS <- numeric(nset) 
    for(k in 1:nset) {
      
      idx <- which(set==k)
      Yk <- Y[idx]
      
      beta.hat <- XV%*%Yk
      r <- Yk - X%*%beta.hat
      RSS[k] <- t(r)%*%Vinv%*%r
    }
    
    #Add extra minus before the log-likelihood because R optimisation functions find minimum on default. 
    #Minimising negative likelihood is equivalent to maximising positive likelihood.
    sigmat.hat <- sum(RSS)/(nset*(n-p)) 
    maxLoglik <- - ( - 0.5*( nset*(n-p)*log(sigmat.hat) + nset*log( det(V.hat) )   + sum(RSS)/sigmat.hat + nset*log( det( XVinvX ) )  ) )
    
    maxLoglik     
    
  }


householder <-
function(x){ 
        m <- length(x)
        alpha <- sqrt(drop(crossprod(x)))
        e <- c(1,rep(0,m-1))
        u <- x - alpha*e
        v <- u/sqrt(drop(crossprod(u)))
        diag(m) - 2*v%*%t(v)
}
signProp <-
function( graph, V ){
  #
  #   Computes the propagation of correlation signs through a network.
  #
  #   graph is a data structure (see package igraph) that describes a network topology.
  #   V is a signed adjacency matrix, i.e. all nonzero elements indicate neighbors, and 
  #   the signs indicate the relation (positive or negative correlation).
  #
  #   Originally by Solve Sabo
  #   Biostatistics, Norwegian University of Life Sciences
  #
  #   Edited by
  #   Lars Snipen, 18. March 2009
  #   Guro Dorum, 19. July 2016
  p <- vcount( graph )
  #allpaths <- shortest.paths( graph, v=0:(p-1), mode = "all" )
  allpaths <- igraph::shortest.paths( graph, mode = "all" )
  allpaths[which( allpaths=="Inf" )] <- p
  max.length <- max( allpaths )
  
  R <- V
  Vprod <- V
  for( j in 2:max.length ){
    
    #Keeping track of signs
    Vprod <- Vprod %*% V
    signs <- sign( Vprod )
    
    #Finding paths from order 2 to max.length
    id <- which( allpaths==j )
    R[id] <- signs[id]    
  }
  diag( R ) <- 1
  return( R )
}

loglik_genedep <- function(theta, Xlist, Zlist, Ylist, TsList, GsList, BsList) {
  #log-likelihood function to be optimised with reml().
  #
  #Input:
  #    - theta: a vector of parameters to be optimised over. The first parameter must be gene variance, 
  #      the second parameter random error variance, the third must be a parameter controlling correlation between time points and the fourth must
  #      be a parameter controlling correlation between genes. The remaining parameters can be variances for other random design factors. 
  #    - X: design matrix with fixed design variables for a gene set (assuming same design for all gene sets).
  #    - Z: design matrix for random design factors for a gene set (assuming same design for all gene sets).
  #    - Y: vector of observed values for all gene sets.
  #    - Ts: structure of time dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - Gs: structure of gene dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - Bs: structure of batch dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - set: vector of indices indicating gene set for each sample.
  #Output:
  #    - minLoglik: minimum negative log-likelihood ( = maximum likelihood).
  
  p <- ncol(Xlist[[1]])
  nset <- length(Xlist)
  
  phi <- theta[3]
  ga <- theta[4]    
  #Time variance, gene variance, random error variance, and variances for other random factors.
  alpha <- c(1,theta[-c(3,4)])
  
  # #This is identical for all gene sets (assuming identical design)
  #    #Time and gene correlations
  #    Rt <- exp(-phi*Ts)
  #    Rg <- exp(-ga*Gs)
  #    
  #    #Make covariance matrix
  #    if( ncol(Z) > 2 ) {
  #        #only works when we have a batch effect. Else use makecov()
  #        V.hat <- alpha[4]*Bs + alpha[1]*Rt + alpha[2]*Rg + diag(n)*alpha[3]
  #    }else V.hat <- makecov_genedep(alpha,Rt=Rt,Rg=Rg,design=Z)  
  #     
  #    Vinv <- solve(V.hat)
  #    XVinv <- t(X)%*%Vinv
  #    XVinvX <- XVinv%*%X
  #    XV <- solve(XVinvX)%*%XVinv
  #    
  #This is different for each gene set (includes Y)
  RSS <- numeric(nset) 
  logDetV <- numeric(nset)
  logDetXVinvX <- numeric(nset)
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
    Rt <- exp(-phi*Ts)
    Rg <- exp(-ga*Gs)
    
    #Make covariance matrix
    if( ncol(Z) > 2 ) {
      #only works when we have a batch effect. Else use makecov()
      V.hat <- alpha[4]*Bs + alpha[1]*Rt + alpha[2]*Rg + diag(n)*alpha[3]
    }else V.hat <- alpha[1]*Rt + alpha[2]*Rg + diag(n)*alpha[3]#V.hat <- makecov_genedep(alpha,Rt=Rt,Rg=Rg,design=Z)  
    
    Vinv <- solve(V.hat)
    XVinv <- t(X)%*%Vinv
    XVinvX <- XVinv%*%X
    XV <- solve(XVinvX)%*%XVinv
    
    beta.hat <- XV%*%Y
    r <- Y - X%*%beta.hat
    
    RSS[k] <- t(r)%*%Vinv%*%r        
    logDetV[k] <- log( det(V.hat) )
    logDetXVinvX[k] <- log( det( XVinvX ) )
  }
  
  
  #Add extra minus before the log-likelihood because R optimisation functions find minimum on default. 
  #Minimising negative likelihood is equivalent to maximising positive likelihood.
  sigmat.hat <- sum(RSS)/(N-nset*p) 
  #minLoglik <- - ( - 0.5*( nset*(n-p)*log(sigmat.hat) + nset*log( det(V.hat) )   + sum(RSS)/sigmat.hat + nset*log( det( XVinvX ) )  ) )
  minLoglik <- - ( - 0.5*( (N-nset*p)*log(sigmat.hat) + sum(logDetV)  + sum(RSS)/sigmat.hat + sum(logDetXVinvX) ) )
  
  minLoglik     
  
}
