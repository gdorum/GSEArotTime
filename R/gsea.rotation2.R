#' GSEA rotation for longitudinal data
#' 
#' GSEA with rotation test for longitudinal data. Assumes common covariance matrix for all genes in the same gene set.
#' @param S Matrix indicating gene set membership with genes as rows and gene sets as columns. 1 indicates that the gene is a member, 0 indicates not a member.
#' @param Y Matrix of gene expression data. Genes are columns and samples are rows
#' @param Xlist Design matrix for fixed factors, each column represent a factor. A list of matrices, one for each gene set.
#' @param Zlist Design matrix for random factors, each column represent a factor. A list of matrices, one for each gene set.
#' @param contrast Matrix of contrasts to be tested. Each column represent one contrast.
#' @param covmatList Covariance matrix. A list of matrices, one for each gene set.
#' @param nrot Number of rotations
#' @param ES.p Weighting parameter for enrichment score (default=1)
#' @return  \item{Ngenes }{Size of gene sets}
#' \item{ES }{Enrichment score}
#' \item{NES }{Normalised enrichment score}
#' \item{p.value }{A p-value per gene set}
#' \item{q.value }{An FDR q-value for each gene set}
#' @references 
#' Dorum, G., Snipen, L., Solheim, M., & Sabo, S. (2014). 
#' Rotation gene set testing for longitudinal expression data. Biometrical Journal, 56(6), 1055-1075.
#' @author Guro Dorum
#' @seealso \code{\link{gsea.rotation}} for GSEA rotation that assumes common covariance matrix for all genes.
#' @importFrom limma squeezeVar
#' @export
#' @examples 
#' #2 gene sets of 3 genes each. 3 time points. 1 replicate at each time point. 
#' #Total number of observations in gene set is n
#' #Note that this example is too small to get good covariance estimates 
#' #(the model is too complex). Increase the number of samples to improve the estimates
#' ngs <- 5
#' ng <- 10
#' nt <- 3
#' nrep <- 3
#' n <- ng*nt*nrep
#'
#'iset <- rep(1:ngs, each=ng) #Indicating which gene set each gene belongs to
#'S <- matrix(0,ng*ngs,ngs)
#'for(i in 1:ngs) {
#'  S[iset==i,i] <- 1
#'}
#'
#'#Find correlations between genes in a gene set by construting a graph for each gene set 
#'#and finding the diffusion matrix. For simplicity, all gene sets are assumed to have 
#'#the same gene dependencies here.
#'G.data <- cbind(1:(ng-1),2:ng)
#'G.graph <- graph.data.frame( G.data, directed=FALSE )
#'G.graph <- addSignMatrix( G.graph )
#'Dmat <- findDiffusion( G.graph )
#'Dmat <- 1/Dmat # Gene dependencies. 'Distance' is reciprocal to diffusion
#'
#'#Design for each gene set
#'Gene <- rep(1:ng,each=nt*nrep)
#'Timer <- rep(1:nt,times=ng*nrep) 
#'Time <- ordered(Timer)
#'#Design matrix for random factors
#' Z <- cbind(Timer,Gene)
#' colnames(Z) <- c("Time","Gene")
#' #Design matrix for fixed factors
#' X <- model.matrix(~Time)[,1:2]
#' colnames(X) <- c("Gene","Linear")
#' 
#' set <- rep(1:ngs,each=n) #Indicating which gene set each sample belong to
#' 
#' #Set values for variance components
#' sigmat <- 2 #Time variance
#' sigmag <- 2 #Gene variance
#' sigmae <- 3 #Random error variance
#' phi <- 0.9 #Time correlation parameter
#' ga <- 0.3 #Gene correlation parameter
#' 
#' #Time and gene correlations
#' Ts <- matrix(0,n,n)
#' Gs <- matrix(0,n,n)
#' for(i in 1:(n-1)) {
#'   for(j in (i+1):n) {
#'     if(Z[i,2]==Z[j,2]){ 
#'       Ts[i,j] <- abs(Z[j,1] - Z[i,1]) 
#'       Gs[i,j] <- Dmat[Z[j,2],Z[i,2]] 
#'     } else Ts[i,j] <- Inf
#'   }
#' }
#' Ts <- Ts + t(Ts)
#' Rt <- exp(-phi*Ts)
#' Gs <- Gs + t(Gs)
#' Rg <- exp(-ga*Gs)
#' 
#' #Make covariance matrix
#' V <- makecov2(theta=c(sigmat,sigmag,sigmae),Rt=Rt,Rg=Rg,design=Z)
#' Vroot <- chol(V)
#' 
#' xz <- vector("list",ngs)
#' Xlist <- lapply(xz,function( x ){X})
#' Zlist <- lapply(xz,function( x ){Z})
#' 
#' #Simulate data with covariance structure
#' Y <- numeric()
#' for(k in 1:ngs) {
#'   Y1 <- matrix(rnorm(n),1,n)
#'   Y <- c(Y,crossprod(t(Y1),Vroot))
#' }
#' # Add fixed gene effects and linear time effecs first gene set
#' betaG <-  2; betaLT <- 3
#' Yeff <- matrix(Y,ncol=ngs)
#' #Gene effects + linear time effect in first gene set
#' Yeff[,1] <- Yeff[,1] + X[,1]*betaG + X[,2]*betaLT
#' 
#' #Estimate covariance matrix
#' #Constraints for optimisation and start values
#' A <- rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1))
#' b <- c(1e-6, 1e-6, 1e-6, 1e-6)
#' alphaStart <- c(sigmag/sigmat,sigmae/sigmat,phi,ga)
#' reml_res <- reml2(alphaStart=alphaStart, X=X, Y=c(Yeff), Z=Z, Ts=Ts, Gs=Gs, set=set, A=A, b=b)
#' sigmat.hat <- reml_res$sigma
#' sigmag.hat <- reml_res$alpha[1]*sigmat.hat; sigmae.hat <- reml_res$alpha[2]*sigmat.hat
#' phi.hat <- reml_res$alpha[3]; ga.hat <- reml_res$alpha[4]
#' #Estimated covariance matrix
#' V.hat <- makecov2(theta=c(sigmat.hat,sigmag.hat,sigmae.hat),
#' Rt=exp(-phi.hat*Ts),Rg=exp(-ga.hat*Gs),design=Z)
#' 
#' if( all(eigen(V.hat)$values>=0) ) { 
#'   a <- vector("list",ngs)
#'   covmatList <- lapply(a,function( x ){V.hat})
#' } else stop("Covariance matrix is not positive definite! Simulate new data")
#' 
#' #Contrasts for testing gene effect and linear time effect
#' contrast <- matrix(0,ncol(X),2)
#' colnames(contrast) <- c("Gene","Linear")
#' contrast[colnames(X)=="Gene",1] <- 1
#' contrast[colnames(X)=="Linear",2] <- 1
#' #GSEA rotation
#' gsea.rotation2(S=S, Y=c(Yeff), Xlist=Xlist, Zlist=Zlist, contrast=contrast, 
#' covmatList=covmatList, nrot=1000, ES.p=1)
gsea.rotation2 <-
function(S, Y, Xlist, Zlist, contrast, covmatList=NULL, nrot=10000, ES.p=1) {
  
  k <- ncol(contrast) #Number of interesting contrasts
  idx <- which( contrast==1, arr.ind=TRUE)[,1] #Get the index of interesting contrasts
  g <- nrow(S) #Number of genes
  n <- length(Y)/g #Number of samples per gene (assuming same number of samples)
  p <- nrow(contrast) #Number of design factors
  d <- n - p #Number of independent samples to calculate variance from
  p0 <- p-k #Number of uninteresting effects
  nset <- ncol(S)
  
  m <- 0
  l <- 0
  eff <- matrix(0,n,g)
  for(i in 1:nset) {     
    
    Xi <- Xlist[[i]]  
    Zi <- Zlist[[i]]
    yi <- Y[(m+1):(m+nrow(Xi))]
    dim(yi) <- c(nrow(Xi),1)
    m <- m + nrow(Xi)
    
    #Transform data to remove correlation between samples
    if( !is.null( covmatList ) ) {
      R <- chol(covmatList[[i]]) 
      Xi <- backsolve(R, Xi, transpose=TRUE)  
      yi <- backsolve(R, yi, transpose=TRUE) 
    }
    
    Xi <- cbind(Xi[,-idx],Xi[,idx]) #Arrange design matrix so interesting columns are last
    
    gi <- length(unique(Zi[,2])) #Number of genes in gene set
    
    for(j in 1:gi) {
      
      l <- l + 1
      geneid <- which(Zi[,2]==j) 
      Xj <- Xi[geneid,]
      yj <- yi[geneid]    
      
      QRdata <- qr(Xj)
      Q <- qr.Q(QRdata,complete=TRUE)  #Adding n-p random orthogonal columns
      eff[,l] <- t(Q)%*%yj  #Project y onto orthonormal basis for X. Projected onto the last n-p random orthogonal columns we achieve independent residuals.
    }
  } #nset
  
  if(p0>0) {
    eff <- eff[-(1:p0),] #Remove first p-k uninteresting rows
  }          
  res <- eff[(k+1):nrow(eff), , drop=FALSE] #n-p independent residuals to calculate variance from
  
  s2 <- colMeans(res^2) #Estimate variance with independent residuals (last n-p columns)
  sv <- squeezeVar(s2, df = d)
  d0 <- sv$df.prior
  s02 <- sv$var.prior
  sd.post <- sqrt(sv$var.post)
  
  B <- eff[1:k,] 
  dim(B) <- c(k,g) 
  sd.post.mat <- matrix(sd.post,nrow(B),ncol(B),byrow=TRUE) 
  modt <- B/sd.post.mat #Modererte t-verdier.
  #Compute F-values if testing more than one contrast, else keep t-values
  if( k > 1 ){
    modf <- apply(modt,2,FUN=function(x) sum(x^2)/length(x)) #En moderert F-verdi for hvert gen (sum av kvadrerte t-verdier delt p? antall kontraster)
  } else modf <- modt
  
  escore <- numeric(nset)
  for (j in 1:nset) {
    escore[j] <- es(modf,S[,j],p=ES.p)
  }
  
  #Rotation test
  modfr <- matrix(0,nrot,g)
  esrot <- matrix(0,nrot,nset)
  for(i in 1:nrot) {
    
    #Rotating y
    Z <- matrix(rnorm((d+k)^2),(d+k),(d+k))
    QRdata <- QR(Z)
    Qr <- QRdata$Q
    effr <- Qr %*% eff
    
    resr <- effr[(k+1):nrow(effr), , drop=FALSE] #n-p independent residuals to calculate variance from
    
    s2r <- colMeans(resr^2) #Estimate variance with independent residuals (last n-p columns)
    if (is.finite(d0)){ 
      sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
    } else{
      sdr.post <- sqrt(s02)
    }
    
    Br <- effr[1:k,] 
    dim(Br) <- c(k,g) 
    sdr.post.mat <- matrix(sdr.post,nrow(Br),ncol(Br),byrow=TRUE)
    modtr <- Br/sdr.post.mat#Modererte t-verdier.
    if( k > 1 ) {
      modfr[i,] <- apply(modtr,2,FUN=function(x) sum(x^2)/length(x)) #En moderert F-verdi for hvert gen (sum av kvadrerte t-verdier delt p? antall kontraster)d
    } else modfr[i,] <- modtr
    
    for (j in 1:nset) {
      esrot[i,j] <- es(modfr[i,],S[,j],p=ES.p)
    }
    
  }  
  
  #genewise p-values
  pvals <- numeric()
  for(j in 1:g) {
    pvals[j] <- sum(modfr[,j]>=modf[j])/nrot
  } 
  
  #Estimate significance
  p.value <- numeric(nset)
  NES <- numeric(nset)
  NES.null <- matrix(0,nrot,nset)
  for( j in 1:nset) {
    sig <- significance(escore[j], esrot[,j])
    p.value[j] <- sig$p.value            #p-value
    NES[j] <- sig$NES                    #Normalised enrichment score
    NES.null[,j] <- sig$NESnull          #Normalised null distribution
  }
  
  #q-values
  #If the test statistic is F, use a "one-tail test", else use a "two-tail test"
  q.value <- numeric(nset)
  if(k > 1) {
    for(j in 1:nset) {
      q.value[j] <- ( ( sum(NES.null[,j] >= NES[j]) + 1)/( nrot + 1 ) ) / ( sum(NES >= NES[j])/nset ) 
    }
    q.value <- ifelse( q.value > 1, 1, q.value)
  } else q.value <- sapply(NES, FUN=fdr, NESobs=NES, NESnull=NES.null)
  
  
  len.set <- colSums(S)
  return(list(NGenes = len.set, ES=escore, NES=NES, p.value=p.value, q.value=q.value))
}
