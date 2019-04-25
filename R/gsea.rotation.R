#' GSEA rotation for longitudinal data
#' 
#' GSEA with rotation test for longitudinal data. Assumes common covariance matrix for all genes.
#' @param S Matrix indicating gene set membership with genes as rows and gene sets as columns. 1 indicates that the gene is a member, 0 indicates not a member.
#' @param y Matrix of gene expression data. Genes are columns and samples are rows
#' @param X Design matrix. Each row represent a sample
#' @param contrast Matrix of contrasts. Each column represent one contrast
#' @param covmat Covariance matrix (optional)
#' @param nrot Number of rotations
#' @param ES.p Weighting parameter for enrichment score (default=1)
#' @return  \item{Ngenes }{Size of gene sets}
#' \item{ES }{Enrichment score}
#' \item{NES }{Normalised enrichment score}
#' \item{p.value }{A p-value per gene set}
#' \item{q.value }{An FDR q-value for each gene set}
#' @references 
#' Dorum, G., Snipen, L., Solheim, M. and Sabo (2009) Rotation Testing in Gene Set Enrichment
#' Analysis for Small Direct Comparison Experiments. \emph{Statistical Applications in Genetics
#'  and Molecular Biology}, \bold{8}(1), article 34.
#'  
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#' Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
#' and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
#' approach for interpreting genome-wide expression profiles, \emph{PNAS}, \bold{102},
#' 15545-15550.
#' @author Guro Dorum
#' @seealso \code{\link{gsea.rotation2}} for GSEA rotation that assumes common covariance matrix for all genes in the same gene set.
#' @importFrom limma squeezeVar
#' @export
#' @examples
#' #2 gene sets of 3 genes each. 3 time points. 1 replicate at each time point. 
#' #Total number of samples per gene is n
#' ngs <- 5
#' ng <- 3
#' nt <- 3
#' nrep <- 2
#' n <- nt*nrep
#' 
#' #Indicating which gene set each gene belongs to
#' iset <- rep(1:ngs, each=ng)
#' S <- matrix(0,ng*ngs,ngs)
#' for(i in 1:ngs) {
#'   S[iset==i,i] <- 1
#' }
#' 
#' #Design matrix per gene
#' Timer <- rep(1:nt,each=nrep) 
#' Time <- ordered(Timer)
#' #Design for random factors
#' Z <- cbind(Timer)
#' colnames(Z) <- c("Time")
#' #Design for fixed factors
#' X <- model.matrix(~Time)
#' X <- X[,-ncol(X)]
#' colnames(X) <- c("Gene","Linear")
#' 
#' #Set values for variance components
#' sigmat <- 2 #Time variance
#' sigmae <- 3 #Random error variance
#' phi <- 0.9 #Time correlation parameter
#' #Correlation between time points for same gene
#' Ts <- matrix(0,n,n)
#' for(i in 1:(n-1)) {
#'   for(j in (i+1):n) {
#'     Ts[i,j] <- abs(Z[j,1] - Z[i,1]) 
#'   }
#' }
#' Ts <- Ts + t(Ts)
#' R <- exp(-phi*Ts)
#' #Covariance matrix
#' subject <- rep(1,each=n)
#' V <- makecov(theta=c(sigmat,sigmae),correlation=R,design=Z,subject=subject)
#' Vroot <- chol(V)
#' 
#' #Simulate data with covariance structure
#' Y1 <- matrix(rnorm(n*ng*ngs),ng*ngs,n)
#' Y <- t(crossprod(t(Y1),Vroot)) 
#' #Add fixed gene effects and linear time effecs in first gene set
#' betaG <-  2; betaLT <- 1
#' Yeff <- Y
#' Yeff[,1:ng] <- Yeff[,1:ng] + matrix(X[,1]*betaG,ncol=ng,nrow=n) 
#' + matrix(X[,2]*betaLT,ncol=ng,nrow=n)
#' 
#' #Estimate covariance matrix
#' #Constraints for optimisation and start values
#' A <- rbind(c(1,0), c(0,1))
#' b <- c(1e-6, 1e-6)
#' alphaStart <- c(sigmae/sigmat,phi)
#' reml_res <- reml( alphaStart=alphaStart, X=X, Y=Yeff, Z=Z, subject=subject, A=A, b=b)
#' sigmat.hat <- reml_res$sigma
#' sigmae.hat <- reml_res$alpha[1]*sigmat.hat
#' phi.hat <- reml_res$alpha[2]
#' R.hat <- exp(-phi.hat*Ts)
#' #Estimated covariance matrix
#' V.hat <- makecov(theta=c(sigmat.hat,sigmae.hat),correlation=R.hat,design=Z,subject=subject)
#' 
#' #Contrasts for testing gene effect and linear time effect
#' contrast <- matrix(0,ncol(X),2)
#' colnames(contrast) <- c("Gene","Linear")
#' contrast[colnames(X)=="Gene",1] <- 1
#' contrast[colnames(X)=="Linear",2] <- 1
#' 
#' #Contrast for testing gene effect
#' contrast <- matrix(0,ncol(X),2)
#' colnames(contrast) <- c("Gene","Linear")
#' contrast[colnames(X)=="Gene",1] <- 1
#' contrast[colnames(X)=="Linear",2] <- 1
#' 
#' #GSEA rotation
#' gsea.rotation(S=S, y=t(Yeff), X=X, contrast=contrast, covmat=V.hat, nrot=1000, ES.p=1)
gsea.rotation <-
function(S, y, X, contrast, covmat=NULL, nrot=10000, ES.p=1) {

  k <- ncol(contrast) #Number of interesting contrasts
  idx <- which( contrast==1, arr.ind=TRUE)[,1] #Get the index of interesting contrasts
  n <- nrow(X)
  p <- ncol(X)
  d <- n - p #Number of independent samples to calculate variance from
  p0 <- p-k #Number of uninteresting effects
  nset <- ncol(S)
  
  #Transform data to remove correlation between samples
  if( !is.null( covmat ) ) {
    R <- chol(covmat)              
    y <- t(backsolve(R, t(y), transpose = TRUE))
    X <- backsolve(R, X, transpose = TRUE)     
  }
  
  X <- cbind(X[,-idx],X[,idx]) #Arrange design matrix so interesting columns are last
  
  QRdata <- qr(X)
  Q <- qr.Q(QRdata,complete=TRUE)  #Adding n-p random orthogonal columns
  #eff <- t(t(Q)%*%t(y))           #Project y onto orthonormal basis for X. Projected onto the last n-p random orthogonal columns we achieve independent residuals.
  eff <- t(Q)%*%t(y) 
  
  if(p0>0) {
    eff <- eff[-(1:p0),] #Remove first p-k uninteresting rows
  }          
  res <- eff[(k+1):nrow(eff), , drop=FALSE] #n-p independent residuals to calculate variance from
  #res <- eff[,(k+1):ncol(eff), drop=FALSE]
  
  s2 <- colMeans(res^2) #Estimate variance with independent residuals (last n-p columns)
  sv <- squeezeVar(s2, df = d)
  d0 <- sv$df.prior
  s02 <- sv$var.prior
  sd.post <- sqrt(sv$var.post)
  
  B <- eff[1:k,] 
  #B <- eff[,1:k] 
  dim(B) <- c(k,nrow(y)) 
  sd.post.mat <- matrix(sd.post,nrow(B),ncol(B),byrow=TRUE) 
  modt <- B/sd.post.mat #Modified t-values
  #Compute F-values if testing more than one contrast, else keep t-values
  if( k > 1 ){
    modf <- apply(modt,2,FUN=function(x) sum(x^2)/length(x)) #one modified F-value per gene (sum of squared t-values fivided by number of contrasts)
  } else modf <- modt
  
  
  escore <- numeric(nset)
  for (j in 1:nset) {
    escore[j] <- es(modf,S[,j],p=ES.p)
  }
  
  #Rotation test
  modfr <- matrix(0,nrot,nrow(y))
  esrot <- matrix(0,nrot,nset)
  
  for(i in 1:nrot) {
    
    #Rotating y
    Z <- matrix(rnorm((d+k)^2),(d+k),(d+k))
    QRdata <- QR(Z)
    Qr <- QRdata$Q
    effr <- Qr %*% eff
    
    #effr <- rotation(eff,method=2)
    
    resr <- effr[(k+1):nrow(effr), , drop=FALSE] #n-p independent residuals to calculate variance from
    
    s2r <- colMeans(resr^2) #Estimate variance with independent residuals (last n-p columns)
    if (is.finite(d0)){ 
      sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
    } else{
      sdr.post <- sqrt(s02)
    }
    
    Br <- effr[1:k,] 
    dim(Br) <- c(k,nrow(y)) 
    sdr.post.mat <- matrix(sdr.post,nrow(Br),ncol(Br),byrow=TRUE)
    modtr <- Br/sdr.post.mat
    if( k > 1 ) {
      modfr[i,] <- apply(modtr,2,FUN=function(x) sum(x^2)/length(x))
    } else modfr[i,] <- modtr
    
    
    for (j in 1:nset) {
      esrot[i,j] <- es(modfr[i,],S[,j],p=ES.p)
    }
    
  }  
  
  #genewise p-values
  pvals <- numeric()
  for(j in 1:nrow(y)) {
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
