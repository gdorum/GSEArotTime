#' Make covariance matrix
#' 
#' Construct a covariance matrix for a model with time series effect and other random effetcs. For balanced and unbalanced designs.
#' @param theta A vector of variances for all random factors, where first element must be time series variance and second 
#' element must be the random error variance. The remaining elements can be variances for any other random design factors
#' @param correlation A matrix of time correlations between samples (usually zero for samples from different subjects)
#' @param design A nxp design matrix for the p random design factors. Must contain time effect in the first column. 
#' The time column must be numeric and the time intervals must be of actual size
#' @param subject A numeric vector indicating which subject each sample belongs to
#' @return Covariance matrix
#' @author Guro Dorum
#' @seealso \code{\link{makecov2}} for construction of a covariance matrix that also includes gene dependencies.
#' @export
makecov <-
function(theta, correlation, design, subject) {

  N <- nrow(design)
  
  V <- matrix(0,N,N)
  #Loop through each element of upper triangle
  for( i in 1:(N-1) ) {
    for( j in (i+1):N ) {
      #Time effect if same subject
      if( subject[i] == subject[j] ) {
        V[i,j] <- theta[1]*correlation[i,j]
      }
      
      #Remaining factors in design matrix
      if( ncol(design) > 1 ) {
        for( k in 2:ncol(design) ) {
          if( design[i,k] == design[j,k] ) {
            V[i,j] <- V[i,j] + theta[k+1]
          }
        }
      }
    }
  }
  V <- V + t(V)
  diag(V) <- sum(theta)
  
  V        
}
