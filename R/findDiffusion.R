#' Find diffusion matrix
#' 
#' Find the diffusion matrix based on a graph.
#' @param graph Graph object
#' @param alpha Parameter that controls speed of diffusion over a network
#' @author Solve Sabo <solve.sabo@nmbu.no>
#' @importFrom igraph as.undirected get.adjacency degree
#' @export
findDiffusion <-
function(graph, alpha=1){
  #graph must be undirected to get symmetric matrix for eigenvalue decomposition
  graph <- as.undirected(graph)
  A <- get.adjacency(graph) #Returns a matrix indicating whether nodes are connected to each other (1) or not (0)
  degs <- degree(graph)
  L <- A - diag(degs)
  S <- eigen(L)
  K <- S$vectors%*%diag(exp(alpha*S$values))%*%t(S$vectors)
  idx <- which( K < 0 ) # Removing negative diffusion
  K[idx] <- 0
  K
}
