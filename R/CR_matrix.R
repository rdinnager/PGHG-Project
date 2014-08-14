#' Function to generate consumer-resource consumption probability matrix with phylogenetic structure
#' @param ctree A phylogenetic tree of class phylo, representing the phylogenetic relationships
#' for the consumers
#' @param rtree A phylogenetic tree of class phylo, representing the phylogenetic relationships
#' for the resources
#' @param sig Total Variance to be generated in the latent variables via the phylogenies
#' @param phi Unit simplex vector (elements between 0 and 1 and sum to 1) describing how the
#' phylogenetic variance is partitioned into different kinds of phylogenetic variability.
#' The elements correspond to: 
#' \itemize{
#'  \item{1) Phylogenetic structure in resource richness}
#'  \item{2) Phylogenetic structure in consumer range}
#'  \item{3) Phylogenetic structure in resources' consumer sets}
#'  \item{4) Phylogenetic structure in consumers' resource sets}
#'  \item{5) Phylogenetic interactive structure between consumer and resource phylogenies}
#' }  
#' @param alpha Intercept for latenet variables (higher values mean probabilities are uniformly
#' larger).
#' @return A consumer-resource consumption probability matrix with consumers in the rows, and 
#' resources in the columns
#' @import MASS boot
#' @export
gen_CR_matrix <- function(ctree = rcoal(15), rtree = rcoal(25), sig = 1, phi = c(0.2, 0.2, 0,2, 0.2, 0.2), alpha = 0) {
  ## generate correlation matrices
  Ac <- vcv(ctree,corr=T)
  Ar <- vcv(rtree,corr=T)
  ## store dimensions
  cdim <- dim(Ac)[1]
  rdim <- dim(Ar)[1]
  ## make J matrice
  Jc <- matrix(rep(1, cdim^2), nrow = cdim, ncol = cdim)
  Jr <- matrix(rep(1, rdim^2), nrow = rdim, ncol = rdim)
  ## make I matrices
  Ic <- diag(cdim)
  Ir <- diag(rdim)
  ## make equation elements
  ele1 <- kronecker(Ar, Jc) ## resources nested in consumers (first go through all hosts for parasite 1, then all hosts for parasite 2, etc.)
  ele2 <- kronecker(Jr, Ac)
  ele3 <- kronecker(Ar, Ic)
  ele4 <- kronecker(Ir, Ac)
  ele5 <- kronecker(Ar, Ac)## set parameters
  ## calculate equation
  fullcor<-sig*(phi[1]*ele1+phi[2]*ele2+phi[3]*ele3+phi[4]*ele4+phi[5]*ele5)
  ## generate latent gaussian variables
  lat_vars <- mvrnorm(1, mu = rep(alpha, cdim*rdim), Sigma = fullcor)
  ## covert to matrix of probabilities
  probs <- inv.logit(lat_vars)
  prob_mat <- matrix(probs, nrow = cdim, ncol = rdim)
  rownames(prob_mat) <- rownames(Ac)
  colnames(prob_mat) <- colnames(Ar)
  return(prob_mat)
}

#' Function to plot a consumer-resource consumption probability matrix, along with 
#' phylogenies
#' @param ctree A phylogenetic tree of class phylo, representing the phylogenetic relationships
#' for the consumers
#' @param rtree A phylogenetic tree of class phylo, representing the phylogenetic relationships
#' for the resources
#' @param prob_mat Consumer-resource consumption probability matrix like that produced 
#' by \code{\link{gen_CR_matrix}}.
#' @export
plot_CR_matrix <- function(ctree, rtree, prob_mat) {
## plot probs
  heatmap(prob_mat,  Rowv = as.dendrogram(as.hclust(ctree)), Colv = as.dendrogram(as.hclust(rtree)), labRow = ctree$tip.label, labCol = rtree$tip.label, scale="none", margins = c(5, 6), font.axis = 3, col = rev(grey(seq(0, 1, length = max(prob_mat)))), cexRow = 1, cexCol = 1) 
}