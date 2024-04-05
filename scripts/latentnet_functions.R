## R file to host latentnet functions ------------------------------------------
library(latentnet)


#' Latentnet Network Generator
#'
#' Function to generate network according to Handcock's RSS latentnet paper 
#' 
#' @param n_nodes A positive integer representing the total number of nodes in 
#' the network
#' @param K A positive integer indicating the true number of clusters
#' @param x_dim A positive integer indicating the number of node-specific
#' covariates to generate.
#' @param betavec Vector of intercept and slope terms
#' @param noise noise/signal ratio (default=1)
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' 
#' @return Network as list object with characteristics
#' 
#' @references Handcock, M. S., A. E. Raftery, and T. J. M. (2007). Model-based 
#' clustering for social networks. Journal of the Royal Statistical Society. 
#' Series A 170 (2), 307–354
#' 
#' @examples
#' links <- generate_latentnet_network(200, 2, 2, c(1, 1))
#' print(find_sbm(links$A, links$z))
#' 
#' @export
#' 
generate_latentnet_network <- function(n_nodes, K, x_dim, betavec, noise = 1, 
                                       directed = FALSE){
  z_tru <- sample(1:K, n_nodes, replace = TRUE)
  centroids <- matrix(stats::rnorm(x_dim * K), K, x_dim)
  X <- centroids[z_tru, ] + stats::rnorm(x_dim * n_nodes, sd = sqrt(noise))
  ## Generate probabilities, and link probability matrix
  ## Using Euclidean distance, calculate true link probabilities
  ZZ <- stats::dist(X, upper = TRUE, diag = TRUE)
  S_ij <- matrix(0, n_nodes, n_nodes)
  S_ij[upper.tri(S_ij)] <- ZZ; S_ij <- as.matrix(Matrix::forceSymmetric(S_ij))
  diag(S_ij) <- 0
  eta <- betavec[1] + betavec[2] * S_ij; diag(eta) <- -Inf
  link_prob_tru <- exp(eta) / (1 + exp(eta))
  ## Finally, generate adjacency matrix and make symmetric if undirected
  A <- matrix(stats::rbinom(n_nodes^2, 1, link_prob_tru), n_nodes, n_nodes)
  if (!directed){ A <- as.matrix(Matrix::forceSymmetric(A)) }
  return (list(A = A, z = z_tru, X = X, lp = link_prob_tru, dis = S_ij))
}


#' Latentnet Cluster Finder
#' 
#' Function to find the optimal number of groups using latentnet package
#' 
#' @param p Number of clusters. The function will iterate from 1 to p
#' @param AA Network object (from network package) containing adjacency and 
#' covariate information
#' @param sample_size Number of non-rejected MCMC samples
#' @param burnin Number of initial MCMC samples to be discarded
#' @param by Number of samples to skip after burn-in period
#' 
#' @return A dataframe with the AIC and BIC values for each K value
#'
#' @references Handcock, M. S., A. E. Raftery, and T. J. M. (2007). Model-based 
#' clustering for social networks. Journal of the Royal Statistical Society. 
#' Series A 170 (2), 307–354
#' 
#' @examples 
#' \dontrun{
#' set.seed(123)
#' net <- network::network(m <- matrix(stats::rbinom(2500, 1, 0.4), 50, 50))
#' find_K_optimal(4, net, 1000, 200, 1)
#' }
#' 
#' @export
#' 
find_K_optimal <- function(p, AA, sample_size, burnin, by){
  ## Set a dataframe for BIC and number of clusters
  best_estimate <- data.frame(K = 1:p, BIC = rep(Inf, p), AIC = rep(Inf, p))
  n <- AA$gal$n
  z_est <- rep(1, p)
  for (i in 1:p){
    samp.fit <- latentnet::ergmm(AA ~ euclidean(d = 2, G = i), 
                                 control = latentnet::control.ergmm(
                                   sample.size = sample_size, burnin = burnin, interval = by), 
                                 verbose = TRUE)
    fit <- summary(samp.fit)$bic$overall
    if (fit < min(best_estimate$BIC)){
      z_est <- samp.fit$mkl$Z.K
    }
    best_estimate$BIC[i] <- fit
    n_params <- (i^2 + i) / 2
    best_estimate$AIC[i] <- fit + 2 * n_params - 
      n_params * log(n * (n - 1) / 2)
    print(c(i, fit))
  }
  return (list(K = which.min(best_estimate$BIC), 
               data = best_estimate, 
               z = z_est))
}








