##
## calfsbm: Covariate-assisted Latent Factor Stochastic block model
## Copyright (C) 2024  Sydney Louit, Jun Yan, and Panpan Zhang
## Jun Yan <jun.yan@uconn.edu>
##
## This file is part of the R package calfsbm.
##
## The R package calfsbm is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package calfsbm is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

#' @importFrom arm bayesglm sim
#' @importFrom boa boa.chain.gandr
#' @importFrom cluster pam
#' @importFrom label.switching aic
#' @importFrom Matrix forceSymmetric
#' @importFrom network network
#' @import nimble
NULL


#' Visualize Link Probability vs. Actual Matrix
#' 
# Using the adjacency matrix, returns a plot of the adjacency matrix, and 
# a raster of the true probabilities
# @param observed Observed adjacency matrix
# @param prob Matrix of true probabilities
# @param z_tru Vector of the true node membership
# @return Graph of the connectivity and graph of true probabilities
# @examples links <- sim_calfsbm(n_nodes = 50, K = 2, m = 2, prob = c(0.5, 0.5),
#beta0 = 1, beta = diag(2) - 3, sigma = 0.3, spat = 0.5)
# plot_conn(links$A, links$lp, links$z)
# @export
#plot_conn <- function(observed, prob, z_tru){
#    n <- length(z_tru)
#    par(mar = c(1, 1, 1, 1))
#    par(mfrow = c(1, 2))
#    B <- observed[order(z_tru), order(z_tru)]
#    xx <- raster::rasterFromXYZ(cbind(expand.grid(1:n, n:1), c(B)))
#    raster::plot(xx, legend=FALSE, axes=FALSE)
#    ## Plot true probabilities
#    B <- prob[order(z_tru), order(z_tru)]
#    xx <- raster::rasterFromXYZ(cbind(expand.grid(1:n, n:1), c(B)))
#    raster::plot(xx, axes=FALSE)
#}


#' Helper function for CALF-SBM Initialization
#' 
#' Get the X matrix for logistic regression by taking clusters, links, and dist.
#' @param initial_z Initial node membership vector
#' @param A Known adjacency matrix
#' @param S_ij Distance matrix
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param offset Logistic regression offset terms
#' 
#' @return List of X matrix and grid of (x, y) indices
#' 
#' @export
gen_factor <- function(initial_z, A, S_ij, directed = FALSE, offset = NULL){
    ## Initialize number of nodes/clusters
    n <- length(initial_z)
    K <- length(unique(initial_z))
    degrees <- rowSums(A)
    eg <- expand.grid(1:n, 1:n)
    cluster <- cbind(initial_z[eg[, 1]], initial_z[eg[, 2]])
    if (!directed){
        wuta <- which(upper.tri(A))
        cluster <- cluster[wuta, ]
        ## Get beta coordinates
        cluster <- t(apply(cluster , 1, sort, decreasing = TRUE))
        cluster <- (cluster[, 1] - 1) * K + cluster[, 2]
        cluster <- data.frame(y = A[wuta], cl = cluster, x = S_ij[wuta])
    } else {
        wuta <- which(upper.tri(A) | lower.tri(A))
        cluster <- cluster[wuta, ]
        ## Get beta coordinates
        cluster <- (cluster[, 1] - 1) * K + cluster[, 2]
        cluster <- data.frame(y = A[wuta], cl = cluster, x = S_ij[wuta])
    }
    if (length(offset) > 0){
        ## Offset terms theta_i and theta_j to be used in logistic model
        edge_offset <- offset[eg[wuta, 1]] + offset[eg[wuta, 2]]
        cluster$offset <- edge_offset
    }
    return (list(cluster = cluster, grid = eg[wuta,]))
}


#' CALF-SBM Network Generation Function
#'
#' Generates network using probabilities, betas, and signal strength
#' @param n_nodes A positive integer representing the total number of nodes in 
#' the network to be generated
#' @param K A positive integer indicating the true number of clusters in the 
#' network to be generated.
#' @param n_covar A positive integer indicating the number of node-specific
#' covariates to generate. As part of the list object returned by the function, 
#' there will be an n_nodes by n_covar matrix representing the covariates.
#' @param prob Vector of weights of being in each group. This should be
#' a numeric vector of length \code{K}, and can be either a vector of weights 
#' or a vector of probabilities
#' @param beta0 Numeric value for the CALF-SBM's intercept term. Lower values 
#' correspond to higher sparsity
#' @param beta \code{K}-by-\code{K} numeric matrix of within and between-cluster 
#' coefficients. Diagonal values correspond to within-cluster effects. If the 
#' network is undirected, this matrix should be symmetric.
#' @param sigma Non-negative number for the standard deviation of the random 
#' effects offset term, set to 0 if not using offset
#' @param spat Non-negative number to represent the ratio of between-cluster 
#' variance to within-cluster variance of the covariates. 0 indicates no signal
#' @param directed Logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @return List of adjacency, membership, link probability, and distance
#' 
#' @examples 
#' n <- 100
#' K_tru <- 2
#' m <- 2
#' sizes <- rep(1 / K_tru, K_tru)
#' beta_0 <- 1
#' beta_mat <- diag(K_tru) - 3
#' hetero <- 0.3
#' sn_ratio <- 1.5
#' links <- sim_calfsbm(n, K_tru, m, sizes, beta_0, beta_mat, hetero, sn_ratio)
#' ## Using adjacency and true clustering, get the density of each block
#' print(find_sbm(links$A, links$z))
#' 
#' @export
sim_calfsbm <- function(n_nodes, K, n_covar, prob, beta0, beta, 
                                sigma, spat, n_dummy = 0, directed = FALSE){
    z_tru <- sample(1:K, n_nodes, replace = TRUE, prob = prob)
    ## Spatial correlation
    X <- matrix(0, nrow = n_nodes, ncol = n_covar)
    ## Centers of new clusters
    angles <- 2 * pi * (runif(1) + 1:K) / K
    initial_mids <- cbind(sin(angles), cos(angles), -sin(angles), -cos(angles))
    ## Initialize midpoints to have variance 1
    mids <- sqrt(spat * 2) * initial_mids[, 1 + 1:n_covar %% K]
    for (i in 1:K){
      cli <- which(z_tru == i)
#      X[cli, ] <- sweep(MASS::mvrnorm(length(cli), mu = rep(0, n_covar), 
#                Sigma = diag(n_covar)), 2, mids[i, ], '+')
      X[cli, ] <- sweep(matrix(rnorm(length(cli) * n_covar), 
                               length(cli)), 2, mids[i, ], '+')
    }
    ## Degree
    theta <- stats::rnorm(n_nodes, 0, sigma)
    ## Generate probabilities, and link probability matrix
    ## Using Euclidean distance, calculate true link probabilities
    ZZ <- stats::dist(X, upper = TRUE, diag = TRUE)
    S_ij <- matrix(0, n_nodes, n_nodes)
    S_ij[upper.tri(S_ij)] <- ZZ; S_ij <- as.matrix(Matrix::forceSymmetric(S_ij))
    diag(S_ij) <- 0
    eta <- beta0 + outer(theta, theta, '+') + beta[z_tru, z_tru] * S_ij
    diag(eta) <- -Inf
    link_prob_tru <- 1 / (1 + exp(-eta))
    if (n_dummy != 0){
      X <- cbind(X, matrix(rnorm(n_nodes * n_dummy, sd = sqrt(spat + 1)), 
                           n_nodes, n_dummy))
      ZZ <- stats::dist(X, upper = TRUE, diag = TRUE)
      S_ij <- matrix(0, n_nodes, n_nodes)
      S_ij[upper.tri(S_ij)] <- ZZ; S_ij <- as.matrix(Matrix::forceSymmetric(S_ij))
      diag(S_ij) <- 0
    }
    ## Finally, generate adjacency matrix and make symmetric if undirected
    A <- matrix(stats::rbinom(n_nodes^2, 1, link_prob_tru), n_nodes, n_nodes)
    if (!directed){ A <- as.matrix(Matrix::forceSymmetric(A)) }
    return (list(A = A, z = z_tru, X = X, lp = link_prob_tru, 
                 dis = S_ij, theta = theta, beta0 = beta0, beta = beta))
}


#' Update Beta in Gibbs Sampler
#' 
#' Helper function to update beta according to adjacency and node membership
#' @param K A positive integer indicating the true number of clusters
#' @param group Model matrix to be fitted on
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param offset Boolean indicating whether to use offset term
#' 
#' @return beta0 and beta, with beta as a matrix
#' 
#' @note Function \code{update_beta} is a helper for the initialization process 
#' in the \code{calf_sbm_nimble} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
#' 
#' @export

update_beta <- function(K, group, directed = FALSE, offset = FALSE){
    mod_mat <- stats::model.matrix(~ 0 + as.factor(cl), group) * group$x 
    mod_mat2 <- as.data.frame(mod_mat)
    mod_mat2$y <- group$y
    ## Fit Bayesian logistic regression to the data
    logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2)
    if (!offset){
      #logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
      #                  prior.mean = 0, prior.scale = 1)
      logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2)
    } else {
      #logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
      #                  prior.mean = 0, prior.scale = 1, offset = group$offset)
      logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2, 
                       offset = group$offset)
    }
    #predicted_beta <- arm::sim(logit_fit, n.sims = 2) 
    #sampled_beta <- stats::coef(predicted_beta)[1, ]
    #sampled_beta <- coef(logit_fit)
    sampled_beta <- logit_fit$coefficients
    ## Convert beta_kl to matrix form
    beta_mat <- matrix(0, K, K)
    if (!directed) {
        beta_mat[upper.tri(beta_mat, diag = TRUE)] <- sampled_beta[-1]
        beta_mat <- as.matrix(Matrix::forceSymmetric(beta_mat))
    } else {
        beta_mat <- matrix(sampled_beta[-1], K, K)
    }
    return (list(beta0 = sampled_beta[1], beta = beta_mat, aic = logit_fit$aic))
}


#' Update Beta in EM algorithm using logistic regression
#' 
#' Helper function to update beta according to adjacency and node membership
#' @param K A positive integer indicating the true number of clusters
#' @param group Model matrix to be fitted on
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param offset Boolean indicating whether to use offset term
#' @param conf Confidence level of beta estimates (default = 0.95)
#' 
#' @return beta0 and beta, with beta as a matrix
#' 
#' @note Function \code{update_beta_bayes} is a helper in the 
#' in the \code{calfsbm_em} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
#' 
#' @export

update_beta_em <- function(K, group, directed = FALSE, offset = FALSE, conf = 0.95){
  mod_mat <- stats::model.matrix(~ 0 + as.factor(cl), group) * group$x 
  mod_mat2 <- as.data.frame(mod_mat)
  mod_mat2$y <- group$y
  ## Fit Bayesian logistic regression to the data
  if (!offset){
    logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2)
    #logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
    #                  prior.mean = 0, prior.scale = 100)
  } else {
    logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2,
                      offset = group$offset)
  }
  sampled_beta <- logit_fit$coefficients
  ## Inverse Fisher Information for uncertainty quantification
  vars <- diag(vcov(logit_fit))
  fac <- -qnorm((1 - conf) / 2)
  beta0_low <- sampled_beta[1] - fac * sqrt(vars[1])
  beta0_high <- sampled_beta[1] + fac * sqrt(vars[1])
  ## Convert beta_kl to matrix form
  beta_mat <- low <- high <- matrix(0, K, K)
  if (!directed) {
    beta_mat[upper.tri(beta_mat, diag = TRUE)] <- sampled_beta[-1]
    low[upper.tri(low, diag = TRUE)] <- sampled_beta[-1] - fac * sqrt(vars[-1])
    high[upper.tri(high, diag = TRUE)] <- sampled_beta[-1] + fac * sqrt(vars[-1])
    beta_mat <- as.matrix(Matrix::forceSymmetric(beta_mat))
    low <- as.matrix(Matrix::forceSymmetric(low))
    high <- as.matrix(Matrix::forceSymmetric(high))
  } else {
    beta_mat <- matrix(sampled_beta[-1], K, K)
    low <- matrix(sampled_beta[-1] - fac * vars[-1], K, K)
    high <- matrix(sampled_beta[-1] + fac * vars[-1], K, K)
  }
  return (list(beta0 = sampled_beta[1], 
               beta = beta_mat, 
               aic = logit_fit$aic,
               low = list(beta0 = beta0_low, beta = low),
               high = list(beta0 = beta0_high, beta = high)))
}


#' Helper function to update theta offset in EM implementation
update_theta <- function(A, S_ij, z, beta0, beta, theta){
  n <- nrow(A)
  new_theta <- numeric(n)
  for(i in 1:n){
    off <- beta0 + beta[z[i], z[-i]] * S_ij[i, -i] + theta[-i]
    new_theta[i] <- unname(glm(as.vector(A[i, -i]) ~ 1, 
                               family = 'binomial', offset = off)$coefficients[1])
  }
  return(new_theta - mean(new_theta))
}


#' EM algorithm implementation to get variational Bayesian distribution
#' 
#' @param network A list object containing adjacency matrix A, 
#' @param K Number of clusters
#' @param offset Boolean to indicate if there should be an offset (default = TRUE)
#' @param verbose Verbosity (default = TRUE)
#' @param conf Confidence level of betas (default = 0.95)
#' @param S Number of conditional maximization iterations per M-step (default = 1)
#' @return List of estimated node membership, betas, and AIC
#' @export
calfsbm_em <- function(network, K, offset = FALSE, verbose = TRUE, conf = 0.95,
                       S = 1){
  start_time <- Sys.time()
  initial_aic <- Inf; gain <- Inf
  ## SET UP PARAMETERS
  n <- nrow(network$A)
  log_prob_mat <- matrix(0, nrow = n, ncol = K)
  ## INITIALIZE Z AND BETA
  z <- mclust::Mclust(network$X, K, verbose = FALSE)$z
  #z <- matrix(rgamma(n * K, 1), n, K)
  #print(paste0('Initial ARI: ', mclust::adjustedRandIndex(apply(z, 1, which.max), 
  #             network$z)))
  if (offset){
    initial_theta <- update_theta(network$A, network$dis, 
                        apply(z, 1, which.max), 0, matrix(0, K, K), numeric(n))
    group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                        S_ij = network$dis, offset = initial_theta)
  } else {
    group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                        S_ij = network$dis, offset = NULL)
  }
  initial_beta <- update_beta_em(K, group$cluster, offset = offset)
  beta0 <- initial_beta$beta0; beta <- initial_beta$beta
  if(verbose){
    print('Parameters Set!')
    print(paste('Initial AIC:', round(initial_beta$aic, 1)))
  }
  ## BEGIN EM ALGORITHM
  while (gain > 0.000001) {
    ## E-STEP: NODE-LEVEL PROBABILITIES
    initial_z <- apply(z, 1, which.max)
    updated_z <- initial_z
    for (i in 1:n){
      for (j in 1:K){
        eta_ij <- beta0 + beta[j, updated_z[-i]] * network$dis[i, -i]
        ## Calculate probabilities
        fit <- 1 / (1 + exp(-eta_ij))
        ## Calculate log-likelihood
        loglik <- sum(log(ifelse(network$A[i, -i] == 1, fit, 1 - fit)))
        log_prob_mat[i, j] <- loglik
      }
      ## Get relative probabilities of each row
      wts <- log_prob_mat[i, ]; #print(wts)
      wts <- wts - max(wts) #  Safety check for extremely low likelihood
      z[i, ] <- exp(wts)/ sum(exp(wts))
      updated_z[i] <- which.max(wts) # Update node assignment individually
    }
    if(verbose){
      #print(z[n, ])
      print('E STEP COMPLETED!')
    }
    ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
    if (offset){
      ## ECM ITERATIONS
      updated_z <- apply(z, 1, which.max)
      for (s in 1:S){
        new_theta <- update_theta(network$A, network$dis, 
                                  updated_z, initial_beta$beta0, 
                            initial_beta$beta, initial_theta)
        initial_theta <- new_theta
        group <- gen_factor(updated_z, A = network$A, 
                            S_ij = network$dis, offset = new_theta)
        beta_new <- update_beta_em(K, group$cluster, offset = offset, conf = conf)
      }
    } else {
      new_theta <- NULL
      group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                          S_ij = network$dis, offset = NULL)
      beta_new <- update_beta_em(K, group$cluster, offset = offset, conf = conf)
    }
    beta0 <- beta_new$beta0
    beta <- beta_new$beta
    new_aic <- beta_new$aic
    gain <- initial_aic - new_aic
    initial_aic <- new_aic
    if(verbose){
      print('M STEP COMPLETED!')
      print(paste('Current AIC:', round(new_aic, 1)))
    }
  }  
  end_time <- Sys.time()
  return(list(z = apply(z, 1, which.max),
              beta0 = beta0, 
              beta = beta,
              low = beta_new$low,
              high = beta_new$high,
              theta = new_theta,
              aic = new_aic,
              time = as.numeric(difftime(end_time, start_time, units = 's'))))
}


#' (Obsolete) Helper for Gibbs Sampler
#'
#' Update node membership using likelihood
#' @param z Initial values of z
#' @param beta0 Draw from posterior distribution
#' @param beta Draw from posterior distribution
#' @param S_ij Similarity or distance matrix
#' @param A Observed adjacency matrix
#' @param K A positive integer indicating the true number of clusters
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' 
#' @return Updated node membership vector
update_z_from_beta <- function(z, beta0, beta, S_ij, A, K, directed = FALSE){
    n <- length(z)
    log_prob_mat <- matrix(0, nrow = n, ncol = K)
    updated_z <- rep(0, n)
    for (i in 1:n){
        for (j in 1:K){
            eta_ij <- beta0 + beta[j, z[-i]] * S_ij[i, -i]
            ## Calculate probabilities
            fit <- 1 / (1 + exp(-eta_ij))
            ## Calculate log-likelihood
            loglik <- sum(log(ifelse(A[i, -i] == 1, fit, 1 - fit)))
            log_prob_mat[i, j] <- loglik
        }
        ## Get relative probabilities of each row
        wts <- log_prob_mat[i, ]
        wts <- wts - max(wts) #  Safety check for extremely low likelihood
        log_prob_mat[i, ] <- exp(wts)/ sum(exp(wts))
        updated_z[i] <- sample(1:K, 1, prob = log_prob_mat[i, ])
    }
    ## Make cluster assignments from 1 to K in case a cluster was deleted
    serial_z <- updated_z; cl <- 1
    for (k in unique(updated_z)){
        updated_z[serial_z == k] <- cl
        cl = cl + 1
    }
    return(list(updated_z = updated_z, probs = log_prob_mat))
}



#' (Obsolete) Gibbs sampler
#' 
#' Gibbs sampler with fixed number of clusters
#' @param K A positive integer indicating the true number of 
#' clusters (should be optimized already)
#' @param alpha Vector of prior probabilities of being in each group
#' @param beta0 Prior intercept for logistic regression
#' @param beta Prior parameters for logistic regression, in matrix form
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @param directed logical (default = \code{FALSE}); if \code{FALSE}, the 
#' network is undirected
#' 
#' @return List containing the estimated beta, node membership, and history
#' 
#' @note Deprecated
#' 
#' @export
gibbs_sampler_fixed_k <- function(K, alpha, beta0, beta, niter, A, S_ij, 
                                  directed = FALSE){
    n <- nrow(A)
    converged <- FALSE
    ## Initialize beta
    initial_beta <- beta
    ## Initialize node membership
    initial_z <- sample(1:K, n, replace = TRUE, prob = alpha)
#    initial_z <- cluster::pam(links$X, K)$clustering
    ## Initialize X matrix for logistic regression
    group <- gen_factor(initial_z, A, S_ij, directed)
    ## Save node membership and fit for all iterations
    all_z <- matrix(0, nrow = niter, ncol = n)
    prob_history <- array(0, dim = c(niter, n, K))
    bic_history <- rep(0, niter)
    intercepts <- rep(0, niter)
    beta_history <- array(0, dim = c(niter, K, K))
    #i <- 1
    for (i in 1:niter){
        ## Update beta through logistic regression
        initial_beta <- update_beta(K, group$cluster, directed)
        ## Update node membership through likelihood and sampling
        z <- update_z_from_beta(initial_z, initial_beta$beta0, 
                                initial_beta$beta, S_ij, A, K, directed)
        initial_z <- z$updated_z
        ## Update input for update_beta() function
        group$cluster$cl <- (initial_z[group$grid[, 1]] - 1) * K + 
            initial_z[group$grid[, 2]]
        if (!directed) {
            lt <- which(initial_z[group$grid[, 1]] < initial_z[group$grid[, 2]])
            group$cluster$cl[lt] <- (initial_z[group$grid[lt, 2]] - 1) * K + 
                initial_z[group$grid[lt, 1]]
        }
        ## Check for convergence
        #if (i > 5) {
        #    z <- all_z[(i-4):i, ]
        #    if (identical(z[1, ], z[2, ], z[3, ], z[4, ], z[5, ])){
        #        ## Post-convergence sampling
        #        print(paste0('Apparent convergence: ', i))
        #        coda <- convergence_procedure(z, group, niter, i, K)
        #        ## Post-convergence history
        #        intercepts[i:niter] <- coda$ints
        #        bic_history[i:niter] <- coda$bics
        #        beta_history[i:niter, , ] <- coda$beta_his
        #        all_z[i:niter, ] <- coda$z_vals
        #        i <- niter
        #        break
        #    }
        #}
        ## Add to history
        all_z[i, ] <- initial_z
        prob_history[i, , ] <- z$probs
        bic_history[i] <- initial_beta$aic - 2 * K^2 + K^2 * log(n * (n - 1) /2)
        intercepts[i] <- initial_beta$beta0
        beta_history[i, , ] <- initial_beta$beta
        ## Reject sample if likelihood of new sample is too low
        #if (i > 1){
        #  if (runif(1) < exp(bic_history[i - 1] - bic_history[i])) {i <- i + 1}}
        #if (i == 1){i <- i + 1}
        #print(c(i, bic_history[i - 1]))
        if (i %% 100 == 0){print(i)} #  Check to see if still running smoothly
    }
    return(list(beta = initial_beta, z = all_z[niter, ], 
                prob_history = prob_history,#z_history = all_z, 
                bic = bic_history, int_history = intercepts, 
                beta_history = beta_history, 
                converged = i))
}


#' Unknown K Gibbs sampler
#'
#' Gibbs sampler with CRP as initial state
#' @param dp Positive number for CRP concentration parameter (recommended 1)
#' @param niter Total number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @param directed logical; if \code{FALSE} (default), the output is from 
#' an undirected network
#' 
#' @return List containing the estimated beta, node membership, and history
#' 
#' @note Not used for now, but may be used in future versions
gibbs_sampler_unknown_k <- function(dp, niter, A, S_ij, directed = FALSE){
    n <- nrow(A)
    converged <- FALSE
    ## Initialize node membership using Chinese Restaurant Process
    initial_z <- nimble::rCRP(n = 1, conc = dp, n)
    K <- length(unique(initial_z))
    print(K)
    ## Initialize beta
    initial_beta <- diag(K) - 2
    ## Initialize X matrix for logistic regression
    group <- gen_factor(initial_z, A, S_ij, directed)
    ## Save node membership and fit for all iterations
    all_z <- matrix(0, nrow = niter, ncol = n)
    bic_history <- rep(0, niter)
    intercepts <- rep(0, niter)
    #beta_history <- array(0, dim = c(niter, K, K))
    for (i in 1:niter){
      ## Update beta through logistic regression
      K <- length(unique(initial_z))
      print(K)
      initial_beta <- update_beta(K, group$cluster, directed)
      ## Update node membership through likelihood and sampling
      z <- update_z_from_beta(initial_z, initial_beta$beta0, 
                          initial_beta$beta, S_ij, A, K, directed)
      initial_z <- z$updated_z
      ## Update input for update_beta() function
      group$cluster$cl <- (initial_z[group$grid[, 1]] - 1) * K + 
        initial_z[group$grid[, 2]]
      if (!directed) {
          lt <- which(initial_z[group$grid[, 1]] < initial_z[group$grid[, 2]])
          group$cluster$cl[lt] <- (initial_z[group$grid[lt, 2]] - 1) * K + 
              initial_z[group$grid[lt, 1]]
      }
      ## Add to history
      all_z[i, ] <- initial_z
      bic_history[i] <- initial_beta$aic - 2 * K^2 + K^2 * log(n * (n - 1) /2)
      intercepts[i] <- initial_beta$beta0
      #beta_history[i, , ] <- initial_beta$beta
      if (i %% 25 == 0){print(c(i, K))} #  Check to see if still running smoothly
    }
    return(list(beta = initial_beta, z = initial_z, bic = bic_history))
                #history = all_z[1:i, ]))
}


#' (Obsolete) K-selection Function
#'
#' Function that uses BIC to find optimal number of clusters
#' 
#' @param p Number of clusters to test (from 1 to p)
#' @param alpha Vector of prior probabilities of being in each group
#' @param beta0 Prior intercept for logistic regression
#' @param beta Prior parameters for logistic regression, in matrix form
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' 
#' @return List containing the estimated beta, node membership, and history
#' 
#' @note Deprecated
find_k_best_bic <- function(p, alpha, beta0, beta, niter, A, S_ij){
    ## Set a dataframe for BIC and number of clusters
    best_estimate <- data.frame(K = 1:p, BIC = rep(Inf, p), AIC = rep(Inf, p))
    n <- nrow(A)[1]
    ## One cluster case/null model
    null_model <- arm::bayesglm(A[upper.tri(A)] ~ S_ij[upper.tri(S_ij)])
    best_estimate[1, ] <- list(K = 1, AIC = null_model$aic, 
                        BIC = null_model$aic - 2 + log(n * (n - 1) / 2))
    for (i in 2:p){
        i_fit <- gibbs_sampler_fixed_k(i, alpha[1:i], beta0, 
                                       beta[1:i, i:i], niter, A, S_ij)
        fit <- i_fit$bic[niter]
        if (fit < min(best_estimate$BIC)){
            z_est <- i_fit$z
        }
        ## Calculate BIC and AIC for each cluster number
        best_estimate$BIC[i] <- fit
        best_estimate$AIC[i] <- fit + 2 * i^2 - i^2 * log(n * (n - 1)/2)
        print(c(i, fit))
    }
    return(list(K = which.min(best_estimate$BIC), 
                data = best_estimate, z = z_est))
}


#' Function for MFM-SBM Clustering
#'
#' Implementation of Algorithm 1 in MFM-SBM for the model
#' @param z node membership vector
#' @param A Known connectivity matrix
#' @param conc Concentration parameter of CRP
#' @param S_ij Distance
#' @param niter Number of iterations
#' 
#' @return List of beta, z, and history of K
#' 
#' @examples 
#' links <- sim_calfsbm(n_nodes = 50, K = 2, n_covar = 2, 
#'                      prob = c(0.5, 0.5), beta0 = 1, 
#'                      beta = diag(2) - 3, sigma = 0.3, spat = 0.5)
#' mfm_sbm(links$z, links$A, 0.8, links$dis)
#' 
#' @export
mfm_sbm <- function(z, A, conc, S_ij, niter = 100){
    n <- nrow(A)
    ## Initialize node membership using Chinese Restaurant Process
    z <- nimble::rCRP(n = 1, conc, n)
    K <- length(unique(z))
    print(K)
    directed <- !isSymmetric(A)
    group <- gen_factor(z, A, S_ij, directed) 
    k_hist <- rep(0, niter)
    Q <- matrix(0.5, n, n)
    for (iter in 1:niter){
        ## Update beta and Q conditional on z
        beta <- update_beta(K, group$cluster, directed)
        for (i in 1:n){
          Q[i, ] <- nimble::expit(beta$beta0 + 
                                    beta$beta[z[i], z[1:n]] * S_ij[i, ])
        }
        diag(Q) <- 0
        #for(i in 1:K){
        #  for (j in 1:K){
        #    c <- which(z == i)
        #    d <- which(z == j)
        #    Q[i, j] <- mean(expit(beta$beta0 + 
        #      beta$beta[z[c], z[d]] * S_ij[c, d]))
        #    #Q[i, j] <- mean(A[c, d])
        #  }
        #}
        ## Update z conditional on Q
        z <- update_z_from_q(z, Q, A, conc)
        K <- length(unique(z))
        k_hist[iter] <- K
        group$cluster$cl <- (z[group$grid[, 1]] - 1) * K + z[group$grid[, 2]]
    }
    return(list(beta = beta, z = z, k_hist = k_hist))
}


#' Helper to update MFM-SBM by iteration
#'
#' Helper function for MFM-SBM Algorithm 1
#' @param z Vector of node membership, length n 
#' @param Q Stochastic Block Model probability matrix
#' @param A Network adjacency matrix, dimension n x n
#' @param conc Concentration parameter
#' @param a Shape parameter of prior beta distribution
#' @param b Shape parameter of prior beta distribution
#' 
#' @return Updated node membership
#' 
#' @note The function \code{update_z_from_q} is used to update the node 
#' membership in the larger \code{mfm_sbm} function
#' 
#' @keywords Internal
#' 

update_z_from_q <- function(z, Q, A, conc, a = 1, b = 1){
    n <- length(z)
    K <- max(z)
    probs <- matrix(0, nrow = n, ncol = K + 1)
    for (i in 1:(n)){
        new_clus <- rep(0, K)
        for (c in 1:K){
            ## Indices
            j <- (i+1):n; k <- 1:(i-1)
            if(i == n){j <- NULL}
            ## Probability of being in each existing cluster
            scale <- length(which(z == c)) + conc
            #p1 <- Q[c, z[j]]^A[i, j] * (1 - Q[c, z[j]])^(1 - A[i, j])
            #p2 <- Q[z[k], c]^A[k, i] * (1 - Q[z[k], c])^(1 - A[k, i])
            p1 <- Q[i, j]^A[i, j] * (1 - Q[i, j])^(1 - A[i, j])
            p2 <- Q[k, i]^A[k, i] * (1 - Q[k, i])^(1 - A[k, i])
            probs[i, c] <- prod(scale, p1, p2)
            ## Component of new cluster probability
            jj <- i + which(z[j] == c); kk <- which(z[k] == c)
            a_star <- sum(A[i, jj], A[kk, i])
            new_clus[c] <- beta(a + a_star, 
                                b + length(which(z == c)) - a_star) / beta(a, b)
        }
        v <- stats::dpois(K + 1, 1)/stats::dpois(K, 1) * conc
        probs[i, K + 1] <- v * prod(new_clus)
        probs[i, ] <- probs[i, ]/sum(probs[i, ])
    }
    z <- apply(probs, 1, nimble::rcat, n = 1)
    print(table(z))
    return(serialize(z))
}

#' Condense Adjacency into SBM with Known Node Membership
#' 
#' Function to calculate the implied stochastic block probabilities
#' @param A The observed adjacency matrix
#' @param z The true node membership (values assumed to be in [1, K])
#' @return A K x K matrix with the observed density of each block
#' 
#' @examples 
#' set.seed(123)
#' links <- sim_calfsbm(n_nodes = 50, K = 2, n_covar = 2, 
#'                      prob = c(0.5, 0.5), beta0 = 1, beta = diag(2) - 3, 
#'                      sigma = 0.3, spat = 0.5)
#' print(find_sbm(links$A, links$z))
#' 
#' @export
find_sbm <- function(A, z){
  K <- max(z)
  sbm <- matrix(0, K, K)
  for (i in 1:K){
    for (j in 1:K){
      sbm[i, j] <- mean(A[z == i, z == j])
      ## Within-cluster density
      if (i == j){
        n_i <- length(which(z == i))
        sbm[i, j] <- sbm[i, j] * n_i / (n_i - 1)
      }
    }
  }
  return(sbm)
}

#' Post-process raw MCMC output sample-by-sample
#' 
#' Take a large matrix of raw MCMC samples from NIMBLE \code{MCMC}, and 
#' reorder each sample according to label-switching constraint
#' @param mcmcSamples Matrix of raw MCMC samples, with each column corresponding 
#' to each beta value
#' @param K Number of groups
#' @param n Total number of nodes. 
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from an 
#' undirected network
#' 
#' @return Entire MCMC samples data reorganized according to constraint
#' 
#' @examples 
#' \dontrun{
#' set.seed(123)
#' links <- sim_calfsbm(n_nodes = 100, K = 2, n_covar = 2, 
#'                      prob = c(0.5, 0.5), beta0 = 1, 
#'                      beta = diag(2) - 3, sigma = 0.3, spat = 1.5)
#' X <- calf_sbm_nimble(links, 1000, 500, 2, 3, 2)
#' post_label_mcmc_samples(X$mcmcSamples, 2, 1)
#' }
#' @export
post_label_mcmc_samples <- function(mcmcSamples, K, n, directed = FALSE){
  base_inds <- matrix(1:K^2, K, K)
  ## First K^2 columns of mcmcSamples are beta matrix
  new_labels <- label.switching::aic(array(mcmcSamples[, diag(base_inds)], 
                                     dim = c(nrow(mcmcSamples), K, 1)), 
                               constraint = 1)$permutations
  #new_labels <- order(colMeans(mcmcSamples[, diag(base_inds)]))
  ## Make lower diagonal elements same as upper diagonal elements
  if (!directed){
    ut <- which(upper.tri(base_inds))
    lt <- t(base_inds)[which(upper.tri(base_inds))]
    mcmcSamples[, lt] <- mcmcSamples[, ut]
  }
  cols <- ncol(mcmcSamples)
  for (i in 1:nrow(mcmcSamples)){
    a <- factor(mcmcSamples[i, (cols - n) + 1:n])
    row_lab <- new_labels[i, ]
    levels(a) <- row_lab
    ## Reorder beta_{ij} according to label-switching constraint
    mcmcSamples[i, 1:K^2] <- mcmcSamples[i, c(base_inds[row_lab, row_lab])]
    ## Relabel z according to constraint
    mcmcSamples[i, (cols - n) + 1:n] <- as.numeric(levels(a)[as.integer(a)])
  }
  return(mcmcSamples)
}


#' Label-Switching According to Constraint 
#'
#' Function to post-label MCMC samples by ascending order of within-cluster 
#' values. This only orders it according the the means of the output. For 
#' label-switching at the sample level, see \code{post_label_mcmc_samples}
#' @param mcmcSamples Matrix of raw MCMC samples, with each column corresponding 
#' to each beta value
#' @param K Number of groups
#' @param directed logical; if \code{FALSE} (default), the MCMC output is 
#' from an undirected network
#' @param lab logical; if \code{TRUE} (default), labels the clusters from 
#' smallest to largest mean values of within-cluster effects beta_{ii}
#' 
#' @return Reorganized dataframe according to apparent beta_{ii} ordering
#' 
#' @export
post_label_mcmc <- function(mcmcSamples, K, directed = FALSE, lab = TRUE){
  base_inds <- matrix(1:K^2, K, K)
  if (lab){
    new_labels <- order(colMeans(mcmcSamples[, diag(base_inds)]))
  } else {
    new_labels <- 1:K
  }
  if (!directed){
    ut <- which(upper.tri(base_inds))
    lt <- t(base_inds)[which(upper.tri(base_inds))]
    mcmcSamples[, lt] <- mcmcSamples[, ut]
  }
  return(mcmcSamples[, c(base_inds[new_labels, new_labels])])
}


#' Implementation of CALF-SBM 
#' 
#' Using MCMC parameters and number of clusters as input, return raw MCMC output
#' @param adj_mat network adjacency matrix, with 1 indicating a connection and 0
#' indicating no connection. The matrix should be in base R matrix type
#' @param simil_mat network similarity matrix, with each non-diagonal entry
#' indicating the similarity between the two corresponding nodes. It should be 
#' non-negative, and the diagonal elements should all be equal to 0. 
#' @param nsim Total number of MCMC iterations per chain
#' @param burnin Number of iterations in each chain to be discarded
#' @param thin Post-burnin thinning parameter
#' @param nchain Number of MCMC chains to run
#' @param K A positive integer indicating the number of clusters to test on.
#' It is recommended for the user to try a range of potential K values, and
#' choose the K value with the best WAIC.
#' @param covariates A matrix with each row corresponding to a node, and each 
#' column corresponding to a variable. Each entry indicates the value of the 
#' node-specific covariates. If \code{NULL}, initial cluster assignment will be 
#' generated at random
#' @param offset logical (default = \code{TRUE}); where \code{TRUE} 
#' indicates to use offset terms theta in the \code{NIMBLE} model
#' @param beta_scale numeric (default = 10) Prior standard deviation 
#' of all beta terms 
#' @param return_gelman logical (default = \code{FALSE}); if \code{TRUE}, 
#' returns the Gelman-Rubin diagnostic for all beta terms
#' 
#' @return List of beta, z, and history of K
#' 
#' @examples
#' \dontrun{
#' links <- sim_calfsbm(n_nodes = 50, K = 2, n_covar = 2, 
#'                     prob = c(0.5, 0.5), beta0 = 1, beta = diag(2) - 3, 
#'                     sigma = 0.3, spat = 1.5)
#' calf_sbm_nimble(adj_mat = links$A, simil_mat = links$dis, nsim = 1000, 
#'                 burnin = 500, thin = 2, nchain = 2, K = 2, 
#'                 covariates = links$dis, offset = FALSE)
#' }
#' @export
calf_sbm_nimble <- function(adj_mat, simil_mat, nsim, burnin, thin, nchain, K, 
                            covariates, offset = TRUE, beta_scale = 10, 
                            return_gelman = FALSE){
  ## Inits
  const <- list(n = nrow(adj_mat), K = K)
  data <- list(A = adj_mat, x = simil_mat)
  directed <- !isSymmetric(adj_mat)
  inits <- list(beta0 = stats::rnorm(1, 0, 5),
                beta = stats::rnorm(const$K^2, 0, 5),
                gamma = matrix(1, const$n, const$K)
  )
  if (is.null(covariates)){
    inits$z <- sample(1:const$K, const$n, replace = TRUE) 
  } else {
    inits$z <- cluster::pam(covariates, const$K)$clustering
  }
  if (offset){
    if (!directed){
      inits$theta <- log(rowSums(adj_mat) * const$n / sum(adj_mat) + 0.0001)
    } else {
      inits$theta_in <- log(rowSums(adj_mat) * const$n / sum(adj_mat) + 0.0001)
      inits$theta_out <- log(colSums(adj_mat) * const$n / sum(adj_mat) + 0.0001)
    }
  }
  ## Initialize betas
  group <- gen_factor(inits$z, adj_mat, simil_mat)
  initial_beta <- update_beta(const$K, group$cluster)
  inits$beta0 <- initial_beta$beta0
  inits$beta <- c(initial_beta$beta)
  const$beta_scale <- beta_scale
  ## NIMBLE code
  monitors <- c('z', 'beta', 'beta0')
  if(offset){monitors <- c(monitors, 'sigma', 'theta')}
  
  code <- nimble::nimbleCode({
    ## Priors for parameter matrix
    beta0 ~ dnorm(mean = 0, sd = beta_scale)
    for (a in 1:K^2){
      beta[a] ~ dnorm(mean = 0, sd = beta_scale)
    }
    ## Priors for offset    
    if (offset){
      if (!directed) {
        for (i in 1:n){
          theta[i] ~ dnorm(mean = 0, var = sigma)
        }
        sigma ~ dinvgamma(1, 1)
      } else {
        for (i in 1:n){
          theta_in[i] ~ dnorm(mean = 0, var = sigma_in)
          theta_out[i] ~ dnorm(mean = 0, var = sigma_out)
        }
        sigma_in ~ dinvgamma(1, 1)
        sigma_out ~ dinvgamma(1, 1)
      }
    }
    ## Node membership
    for (i in 1:n){
      z[i] ~ dcat(alpha[i, 1:K])
      alpha[i, 1:K] ~ ddirch(gamma[i, 1:K])
    }
    ## Adjacency matrix from fitted values
    for (i in 1:n){
      ## Undirected network
      if (!directed){
        for (j in (i+1):n){
          if (offset){
            A[i, j] ~ dbin(
              expit(beta0 + theta[i] + theta[j] +
                    beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
          } else {
            A[i, j] ~ dbin(
              expit(beta0 + 
                   beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
          }
        }
      } else {
        for (j in 1:n){
          A[i, j] ~ dbin(expit(beta0 + 
                                 theta_in[i] + theta_out[j] + 
                                 beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
        }
      }
    }
  })
  print(const); print(inits)
  ## Compile model
  model <- nimble::nimbleModel(code, 
                               constants = const, 
                               data = data, 
                               inits = inits, 
                               check = FALSE)
  cmodel <- nimble::compileNimble(model)
  ## Compile MCMC sampler
  modelConf <- nimble::configureMCMC(model, monitors = monitors, enableWAIC = TRUE)
  modelMCMC <- nimble::buildMCMC(modelConf)
  cmodelMCMC <- nimble::compileNimble(modelMCMC, project = model) #1 min
  ## Multiple chains runner
  mcmcSamples <- nimble::runMCMC(cmodelMCMC, niter = nsim, nburnin = burnin, 
                                 thin = thin, nchains = nchain)
  ## Return Gelman-Rubin if selected
  if (return_gelman){
    if (nchain == 1) {
      stop('Calculation of Gelman statistic requires more than 1 chain')
    }
    param_names <- colnames(mcmcSamples$chain1)
    gelman.diag <- boa::boa.chain.gandr(
      mcmcSamples, 
      list(mcmcSamples$chain1 - Inf, mcmcSamples$chain1 + Inf), 
      alpha = 0.05, pnames = param_names[grep('beta', param_names)])
  }
  ## Combine mcmcSamples into one matrix
  if (nchain > 1){
    mcmcSamples <- do.call('rbind', mcmcSamples)
  } else {
    mcmcSamples <- as.matrix(mcmcSamples)
  }
  ## Post-process samples using label.switching library
  mcmcSamples <- post_label_mcmc_samples(mcmcSamples, const$K, const$n, directed)
  if (return_gelman){
    return(list(mcmcSamples = mcmcSamples, 
                gelman.diag = gelman.diag, 
                WAIC = cmodelMCMC$getWAIC()))
  } else {
    return(list(mcmcSamples = mcmcSamples, 
                WAIC = cmodelMCMC$getWAIC()))
  }
}


#' Serialize Cluster Labels
#'
#' Number clusters from 1-K in case a cluster was deleted
#' @param z Vector of node membership to be renumbered from 1-K
#' @return Vector of serialized node membership
#' 
#' @examples 
#' x <- serialize(c(2, 2, 3, 3))
#' print(x)
#' 
#' @export
serialize <- function(z){
    z <- as.factor(z)
    levels(z) <- 1:length(levels(z))
    return(as.numeric(levels(z))[z])
}
