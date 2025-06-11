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
#' @importFrom mclust Mclust
#' @importFrom network network
#' @import nimble
NULL


#' Cross-validation function for CALF-SBM variable selection
#' @param network List object containing adjacency, similarity, and covariates
#' @param K Number of clusters
#' @param folds Number of folds (default = 4)
#' @param cluster_intercept Logical (default = \code{FALSE}) indicating 
#' whether the intercept term should be a matrix or a single value
#'  
#' @return BIC of cross-validated model
calfsbm_cv <- function(network, K, folds = 4, cluster_intercept = FALSE){
  n <- nrow(network$X)
  z_est <- numeric(n)
  cv_sample <- sample(1:n)
  ll_total <- 0
  for (i in 1:folds){
    ## Partition of the network
    min_ind <- round((i - 1) * n / folds) + 1
    max_ind <- round(i * n / folds)
    #print(c(min_ind, max_ind))
    train_inds <- cv_sample[-(min_ind:max_ind)]
    #print(train_inds)
    test_inds <- cv_sample[min_ind:max_ind]
    #print(test_inds)
    train_network <- list(A = network$A[train_inds, train_inds], 
                          dis = network$dis[train_inds, train_inds], 
                          X = network$X[train_inds, ])
    test_network <- list(A = network$A[test_inds, test_inds], 
                         dis = network$dis[test_inds, test_inds], 
                         X = network$X[test_inds, ])
    fit <- calfsbm_em(train_network, K, verbose = FALSE, 
                      cluster_intercept = cluster_intercept)
    #print(fit$z); print(fit$beta)
    #print(fit$aic)
    test_z <- c()
    for (i in test_inds){
      z_prob <- numeric(K)
      for (j in 1:K){
        if (!cluster_intercept){
          eta_ij <- fit$beta0 + fit$beta[j, fit$z] * 
            network$dis[i, train_inds]
        } else {
          eta_ij <- fit$beta0[j, fit$z] + fit$beta[j, fit$z] * 
            network$dis[i, train_inds]
        }
        p_hat <- 1 / (1 + exp(-eta_ij))
        loglik <- sum(log(ifelse(network$A[i, train_inds] == 1, 
                                 p_hat, 1 - p_hat)))
        z_prob[j] <- loglik
      }
      #print(z_prob)
      test_z <- c(test_z, which.max(z_prob))
    }
    z_est[test_inds] <- test_z
    m <- expand.grid(test_z, test_z)
    ti <- expand.grid(test_inds, test_inds)
    wuta <- which(ti[, 1] > ti[, 2])
    m <- m[wuta, ]; ti <- ti[wuta, ]
    #print(head(m)); print(head(ti))
    ## Predictions on the holdout
    if (!cluster_intercept){
      etas <- fit$beta0 + fit$beta[m[, 1] * (K - 1) + m[, 2]] * 
        network$dis[ti[, 1] * (n - 1) + ti[, 2]]
    } else {
      etas <- fit$beta0[m[, 1] * (K - 1) + m[, 2]] + 
        fit$beta[m[, 1] * (K - 1) + m[, 2]] * network$dis[ti[, 1] * (n - 1) + ti[, 2]]
    }
    p_hat <- 1 / (1 + exp(-etas))
    ll_pred <- sum(log(ifelse(network$A[ti[, 1] * (n - 1) + ti[, 2]] == 1, 
                              p_hat, 1 - p_hat)))
    ll_total <- ll_total + ll_pred
  }
  ## Return BIC on out-of-sample
  return(-2 * ll_total + (K * (K + 1)) * log(n * (n - 1) / 2))
}


#' EM algorithm implementation to get variational Bayesian distribution
#' 
#' @param network A list object containing adjacency matrix A, 
#' @param K Number of clusters
#' @param cluster_intercept Logical (default = \code{FALSE}) indicating 
#' whether the intercept term should be a matrix or a single value
#' @param verbose Verbosity (default = \code{TRUE})
#' @param find_cov Logical (default = \code{FALSE}) indicating whether the covariance
#' matrix accounting for the uncertainty in the z should be found, will return 
#' observed-data covariance otherwise. The complete variance estimation 
#' can be computationally intensive if chosen
#' @param tol Positive numeric indicating the gain in log-likelihood for which
#' the algorithm will terminate if the improvement is less than
#' 
#' @return List of estimated node membership, betas, and AIC
#' @export
calfsbm_em <- function(network, 
                       K, 
                       cluster_intercept = FALSE,
                       verbose = TRUE, 
                       find_cov = FALSE,
                       tol = 1e-06
                       ){
  start_time <- Sys.time()
  initial_ll <- -Inf; gain <- Inf
  ## SET UP PARAMETERS
  n <- nrow(network$A)
  log_prob_mat <- matrix(0, nrow = n, ncol = K)
  ## INITIALIZE Z AND BETA
  z <- mclust::Mclust(network$X, K, verbose = FALSE)$z
  group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                        S_ij = network$dis, offset = NULL)
  initial_beta <- update_beta_em(K, group$cluster, 
                                 cluster_intercept = cluster_intercept)
  beta0 <- initial_beta$beta0; beta <- initial_beta$beta
  if(verbose){
    print('Parameters Set!')
  }
  ## BEGIN EM ALGORITHM
  while (gain > tol) {
    ## E-STEP: NODE-LEVEL PROBABILITIES
    initial_z <- apply(z, 1, which.max)
    updated_z <- initial_z
    for (i in 1:n){
      for (j in 1:K){
        if (!cluster_intercept){ 
          eta_ij <- beta0 + beta[j, updated_z[-i]] * network$dis[i, -i]
        } else {
          eta_ij <- beta0[j, updated_z[-i]] + 
            beta[j, updated_z[-i]] * network$dis[i, -i]
        }
        ## Calculate probabilities
        fit <- 1 / (1 + exp(-eta_ij))
        ## Calculate log-likelihood
        loglik <- sum(log(ifelse(network$A[i, -i] == 1, fit, 1 - fit)))
        log_prob_mat[i, j] <- loglik
      }
      ## Get relative probabilities of each row
      wts <- log_prob_mat[i, ]; #print(wts)
      wts <- wts - max(wts) #  Safety check for extremely low likelihood
      z[i, ] <- exp(wts) / sum(exp(wts))
      updated_z[i] <- which.max(wts) # Update node assignment individually
    }
    new_ll <- sum(log_prob_mat * z)
    gain <- new_ll - initial_ll
    initial_ll <- new_ll
    if(verbose){
      print('E STEP COMPLETED!')
    }
    ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
    new_theta <- NULL
    group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                          S_ij = network$dis, offset = NULL)
    beta_new <- update_beta_em(K, group$cluster, 
                                 cluster_intercept = cluster_intercept)
    beta0 <- beta_new$beta0
    beta <- beta_new$beta
    if(verbose){
      print('M STEP COMPLETED!')
      print(paste('Current Log-Likelihood:', round(new_ll, 1)))
    }
  }
  if (find_cov) {
    cov_mat <- unname(sem_richardson(network, K, 
      c(mat2vec(beta_new$beta0), mat2vec(beta_new$beta)), z, 
        solve(stats::vcov(beta_new$model)),
        cluster_intercept = cluster_intercept))
  } else {
    cov_mat <- stats::vcov(beta_new$model)
  }
  end_time <- Sys.time()
  return(list(z = apply(z, 1, which.max),
              z_soft = z,
              beta0 = beta0, 
              beta = beta,
              var = cov_mat,
              loglik = new_ll,
              model = beta_new$model,
              time = as.numeric(difftime(end_time, start_time, units = 's'))))
}


#' Function to perform an exact EM algorithm on network data
#' @param network List object containing similarity matrix dis, 
#' @param K Number of clusters
#' @param cluster_intercept Logical (default = \code{FALSE}) indicating 
#' whether the intercept term should be a matrix or a single value
#' @param verbose Logical indicating the verbosity of the function
#' @param find_cov Logical (default = \code{FALSE}) indicating whether the covariance
#' matrix accounting for the uncertainty in the z should be found, will return 
#' observed-data covariance otherwise. The complete variance estimation 
#' can be computationally intensive if chosen
#' @param tol Positive numeric indicating the gain in log-likelihood for which
#' the algorithm will terminate if the improvement is less than
#' 
#' @return Soft clustering estimate z, beta0 and beta estimates, 
#' log-likelihood, GLM output, and computing time
#' @export
calfsbm_em_exact <- function(network, K, cluster_intercept = FALSE, verbose = TRUE, 
                             find_cov = FALSE, tol = 1e-06){
  start_time <- Sys.time()
  initial_ll <- -Inf; gain <- Inf
  ## SET UP PARAMETERS
  n <- nrow(network$A)
  log_prob_mat <- matrix(0, nrow = n, ncol = K)
  ## INITIALIZE Z AND BETA
  z <- mclust::Mclust(network$X, K)$z
  group <- gen_factor_soft(z, A = network$A, 
                           S_ij = network$dis, offset = NULL)
  initial_beta <- update_beta_exact(K, group$cluster, 
                                    cluster_intercept = cluster_intercept)
  beta0 <- initial_beta$beta0; beta <- initial_beta$beta
  if(verbose){
    print('Parameters Set!')
  }
  ## BEGIN EM ALGORITHM
  while (gain > tol) {
    ## E-STEP: NODE-LEVEL PROBABILITIES
    for (i in 1:n){
      loglik <- matrix(0, K, K)
      for (k in 1:K){
        res <- matrix(0, n - 1, K)
        ## Average over all possible configurations
        for (l in 1:K){
          if (cluster_intercept){
            eta_ijkl <- 1 / (1 + exp(-(beta0[k, l] + beta[k, l] * network$dis[i, -i])))
          } else {
            eta_ijkl <- 1 / (1 + exp(-(beta0 + beta[k, l] * network$dis[i, -i])))
          }
          ## Calculate log-likelihood
          loglik[k, l] <- #z[i, k] * 
            sum(log(ifelse(network$A[i, -i] == 1, eta_ijkl, 1 - eta_ijkl)) * z[-i, l])
          res[, l] <- ifelse(network$A[i, -i] == 1, eta_ijkl, 1 - eta_ijkl) * z[-i, l]
          #network$A[i, -i] * log(eta_ijkl) + (1 - network$A[i, -i]) * log(1 - eta_ijkl)
        }
        #log_prob_mat[i, k] <- sum(loglik[k, ])
        log_prob_mat[i, k] <- sum(log(rowSums(res) + 1e-100)) #+ log(z[i, j] + 1e-100)
      }
    }
    for (i in 1:n){
      ## Get relative probabilities of each row
      wts <- log_prob_mat[i, ]
      wts <- wts - max(wts) #  Safety check for extremely low likelihood
      z[i, ] <- exp(wts) / sum(exp(wts))
    }
    print(min(apply(z, 1, max)))
    if(verbose){
      print('E STEP COMPLETED!')
    }
    ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
    group <- gen_factor_soft(z, A = network$A, 
                             S_ij = network$dis, offset = NULL)
    beta_new <- update_beta_exact(K, group$cluster, 
                    init = c(mat2vec(beta0), mat2vec(beta)), 
                    cluster_intercept = cluster_intercept)
    beta0 <- beta_new$beta0
    beta <- beta_new$beta
    ## Update expected log-likelihood
    new_ll <- sum(log_prob_mat * z) + sum(z * log(z + 1e-100))
    gain <- new_ll - initial_ll
    initial_ll <- new_ll
    if(verbose){
      print('M STEP COMPLETED!')
      print(paste('Current Log-Likelihood:', round(new_ll, 1)))
    }
  }
  if (find_cov) {
    cov_mat <- unname(sem_richardson(network, K, 
        c(mat2vec(beta_new$beta0), mat2vec(beta_new$beta)), z, 
        solve(stats::vcov(beta_new$model)), cluster_intercept = cluster_intercept))
  } else {
    cov_mat <- stats::vcov(beta_new$model)
  }
  end_time <- Sys.time()
  return(list(z = z,
              beta0 = beta0, 
              beta = beta,
              loglik = new_ll,
              var = cov_mat,
              time = as.numeric(difftime(end_time, start_time, units = 's')),
              model = beta_new$model))
}


#' Fixed point EM iteration to be used for squarem
#' @param beta description
#' @param network List object containing the following attributes:
#' dis: Similarity matrix
#' A: Adjacency matrix
#' @param z Clustering vector
#' @param cluster_intercept Logical indicating whether the intercept should be 
#' a single value or a matrix corresponding to cluster label
#' @return Vector of updated betas
#' @export
calfsbm_q <- function(beta, network, z, cluster_intercept = FALSE){
  ## E-STEP: NODE-LEVEL PROBABILITIES
  if (!cluster_intercept){
    beta0 <- beta[1]
    beta <- vec2mat(beta[-1])
  } else {
    interc <- 1:(length(beta) / 2)
    beta0 <- vec2mat(beta[interc])
    beta <- vec2mat(beta[-interc])
  }
  initial_z <- apply(z, 1, which.max)
  updated_z <- initial_z
  n <- nrow(network$A); K <- nrow(beta)
  log_prob_mat <- matrix(0, n, K)
  for (i in 1:n){
    for (j in 1:K){
      if (!cluster_intercept){
        eta_ij <- beta0 + 
          beta[j, updated_z[-i]] * network$dis[i, -i]
      } else {
        eta_ij <- beta0[j, updated_z[-i]] + 
          beta[j, updated_z[-i]] * network$dis[i, -i]
      }
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
  ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
  group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                      S_ij = network$dis, offset = NULL)
  beta_new <- update_beta_em(K, group$cluster, 
                             cluster_intercept = cluster_intercept)
  return(c(beta_new$beta0, mat2vec(beta_new$beta)))
}


#' Exact fixed-point EM iteration to be used for squarem
#' @param beta Vector containing initial beta parameter estimates
#' @param network List object containing the following attributes:
#' dis: Similarity matrix
#' A: Adjacency matrix
#' @param z Soft clustering matrix, with each row corresponding to a node,
#' and each column corresponding to a cluster
#' @param cluster_intercept Logical indicating whether the intercept should be 
#' a single value or a matrix corresponding to cluster label
#' 
#' @return Vector of updated betas
#' @export
calfsbm_q_exact <- function(beta, network, z, cluster_intercept = FALSE){
  ## E-STEP: NODE-LEVEL PROBABILITIES
  if (!cluster_intercept){
    beta0 <- beta[1]
    beta <- vec2mat(beta[-1])
  } else {
    interc <- 1:(length(beta) / 2)
    beta0 <- vec2mat(beta[interc])
    beta <- vec2mat(beta[-interc])
  }
  n <- nrow(z); K <- nrow(beta)
  log_prob_mat <- matrix(0, n, K)
  for (i in 1:n){
    loglik <- matrix(0, K, K)
    for (j in 1:K){
      ## Average over all possible configurations
      for (k in 1:K){
        if (!cluster_intercept){
          eta_ij <- beta0 + beta[j, k] * network$dis[i, -i]
        } else {
          eta_ij <- beta0[j, k] + beta[j, k] * network$dis[i, -i]
        }
        ## Calculate probabilities
        fit <- 1 / (1 + exp(-eta_ij))
        ## Calculate log-likelihood
        loglik[j, k] <- sum(log(ifelse(network$A[i, -i] == 1, fit, 1 - fit)) * z[-i, k])
      }
      log_prob_mat[i, j] <- sum(loglik[j, ])
    }
    wts <- log_prob_mat[i, ]; #print(wts)
    wts <- wts - max(wts) #  Safety check for extremely low likelihood
    z[i, ] <- exp(wts) / sum(exp(wts))
  }
  ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
  group <- gen_factor_soft(z, A = network$A, 
                           S_ij = network$dis, offset = offset)
  if (!cluster_intercept){
    beta_new <- update_beta_exact(K, group$cluster, init = c(beta0, mat2vec(beta)), 
                                  cluster_intercept = cluster_intercept)
  } else {
    beta_new <- update_beta_exact(K, group$cluster, 
      init = c(mat2vec(beta0), mat2vec(beta)), cluster_intercept = cluster_intercept)
  }
  beta0 <- beta_new$beta0
  beta <- beta_new$beta
  #print(beta0); print(mat2vec(beta))
  return(c(mat2vec(beta_new$beta0), mat2vec(beta_new$beta)))
}


#' Function to perform cross-validation on the CAMM-SBM model. This function can  
#' be used to select the number of clusters or to perform feature selection
#' @param network List object containing attributes
#' @param K List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param folds Number of folds of cross-validation (default = 5)
#' @param mcmc Boolean indicating whether to run the model in MCMC format.
#' The model will run much slower, but may be more accurate
#' @return Combined BIC of all the test sets
#' @export
cammsbm_cv <- function(network, K, folds = 5, mcmc = FALSE){
  n <- nrow(network$X)
  z_est <- numeric(n)
  cv_sample <- sample(1:n)
  ll_total <- 0
  for (i in 1:folds){
    ## Partition of the network
    min_ind <- round((i - 1) * n / folds) + 1
    max_ind <- round(i * n / folds)
    #print(c(min_ind, max_ind))
    train_inds <- cv_sample[-(min_ind:max_ind)]
    test_inds <- cv_sample[min_ind:max_ind]
    train_network <- list(A = network$A[train_inds, train_inds], 
                          dis = network$dis[train_inds, train_inds], 
                          X = network$X[train_inds, ])
    test_network <- list(A = network$A[test_inds, test_inds], 
                         dis = network$dis[test_inds, test_inds], 
                         X = network$X[test_inds, ])
    if (!mcmc){
      fit <- calfsbm_em(train_network, K, verbose = FALSE, 
                        cluster_intercept = TRUE)
    } else {
      fit <- cammsbm_nimble(train_network, K)
    }
    test_z <- c()
    for (i in test_inds){
      z_prob <- numeric(K)
      for (j in 1:K){
        if (!cluster_intercept){
          eta_ij <- fit$beta0 + fit$beta[j, fit$z] * 
            network$dis[i, train_inds]
        } else {
          eta_ij <- fit$beta0[j, fit$z] + fit$beta[j, fit$z] * 
            network$dis[i, train_inds]
        }
        p_hat <- 1 / (1 + exp(-eta_ij))
        loglik <- sum(log(ifelse(network$A[i, train_inds] == 1, 
                                 p_hat, 1 - p_hat)))
        z_prob[j] <- loglik
      }
      test_z <- c(test_z, which.max(z_prob))
    }
    z_est[test_inds] <- test_z
    m <- expand.grid(test_z, test_z)
    ti <- expand.grid(test_inds, test_inds)
    wuta <- which(ti[, 1] > ti[, 2])
    m <- m[wuta, ]; ti <- ti[wuta, ]
    ## Predictions on the holdout
    if (!cluster_intercept){
      etas <- fit$beta0 + fit$beta[m[, 1] * (K - 1) + m[, 2]] * 
        network$dis[ti[, 1] * (n - 1) + ti[, 2]]
    } else {
      etas <- fit$beta0[m[, 1] * (K - 1) + m[, 2]] + 
        fit$beta[m[, 1] * (K - 1) + m[, 2]] * network$dis[ti[, 1] * (n - 1) + ti[, 2]]
    }
    p_hat <- 1 / (1 + exp(-etas))
    ll_pred <- sum(log(ifelse(network$A[ti[, 1] * (n - 1) + ti[, 2]] == 1, 
                              p_hat, 1 - p_hat)))
    ll_total <- ll_total + ll_pred
  }
  ## Return BIC on out-of-sample
  return(-2 * ll_total + (K * (K + 1)) * log(n * (n - 1) / 2))
}


#' Function to run MCMC version of CAMM-SBM model
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param K Number of clusters
#' @param nsim Number of simulations
#' @param burnin Number of burn-in MCMC iterations to discard
#' @param thin Thinning rate of post-burn-in iterations
#' @param nchain Number of chains to run
#' @param beta_scale Standard deviation for the beta parameters
#' @return WAIC of the fitted model
cammsbm_nimble <- function(network, K,
                           nsim = 2000, burnin = 1000, thin = 2, 
                           nchain = 2, beta_scale = 10){
  adj_mat <- network$A; simil_mat <- network$dis; covariates <- network$X
  ## Inits
  const <- list(n = nrow(adj_mat), K = K)
  data <- list(A = adj_mat, x = simil_mat)
  directed <- !isSymmetric(adj_mat)
  inits <- list(beta0 = stats::rnorm(const$K^2, 0, 5),
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
    for (a in 1:K^2){
      beta0[a] ~ dnorm(mean = 0, sd = beta_scale)
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
              expit(beta0[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] + theta[i] + theta[j] +
                      beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
          } else {
            A[i, j] ~ dbin(
              expit(beta0[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] + 
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
  
  ## Combine mcmcSamples into one matrix
  if (nchain > 1){
    mcmcSamples <- do.call('rbind', mcmcSamples)
  } else {
    mcmcSamples <- as.matrix(mcmcSamples)
  }
  ## Post-process samples using label.switching library
  mcmcSamples <- post_label_mcmc_samples(mcmcSamples, const$K, const$n, directed)
  beta0 <- Matrix::forceSymmetric(matrix(colMeans(mcmcSamples[, 1:K^2]), K, K))
  beta <- Matrix::forceSymmetric(matrix(colMeans(mcmcSamples[, K^2 + 1:K^2]), K, K))
  z <- numeric(const$n)
  for (i in 1:const$n){
    z[i] <- which.max(table(mcmcSamples[, (2 * K^2) + i]))
  }
  return(list(mcmcSamples = mcmcSamples, 
              WAIC = cmodelMCMC$getWAIC(),
              beta0 = beta0,
              beta = beta,
              z = z))
}


#' Variational EM implementation for CAMM-SBM
#' @param network List object containing adjacency (A), similarity (dis) and 
#' covariates (X)
#' @param K Number of clusters to include
#' @param verbose Boolean for verbosity
#' @param tol ELBO gain threshold (less than this
#' gain will stop the gain)
#' @return Clustering probabilities, beta parameters, along with 
#' diagnostics like iterations, ELBO values, and computing time.
#' @export
cammsbm_vem <- function(network, K, verbose = TRUE, tol = 1e-03){
  start_time <- Sys.time()
  initial_elbo <- -Inf; gain <- Inf
  elbo_values <- c()
  ## SET UP PARAMETERS
  n <- nrow(network$A)
  log_prob_mat <- matrix(0, nrow = n, ncol = K)
  ## INITIALIZE Z AND BETA
  z <- mclust::Mclust(network$X, K)$z
  group <- gen_factor_soft(z, A = network$A, S_ij = network$dis, offset = NULL)
  initial_beta <- update_beta_exact(K, group$cluster, 
                                    cluster_intercept = TRUE)
  #pseudo_em_result <- calfsbm_em(network, K, cluster_intercept = TRUE)
  #z <- pseudo_em_result$z_soft
  beta0 <- initial_beta$beta0; beta <- initial_beta$beta
  if(verbose){
    print('Parameters Set!')
  }
  ## BEGIN EM ALGORITHM
  while (gain > tol) {
    ## E-STEP: NODE-LEVEL PROBABILITIES
    for (i in 1:n){
      loglik <- matrix(0, K, K)
      for (k in 1:K){
        res <- matrix(0, n - 1, K)
        ## Average over all possible configurations
        for (l in 1:K){
          if (cluster_intercept){
            eta_ijkl <- 1 / (1 + exp(-(beta0[k, l] + beta[k, l] * network$dis[i, -i])))
          } else {
            eta_ijkl <- 1 / (1 + exp(-(beta0 + beta[k, l] * network$dis[i, -i])))
          }
          ## Calculate log-likelihood
          loglik[k, l] <-
            sum(log(ifelse(network$A[i, -i] == 1, eta_ijkl, 1 - eta_ijkl)) * z[-i, l])
          res[, l] <- ifelse(network$A[i, -i] == 1, eta_ijkl, 1 - eta_ijkl) * z[-i, l]
        }
        log_prob_mat[i, k] <- sum(log(rowSums(res) + 1e-100)) #+ log(z[i, j] + 1e-100)
      }
      wts <- log_prob_mat[i, ]
      wts <- wts - max(wts) #  Safety check for extremely low likelihood
      z[i, ] <- exp(wts) / sum(exp(wts))
    }
    #print(min(apply(z, 1, max)))
    if(verbose){
      print('E STEP COMPLETED!')
      #print(find_elbo(c(mat2vec(beta0), mat2vec(beta)), network, z, cluster_intercept))
    }
    ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
    group <- gen_factor_soft(z, A = network$A, 
                             S_ij = network$dis, offset = NULL)
    beta_new <- update_beta_exact(K, group$cluster, 
                                  init = c(mat2vec(beta0), mat2vec(beta)), 
                                  cluster_intercept = cluster_intercept)
    beta0 <- beta_new$beta0
    beta <- beta_new$beta
    ## Update expected log-likelihood
    new_elbo <- -find_elbo(c(mat2vec(beta0), mat2vec(beta)), network, z, cluster_intercept)
    gain <- new_elbo - initial_elbo
    initial_elbo <- new_elbo
    elbo_values <- c(elbo_values, new_elbo)
    if(verbose){
      print('M STEP COMPLETED!')
      print(paste('Current ELBO:', round(new_elbo, 3)))
      #print(paste('Current ELBO:', round(new_elbo, 1)))
    }
  }
  if (find_cov) {
    cov_mat <- unname(sem_richardson(network, K, 
                                     c(mat2vec(beta_new$beta0), mat2vec(beta_new$beta)), z, 
                                     solve(stats::vcov(beta_new$model)), cluster_intercept = cluster_intercept))
  } else {
    cov_mat <- stats::vcov(beta_new$model)
  }
  end_time <- Sys.time()
  return(list(z = z,
              beta0 = beta0, 
              beta = beta,
              elbo = new_elbo,
              var = cov_mat,
              time = as.numeric(difftime(end_time, start_time, units = 's')),
              iter = length(elbo_values),
              evals = elbo_values,
              model = beta_new$model))
}


#' Function to find the ELBO for the CAMM-SBM model. As currently constructed,
#' the function will actually report the negative of the ELBO to be compatible
#' with the optim() function, which only minimizes objective functions
#' @param beta_vec beta_0 an beta parameters flattened into 
#' a combined vector
#' @param network 
#' @param z Matrix of clustering probabilities found in the E-Step
#' @return Observed ELBO with the given parameters
find_elbo <- function(beta_vec, network, z, cluster_intercept = TRUE){
  K <- ncol(z); n <- nrow(z)
  interc <- 1:(length(beta_vec) / 2)
  beta0 <- vec2mat(beta_vec[interc])
  beta <- vec2mat(beta_vec[-interc])
  elbo <- 0
  for (i in 1:n){
    loglik <- matrix(0, K, K)
    for (k in 1:K){
      res <- matrix(0, n - 1, K)
      ## Average over all possible configurations
      for (l in 1:K){
        if (cluster_intercept){
          eta_ijkl <- 1 / (1 + exp(-(beta0[k, l] + beta[k, l] * network$dis[i, -i])))
        } else {
          eta_ijkl <- 1 / (1 + exp(-(beta0 + beta[k, l] * network$dis[i, -i])))
        }
        elbo <- elbo + 
          z[i, k] * sum(z[-i, l] * log(ifelse(network$A[i, -i] == 1, eta_ijkl, 1 - eta_ijkl)))
      }
    }
  }
  return(-elbo)
}


#' Function to perform forward stepwise variable selection
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param K Number of clusters
#' @param cluster_intercept Whether to use a matrix of beta_0 values 
#' (default = FALSE)
#' @param cv Whether to use cross-validation (default = TRUE)
#' 
#' @return Vector of covariate indices which correspond to the optimal model
#' @export
forward_stepwise_em <- function(network, K, cluster_intercept = FALSE, 
                                cv = TRUE) {
  current_bic <- Inf
  remaining_predictors <- 1:ncol(network$X)  # List of all predictor names
  selected_predictors <- c(1)          # List to store selected predictors
  while(length(remaining_predictors) > 0) {
    bic_values <- numeric(length(remaining_predictors))
    print(remaining_predictors)
    for (i in 1:length(remaining_predictors)) {
      predictor <- remaining_predictors[i]
      S_ij <- dist(network$X[, c(selected_predictors[-1], predictor)], 
                   upper = TRUE, diag = TRUE)
      S_ij <- as.matrix(Matrix::forceSymmetric(as.matrix(S_ij)))
      diag(S_ij) <- 0
      network$dis <- S_ij
      if (cv){
        bic_values[i] <- calfsbm_cv(network, K, 
                                    cluster_intercept = cluster_intercept)
      } else {
        n <- nrow(network$A)
        bic_values[i] <- calfsbm_em(network, K, 
                                    cluster_intercept = cluster_intercept, 
                                    verbose = FALSE)$loglik + 
          (K^2 + K) * log(n * (n - 1) / 2)
      }
    }
    print(bic_values)
    best_index <- which.min(bic_values)
    if (bic_values[best_index] < current_bic) {
      selected_predictors <- c(selected_predictors, remaining_predictors[best_index])
      remaining_predictors <- remaining_predictors[-best_index]
      current_bic <- bic_values[best_index]  # Update the current BIC
    } else {
      break
    }
  }
  return(selected_predictors[-1])
}


#' Function to set up Q function in exact EM form
#' @param initial_z n-by-K matrix of node-level clustering probabilities
#' @param A n-by-n adjacency matrix
#' @param S_ij n-by-n similarity matrix
#' @param directed Logical indicating whether the network should be directed
#' @param offset Logical indicating whether to include node-level 
#' random-effects terms (default = FALSE)
#' 
#' @return Data frame with weights, combinations of connections and clusters
#' @note Internal helper for \code{calfsbm_em_exact} and \code{calfsbm_q_exact}
gen_factor_soft <- function(initial_z, A, S_ij, directed = FALSE, offset = FALSE){
  ## Initialize number of nodes/clusters
  n <- nrow(initial_z)
  K <- ncol(initial_z)
  eg <- expand.grid(1:n, 1:n, 1:K, 1:K)
  ## Organize outcomes, distance, and weights into dataframe
  cluster <- data.frame(y = A[(eg[, 1] - 1) * n + eg[, 2]],
                        x = S_ij[(eg[, 1] - 1) * n + eg[, 2]],
                        wts = initial_z[(eg[, 3] - 1) * n + eg[, 1]] * 
                              initial_z[(eg[, 4] - 1) * n + eg[, 2]],
                        cl = (eg[, 3] - 1) * K + eg[, 4]
  )
  wuta <- which((eg[, 1] > eg[, 2]) & (eg[, 3] >= eg[, 4])) 
  #print(head(cluster$wts, 50))
  return (list(cluster = cluster[wuta, ], grid = eg[wuta, ]))
}


#' Helper function to switch from a vector to a matrix representation of beta
#' @param beta Vector of entries
#' @param directed Logical indicating whether network is directed
#'
#' @return Matrix of entries
#' @note Small helper function
mat2vec <- function(beta, directed = FALSE){
  if (directed){ return(as.vector(beta[upper.tri(beta, diag = TRUE)])) }
  else { return(beta[upper.tri(beta, diag = TRUE)]) }
}


#' Function to run the SEM Algorithm to estimate variance of parameters
#' 
#' @param network List object containing similarity matrix dis, 
#' adjacency matrix A, and matrix of covariates X
#' @param K Number of clusters
#' @param beta Fitted beta values from the EM algorithm
#' @param z Fitted matrix of z probabilities
#' @param obs_fisher Observed Fisher information from the EM fit
#' @param increment Small increment to calculate numerical Jacobian 
#' @return Covariance matrix of betas
#' @note Internal helper
sem <- function(network, K, beta, z, obs_fisher, increment = 0.0001){
  m <- length(beta)
  DM <- matrix(0, m, m)
  ## ITERATE OVER BETA_T VECTOR
  for (i in 1:m){
    u_i <- ifelse(1:m == i, increment, 0)
    ## FULL ITERATION OF EM FROM MODIFIED BETA VALUES
    new_beta <- calfsbm_q_exact(beta + u_i, network, z)
    #print(new_beta)
    DM[i, ] <- (new_beta - beta) / (increment)
  }
  V <- solve(obs_fisher) %*% solve(diag(m) - DM)
  return(V)
}


#' SEM Algorithm using Richardson extrapolation
#' 
#' @param network List object containing similarity matrix dis, 
#' adjacency matrix A, and matrix of covariates X
#' @param K Number of clusters
#' @param beta Fitted beta values from the EM algorithm, as vector
#' @param z Fitted matrix of z probabilities
#' @param obs_fisher Observed Fisher information from the EM fit
#' @param cluster_intercept Logical indicating whether the intercept should be 
#' a single value or a matrix corresponding to cluster label
#' @param increment Small increment to calculate numerical Jacobian 
#' @return Covariance matrix of betas
#' @note Internal helper
sem_richardson <- function(network, K, beta, z, obs_fisher, 
                           cluster_intercept = FALSE, increment = 0.0001){
  m <- length(beta)
  DM <- matrix(0, m, m)
  ## ITERATE OVER BETA_T VECTOR
  for (j in 1:m){
    u_j <- ifelse(1:m == j, increment, 0)
    ## FULL ITERATION OF EM FROM MODIFIED BETA VALUES
    beta_h <- calfsbm_q_exact(beta + u_j, network, z, cluster_intercept)
    beta_2h <- calfsbm_q_exact(beta + 2 * u_j, network, z, cluster_intercept)
    beta_nh <- calfsbm_q_exact(beta - u_j, network, z, cluster_intercept)
    beta_n2h <- calfsbm_q_exact(beta - 2 * u_j, network, z, cluster_intercept)
    DM[j, ] <- (beta_n2h - 8 * beta_nh + 8 * beta_h - beta_2h) / 
      (12 * increment)
  }
  V <- solve(obs_fisher) %*% solve(diag(m) - DM)
  return(V)
}


#' Stochastic derivative code, tailored for the calfsbm model
#' @param g Fixed-point function
#' @param x Parameters
#' @param n Noise
#' @param B Number of iterations
#' @param s Noise parameter
#' @param network Network
#' @param z_soft Soft clustering matrix
#' 
#' @return Sampled Jacobian matrix
#' @note Helper
stoxd <- function(g, x, n = 100, B = 100, s = 100, network, z_soft) {
  p <- length(x)
  z <- matrix(rnorm(p * B) / sqrt(n) / s, B, p)
  gxz <- y <- z
  gx <- g(x, network = network, z = z_soft)
  for (i in 1:B) {
    gxz[i, ] <- g(x + z[i, ], network = network, z = z_soft)
    y[i, ] <- gxz[i, ] - gx
  }
  der <- matrix(0, p, p)
  for (i in 1:p) {
    der[,i] <- coef(lm.fit(z, y[,i]))
  }
  der
}


#' Helper function to switch from a vector to a matrix representation of beta
#' @param beta Vector of entries
#' @param directed Logical indicating whether network is directed
#'
#' @return Matrix of entries
vec2mat <- function(beta, directed = FALSE){
  if (directed) {
    K <- floor(sqrt(length(beta)))
    return (matrix(beta, K, K))
  } else {
    K <- as.integer(sqrt(2 * length(beta) + 0.25) - 0.5)
    mat <- matrix(0, K, K)
    mat[upper.tri(mat, diag = TRUE)] <- beta
    return (as.matrix(Matrix::forceSymmetric(mat)))
  }
}


#' Update Beta in EM algorithm using logistic regression
#' 
#' Helper function to update beta according to adjacency and node membership
#' @param K A positive integer indicating the true number of clusters
#' @param group Model matrix to be fitted on
#' @param init Starting values for logistic regression (default = \code{NULL})
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param cluster_intercept Logical indicating whether the intercept term
#' should be a matrix or a single value
#' 
#' @return beta0 and beta, with beta as a matrix
#' 
#' @note Function \code{update_beta_em} is a helper in the 
#' in the \code{calfsbm_em} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
update_beta_em <- function(K, group, init = NULL, directed = FALSE, 
                           cluster_intercept = FALSE){
  mod_mat <- stats::model.matrix(~ 0 + as.factor(cl), group) * group$x 
  if (cluster_intercept){
    int_mat <- unname(stats::model.matrix(~ 0 + as.factor(cl), group))
    mod_mat <- cbind(int_mat, mod_mat)
  }
  mod_mat2 <- as.data.frame(mod_mat)
  mod_mat2$y <- group$y
  if (!cluster_intercept){
    logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2, start = init)
    
  } else {
    logit_fit <- glm(y ~ 0 + . , family = 'binomial', data = mod_mat2, start = init)
  }
  sampled_beta <- logit_fit$coefficients
  ## Intercept term indexing
  if (cluster_intercept) {interc <- 1:(K * (K + 1) / 2); } 
  else {interc <- 1}
  ## Convert beta_kl to matrix form
  if (!directed) {
    beta_mat <- vec2mat(sampled_beta[-interc])
  } else {
    beta_mat <- matrix(sampled_beta[-interc], K, K)
  }
  if(cluster_intercept){
    intercept <- vec2mat(sampled_beta[interc])
  } 
  else {intercept <- sampled_beta[1]}
  return (list(beta0 = intercept, 
               beta = beta_mat, 
               aic = logit_fit$aic,
               model = logit_fit))
}


#' Update Beta in exact EM form
#' 
#' Helper function to update beta according to adjacency and node membership
#' @param K A positive integer indicating the true number of clusters
#' @param group Model matrix to be fitted on
#' @param init Starting values for logistic regression (default = \code{NULL})
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param offset Logical indicating whether to use offset term
#' @param cluster_intercept Logical indicating whether the intercept term
#' should be a matrix or a single value
#' 
#' @return List object with fitted beta0 and beta, along with model diagnostics
#' 
#' @note Function \code{update_beta_exact} is a helper in the 
#' in the \code{calfsbm_em_exact} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
update_beta_exact <- function(K, group, init = NULL, directed = FALSE, 
                              offset = FALSE, cluster_intercept = FALSE){
  mod_mat <- stats::model.matrix(~ 0 + as.factor(cl), group) * group$x 
  if (cluster_intercept){
    int_mat <- unname(stats::model.matrix(~ 0 + as.factor(cl), group))
    mod_mat <- cbind(int_mat, mod_mat)
  }
  mod_mat2 <- as.data.frame(mod_mat)
  mod_mat2$y <- group$y
  mod_mat2$wts <- group$wts
  ## Maximize likelihood through logistic regression
  if (cluster_intercept) {
    iter <- 1
    intercept <- beta_mat <- rep(0, K * (K + 1) / 2)
    ## Iterate sequentially over all beta0 and beta terms
    for (i in unique(group$cl)){
      cols <- c(iter, iter + K * (K + 1) / 2, K * (K + 1) + 1:2)
      subset <- mod_mat2[which(mod_mat2[, iter] == 1), cols]
      logit_fit <- suppressWarnings(stats::glm(y ~ 0 + . - wts, family = 'binomial', 
            data = subset, weights = wts, start = init[cols[1:2]]))
      intercept[iter] <- logit_fit$coefficients[1]
      beta_mat[iter] <- logit_fit$coefficients[2]
      iter <- iter + 1
    }
    intercept <- vec2mat(intercept)
    beta_mat <- vec2mat(beta_mat)
  } else {
    logit_fit <- suppressWarnings(stats::glm(y ~ . - wts, family = 'binomial', 
                                             data = mod_mat2, weights = wts, start = init))
    sampled_beta <- logit_fit$coefficients
    beta_mat <- low <- high <- matrix(0, K, K)
    interc <- 1
    beta_mat <- matrix(sampled_beta[-interc], K, K)
    intercept <- sampled_beta[1]
  }
  return (list(beta0 = intercept,
               beta = beta_mat, 
               aic = logit_fit$aic,
               model = logit_fit))
}

