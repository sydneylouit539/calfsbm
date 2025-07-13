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
#' @importFrom mclust Mclust mclustBIC
#' @importFrom network network
#' @import nimble
#' @importFrom stats binomial coef dist lm.fit offset rnorm runif
#' @importFrom utils flush.console
NULL


#' Cross-validation function for CAMM-SBM variable selection
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param K Number of clusters
#' @param folds Number of folds (default = 5)
#'  
#' @return Predictive pseudo-likelihood of cross-validated model
cammsbm_cv <- function(network, K, folds = 5){
  n <- nrow(network$X)
  z_est <- numeric(n)
  cv_sample <- sample(1:n)
  ll_total <- 0
  for (i in 1:folds){
    ## Partition of the network
    min_ind <- round((i - 1) * n / folds) + 1
    max_ind <- round(i * n / folds)
    train_inds <- cv_sample[-(min_ind:max_ind)]
    test_inds <- cv_sample[min_ind:max_ind]
    train_network <- list(A = network$A[train_inds, train_inds], 
                          dis = network$dis[train_inds, train_inds], 
                          X = network$X[train_inds, ])
    test_network <- list(A = network$A[test_inds, test_inds], 
                         dis = network$dis[test_inds, test_inds], 
                         X = network$X[test_inds, ])
    fit <- cammsbm_plem(train_network, K, verbose = FALSE)
    test_z <- c()
    for (i in test_inds){
      z_prob <- numeric(K)
      for (j in 1:K){
        eta_ij <- fit$beta0[j, fit$z] + fit$beta[j, fit$z] * 
            network$dis[i, train_inds]
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
    etas <- fit$beta0[m[, 1] * (K - 1) + m[, 2]] + 
      fit$beta[m[, 1] * (K - 1) + m[, 2]] * network$dis[ti[, 1] * (n - 1) + ti[, 2]]
    p_hat <- 1 / (1 + exp(-etas))
    p_hat <- min(max(p_hat, 1e-8), 1 - 1e-8)
    ll_pred <- sum(log(ifelse(network$A[ti[, 1] * (n - 1) + ti[, 2]] == 1, 
                              p_hat, 1 - p_hat)))
    ll_total <- ll_total + ll_pred
  }
  ## Return predictive out-of-sample log-pseudo-likelihood
  return(ll_total)
}


#' Pseudo-likelihood EM algorithm implementation
#' 
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param K Number of clusters
#' @param max_iter Numeric (default = 20) indicating 
#' the maximum number of iterations to run (in the case of clustering 
#' configuration alternating between iterations)
#' @param tol Positive numeric indicating the gain in log-likelihood for which
#' the algorithm will terminate if the improvement is less than
#' @param verbose Verbosity (default = \code{TRUE})
#' 
#' @return List of estimated node membership, betas, and log-pseudo-likelihood
#' @export
cammsbm_plem <- function(network, K, max_iter = 20, tol = 1e-4, verbose = TRUE) {
  n <- nrow(network$A)
  A <- network$A
  S <- network$dis
  ## Initialize cluster assignments randomly
  z <- apply(mclust::Mclust(network$X, G = K, verbose = FALSE)$z, 1, which.max)
  
  ## Initialize beta parameters
  beta_0 <- matrix(0, K, K)
  beta_1 <- matrix(0, K, K)
  prev_lpl <- -Inf  
  
  start_time <- Sys.time() # Timer
  for (iter in 1:max_iter) {
    if (verbose){ cat("Iteration", iter, "\n")}
    ## M-step: Fit logistic regressions
    for (k in 1:K) {
      for (l in k:K) {
        i_idx <- which(z == k)
        j_idx <- which(z == l)
        if (length(i_idx) == 0 || length(j_idx) == 0) {next}
        ## Get unique unordered pairs (i < j)
        pairs <- t(combn(c(i_idx, j_idx), 2))
        pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
        valid_pairs <- pairs[
          (z[pairs[, 1]] == k & z[pairs[, 2]] == l) |
            (z[pairs[, 1]] == l & z[pairs[, 2]] == k), , drop = FALSE
        ]
        if (nrow(valid_pairs) == 0) {next}
        y <- A[valid_pairs]
        x <- S[valid_pairs]
        df <- data.frame(y = y, x = x)
        start_vals <- c(beta_0[k, l], beta_1[k, l]) # Warm start
        fit <- tryCatch(
          glm(y ~ x, family = binomial(), data = df, start = start_vals), 
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          beta_0[k, l] <- coef(fit)[1]
          beta_1[k, l] <- coef(fit)[2]
          beta_0[l, k] <- beta_0[k, l]
          beta_1[l, k] <- beta_1[k, l]
        } else {
          beta_0[k, l] <- beta_1[k, l] <- 0
          beta_0[l, k] <- beta_1[l, k] <- 0
        }
      }
    }
    ## E-step: Update cluster assignments
    z_new <- z
    for (i in 1:n) {
      log_likelihoods <- rep(0, K)
      for (k in 1:K) {
        ll <- 0
        for (j in setdiff(1:n, i)) {
          l <- z[j]
          eta <- plogis(beta_0[k, l] + beta_1[k, l] * S[i, j])
          eta <- min(max(eta, 1e-6), 1 - 1e-6)  # numerical safety
          ll <- ll + A[i, j] * log(eta) + (1 - A[i, j]) * log(1 - eta)
        }
        log_likelihoods[k] <- ll
      }
      z_new[i] <- which.max(log_likelihoods)
    }
    z <- z_new
    ## LPL computation
    lpl <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- z[i]; l <- z[j]
        eta <- plogis(beta_0[k, l] + beta_1[k, l] * S[i, j])
        eta <- min(max(eta, 1e-8), 1 - 1e-8)
        lpl <- lpl + A[i, j] * log(eta) + (1 - A[i, j]) * log(1 - eta)
      }
    }
    delta_lpl <- abs(lpl - prev_lpl)
    ## Check for NaN/NA or non-finite values
    if (iter == 1) {
      delta_lpl <- Inf  # no convergence check on first iteration
    } else {
      delta_lpl <- abs(lpl - prev_lpl)
    }
    if (!is.finite(lpl) || (!is.finite(delta_lpl) && iter > 1)) {
      warning("Non-finite LPL or delta encountered - exiting early.")
      break
    }
    if (verbose) {
      cat("  Log-pseudo-likelihood:", round(lpl, 4), "\n")
    }
    if (iter > 1 && delta_lpl < tol) {
      if (verbose) {cat("Converged based on LPL.\n")}
      break
    }
    prev_lpl <- lpl
  }
  end_time <- Sys.time()
  return(list(z = z, beta0 = beta_0, beta = beta_1, loglik = lpl,
              time = as.numeric(difftime(end_time, start_time, units = 's'))))
}


#' Function to perform an exact EM algorithm on network data
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X) 
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
#' @note Internal
cammsbm_plem_exact <- function(network, K, cluster_intercept = FALSE, verbose = TRUE, 
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


#' Function to perform cross-validation on the CAMM-SBM model. This function can  
#' be used to select the number of clusters or to perform feature selection
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param K List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param folds Number of folds of cross-validation (default = 5)
#' @param mcmc Boolean indicating whether to run the model in MCMC format.
#' The model will run much slower, but may be more accurate
#' @return Combined BIC of all the test sets
#cammsbm_cv <- function(network, K, folds = 5, mcmc = FALSE){
#  n <- nrow(network$X)
#  z_est <- numeric(n)
#  cv_sample <- sample(1:n)
#  ll_total <- 0
#  for (i in 1:folds){
#    ## Partition of the network
#    min_ind <- round((i - 1) * n / folds) + 1
#    max_ind <- round(i * n / folds)
#    train_inds <- cv_sample[-(min_ind:max_ind)]
#    test_inds <- cv_sample[min_ind:max_ind]
#    train_network <- list(A = network$A[train_inds, train_inds], 
#                          dis = network$dis[train_inds, train_inds], 
#                          X = network$X[train_inds, ])
#    test_network <- list(A = network$A[test_inds, test_inds], 
#                         dis = network$dis[test_inds, test_inds], 
#                         X = network$X[test_inds, ])
#    if (!mcmc){
#      fit <- cammsbm_vem(train_network, K, verbose = FALSE)
#      #fit <- cammsbm_plem(train_network, K, verbose = FALSE, 
#      #                  cluster_intercept = TRUE)
#    } else {
#      fit <- cammsbm_nimble(train_network, K)
#    }
#    full_gamma <- matrix(0, n, K)
#    full_gamma[train_inds, ] <- fit$gamma
#    ## Estimate gamma across test nodes
#    for (i in test_inds) {
#      log_prob <- rep(0, K)
#      for (k in 1:K) {
#        acc <- 0
#        for (j in train_inds) {
#          for (l in 1:K) {
#            eta <- 1 / (1 + 
#                  exp(-(fit$beta0[k, l] + fit$beta[k, l] * network$dis[i, j])))
#            acc <- acc + full_gamma[j, l] * (
#              network$A[i, j] * log(eta + 1e-100) + 
#              (1 - network$A[i, j]) * log(1 - eta + 1e-100)
#            )
#          }
#        }
#        log_prob[k] <- acc
#      }
#      log_prob <- log_prob - max(log_prob)
#      full_gamma[i, ] <- exp(log_prob) / sum(exp(log_prob))
#    }
#    ## Edge prediction
#    pred_test_test <- matrix(NA, n, n)
#    for (i in test_inds) {
#      for (j in test_inds) {
#        if (i < j) {
#          eta_kl <- outer(1:K, 1:K, 
#              Vectorize(function(k, l) 
#                fit$beta0[k, l] + fit$beta[k, l] * network$dis[i, j]))
#          prob <- sum(full_gamma[i, ] %*% (1 / (1 + exp(-eta_kl))) * full_gamma[j, ])
#          pred_test_test[i, j] <- pred_test_test[j, i] <- prob
#        }
#      }
#    }
#    ## Evaluate log-likelihood on test network
#    true_tt <- test_network$A
#    pred_tt <- pred_test_test[test_inds, test_inds]
#    mask_tt <- upper.tri(true_tt) & !is.na(pred_tt)
#    ll_pred <- sum(true_tt[mask_tt] * log(pred_tt[mask_tt] + 1e-100) +
#                  (1 - true_tt[mask_tt]) * log(1 - pred_tt[mask_tt] + 1e-100))
#    ll_total <- ll_total + ll_pred
#  }
#  ## Return sum of BICs on out-of-sample networks
#  return(-2 * ll_total + (K * (K + 1)) * log(n * (n - 1) / 2))
#}


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


#' Fixed-point variational EM iteration to be used for squarem
#' @param beta_gamma_vec Vector containing initial beta parameter estimates
#' and gamma estimates, flattened into a single vector
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param K Number of clusters
#' 
#' @return Vector of updated betas and gammas
#' @export
cammsbm_q_vem <- function(beta_gamma_vec, network, K){
  interc <- (K^2 + K) / 2
  n <- nrow(network$A)
  beta0 <- vec2mat(beta_gamma_vec[1:interc])
  beta <- vec2mat(beta_gamma_vec[interc + 1:interc])
  gamma <- matrix(beta_gamma_vec[-(1:(2 * interc))], n, K)
  ## E-STEP
  gamma <- vem_e(network, gamma, beta0, beta)
  ## M-STEP
  beta_new <- vem_m(network, gamma, beta0, beta)
  return(c(mat2vec(beta_new$beta0), mat2vec(beta_new$beta), gamma))
}


#' Variational EM implementation for CAMM-SBM
#' @param network List object containing adjacency (A), similarity (dis) and 
#' covariates (X)
#' @param K Number of clusters to include
#' @param verbose Boolean for verbosity
#' @param tol ELBO gain threshold (default = 0.001). Less than this
#' gain will stop the algorithm)
#' @param gamma Initial values of gamma (n x K matrix) for a warm start
#' @param beta0 Initial values of beta0 (K x K matrix) for a warm start
#' @param beta Initial values of beta (K x K matrix) for a warm start
#' @return Clustering probabilities, beta parameters, along with 
#' diagnostics like iterations, ELBO values, and computing time
#' @export
cammsbm_vem <- function(network, K, verbose = TRUE, tol = 1e-03,
                        gamma = NULL, beta0 = NULL, beta = NULL){
  start_time <- Sys.time()
  initial_elbo <- -Inf; gain <- Inf
  elbo_values <- c()
  ## INITIALIZE GAMMA AND BETA USING MIXTURE NORMAL
  if (is.null(gamma)){
    gamma <- mclust::Mclust(network$X, K, verbose = FALSE)$z
  }
  if (is.null(beta0)){
    initial_beta <- vem_m(network, gamma)
    beta0 <- initial_beta$beta0; beta <- initial_beta$beta
  }
  if(verbose){
    message('Initial Parameters Set...')
  }
  iter <- 1
  ## BEGIN EM ALGORITHM
  while (gain > tol) {
    ## E-STEP: NODE-LEVEL PROBABILITIES
    gamma <- vem_e(network, gamma, beta0, beta)
    ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
    beta_new <- vem_m(network, gamma, beta0, beta)
    beta0 <- beta_new$beta0; beta <- beta_new$beta
    se_beta0 <- beta_new$se_beta0; se_beta <- beta_new$se_beta
    ## Update expected log-likelihood
    new_elbo <- find_elbo(c(mat2vec(beta0), mat2vec(beta)), network, gamma, 
                          cluster_intercept = TRUE)
    gain <- new_elbo - initial_elbo
    initial_elbo <- new_elbo
    elbo_values <- c(elbo_values, new_elbo)
    if (verbose){
      cat(sprintf("\r[EM] Iteration %3d | ELBO: %10.4f", iter, new_elbo))
      flush.console()
    }
    iter <- iter + 1
  }
  cat("\n")
  end_time <- Sys.time()
  result <- list(gamma = gamma,
                 beta0 = beta0, 
                 beta = beta,
                 se_beta0 = se_beta0,
                 se_beta = se_beta,
                 elbo = new_elbo,
                 time = as.numeric(difftime(end_time, start_time, units = 's')),
                 iter = length(elbo_values),
                 evals = elbo_values
  )
  class(result) <- 'cammsbm'
  return(result)
}


#' Function to find the ELBO for the CAMM-SBM model.
#' @param beta_vec beta_0 an beta parameters flattened into 
#' a combined vector
#' @param network List object containing adjacency (A), similarity (dis) and 
#' covariates (X)
#' @param gamma Matrix of clustering probabilities found in the E-Step
#' @param cluster_intercept Whether to have the intercept
#' term represented as a matrix (default = \code{TRUE})
#' @return Observed ELBO with the given parameters
#' @export
find_elbo <- function(beta_vec, network, gamma, cluster_intercept = TRUE){
  K <- ncol(gamma); n <- nrow(gamma)
  interc <- 1:(length(beta_vec) / 2)
  beta0 <- vec2mat(beta_vec[interc])
  beta <- vec2mat(beta_vec[-interc])
  elbo <- 0
  for (k in 1:K) {
    for (l in k:K) {
      eta_mat <- 1 / (1 + exp(-(beta0[k, l] + beta[k, l] * network$dis)))
      log_p <- network$A * log(eta_mat + 1e-100) + 
        (1 - network$A) * log(1 - eta_mat + 1e-100)
      ## Outer product of gamma[, k] and gamma[, l]
      gamma_kl <- gamma[, k] %*% t(gamma[, l])
      contrib <- gamma_kl * log_p
      if (k == l) {
        elbo <- elbo + sum(contrib[upper.tri(contrib)])
      } else {
        elbo <- elbo + sum(contrib[upper.tri(contrib)] + t(contrib)[upper.tri(contrib)])
      }
    }
  }
  elbo <- elbo - sum(gamma * log(gamma + 1e-100))
  return(elbo)
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
#' @note Internal helper for \code{cammsbm_plem_exact} and \code{calfsbm_q_exact}
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
#' @export
#' @note Small helper function
mat2vec <- function(beta, directed = FALSE){
  if (directed){ return(as.vector(beta[upper.tri(beta, diag = TRUE)])) }
  else { return(beta[upper.tri(beta, diag = TRUE)]) }
}


#' Function to run the SEM Algorithm to estimate variance of parameters
#' 
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
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
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
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


#' Helper function to switch from a vector to a matrix representation of beta
#' @param beta Vector of entries
#' @param directed Logical indicating whether network is directed
#'
#' @return Matrix of entries
#' @export
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


#' Helper E-Step for VEM algorithm
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param gamma n x K matrix of clustering probabilities
#' @param beta0 K x K intercept matrix
#' @param beta K x K matrix
#' @return Fitted n x K matrix of clustering probabilities
#' @export
vem_e <- function(network, gamma, beta0, beta) {
  n <- nrow(network$A)
  K <- ncol(gamma)
  log_gamma_new <- matrix(0, n, K)
  
  for (k in 1:K) {
    log_prob_k <- matrix(0, n, n)
    
    for (l in 1:K) {
      ## Compute eta for all i,j at once
      eta_mat <- 1 / (1 + exp(-(beta0[k, l] + beta[k, l] * network$dis)))
      
      ## Compute log-likelihood matrix
      log_eta <- log(eta_mat + 1e-100)
      log_1_eta <- log(1 - eta_mat + 1e-100)
      log_p_mat <- network$A * log_eta + (1 - network$A) * log_1_eta
      
      diag(log_p_mat) <- 0 # No self-loops
      
      ## Weight by gamma[k, l]
      weighted_log_p <- t(t(log_p_mat) * gamma[, l])
      
      ## Sum over all j
      log_prob_k <- log_prob_k + weighted_log_p
    }
    
    log_gamma_new[, k] <- rowSums(log_prob_k)  # Sum over j
  }
  
  ## Stability check and probabilities
  log_gamma_new <- log_gamma_new - apply(log_gamma_new, 1, max)
  gamma_new <- exp(log_gamma_new)
  gamma_new <- gamma_new / rowSums(gamma_new)
  
  return(gamma_new)
}


#' Helper M-Step for VEM algorithm
#' @param network List object containing adjacency (A), similarity (dis), 
#' and covariates (X)
#' @param gamma n x K matrix of clustering probabilities
#' @param beta0 Optional K x K intercept matrix (default = \code{NULL})
#' @param beta K x K matrix (default = \code{NULL})
#' @return List of fitted beta0 and beta matrices
#' @export
vem_m <- function(network, gamma, beta0 = NULL, beta = NULL) {
  n <- nrow(network$A)
  K <- ncol(gamma)
  ## Get upper triangular indices (i < j)
  idx <- which(upper.tri(network$A), arr.ind = TRUE)
  i_vec <- idx[, 1]
  j_vec <- idx[, 2]
  d_vec <- network$dis[idx]
  A_vec <- network$A[idx]
  
  ## Initialize betas if warm start not provided
  if (is.null(beta0)){
    beta0 <- matrix(0, K, K)
    beta <- matrix(0, K, K)
  }
  ## Initialize standard errors
  se_beta0 <- se_beta <- matrix(0, K, K)
  ## Iterate over all valid cluster pairs
  for (k in 1:K) {
    for (l in k:K) {
      gamma_i_k <- gamma[i_vec, k]
      gamma_j_l <- gamma[j_vec, l]
      if (k == l) {
        weights <- gamma_i_k * gamma_j_l
      } else {
        gamma_i_l <- gamma[i_vec, l]
        gamma_j_k <- gamma[j_vec, k]
        weights <- gamma_i_k * gamma_j_l + gamma_i_l * gamma_j_k
      }
      df <- data.frame(y = A_vec, x = d_vec, w = weights)
      ## Logistic fit with warm start
      fit <- suppressWarnings(glm(y ~ x, weights = w, family = binomial(), 
                                  data = df, start = c(beta0[k, l], beta[k, l])))
      beta0[k, l] <- coef(fit)[1]
      beta[k, l] <- coef(fit)[2]
      se_kl <- sqrt(diag(vcov(fit)))
      se_beta0[k, l] <- se_kl[1]
      se_beta[k, l] <- se_kl[2]
      
      if (k != l) {
        beta0[l, k] <- beta0[k, l]
        beta[l, k] <- beta[k, l]
        se_beta0[l, k] <- se_beta0[k, l]
        se_beta[l, k] <- se_beta[k, l]
      }
    }
  }
  
  return(list(beta0 = beta0, beta = beta, se_beta0 = se_beta0, se_beta = se_beta))
}


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
#' in the \code{cammsbm_plem} function, deriving an estimate for beta using 
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
#' in the \code{cammsbm_plem_exact} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
update_beta_exact <- function(K, group, init = NULL, directed = FALSE, 
                              offset = FALSE, cluster_intercept = TRUE){
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

