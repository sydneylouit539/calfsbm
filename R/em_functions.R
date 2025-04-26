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
#' 
#' @return AIC of cross-validated model
calfsbm_cv <- function(network, K, folds = 4){
  n <- nrow(network$X)
  z_est <- numeric(n)
  set.seed(1)
  #  cv_sample <- sample(1:n)
  cv_sample <- 1:n
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
    fit <- calfsbm_em(train_network, K, verbose = FALSE)
    print(fit$z); print(fit$beta)
    #print(fit$aic)
    test_z <- c()
    for (i in test_inds){
      z_prob <- numeric(K)
      for (j in 1:K){
        eta_ij <- fit$beta0 + fit$beta[j, fit$z] * 
          network$dis[i, train_inds]
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
    etas <- fit$beta0 + fit$beta[m[, 1] * (K - 1) + m[, 2]] * 
      network$dis[ti[, 1] * (n - 1) + ti[, 2]]
    p_hat <- 1 / (1 + exp(-etas))
    ll_pred <- sum(log(ifelse(network$A[ti[, 1] * (n - 1) + ti[, 2]] == 1, p_hat, 1 - p_hat)))
    ll_total <- ll_total + ll_pred
  }
  #print(z_est)
  return(-2 * ll_total + (K * (K + 1) + 2))
}


#' EM algorithm implementation to get variational Bayesian distribution
#' 
#' @param network A list object containing adjacency matrix A, 
#' @param K Number of clusters
#' @param offset Boolean to indicate if there should be an offset (default = TRUE)
#' @param verbose Verbosity (default = TRUE)
#' @param conf Confidence level of betas (default = 0.95)
#' @param S Number of conditional maximization iterations per M-step (default = 1)
#' 
#' @return List of estimated node membership, betas, and AIC
#' @export
calfsbm_em <- function(network, 
                       K, 
                       offset = FALSE, 
                       cluster_intercept = FALSE,
                       verbose = TRUE, 
                       conf = 0.95,
                       S = 1){
  start_time <- Sys.time()
  initial_ll <- -Inf; gain <- Inf
  ## SET UP PARAMETERS
  n <- nrow(network$A)
  log_prob_mat <- matrix(0, nrow = n, ncol = K)
  ## INITIALIZE Z AND BETA
  z <- mclust::Mclust(network$X, K, verbose = FALSE)$z
  if (offset){
    initial_theta <- update_theta(network$A, network$dis, 
                                  apply(z, 1, which.max), 0, matrix(0, K, K), numeric(n))
    group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                        S_ij = network$dis, offset = initial_theta)
  } else {
    group <- gen_factor(apply(z, 1, which.max), A = network$A, 
                        S_ij = network$dis, offset = NULL)
  }
  #initial_beta <- update_beta_em(K, group$cluster, offset = offset)
  initial_beta <- update_beta_em(K, group$cluster, cluster_intercept = cluster_intercept)
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
    print(new_ll)
    gain <- new_ll - initial_ll
    initial_ll <- new_ll
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
      #beta_new <- update_beta_em(K, group$cluster, offset = offset, conf = conf)
      beta_new <- update_beta_em(K, group$cluster, 
                                 cluster_intercept = cluster_intercept, conf = conf)
    }
    beta0 <- beta_new$beta0
    beta <- beta_new$beta
    #new_aic <- beta_new$aic
    #gain <- initial_aic - new_aic
    #initial_aic <- new_aic
    if(verbose){
      print('M STEP COMPLETED!')
      #print(paste('Current AIC:', round(new_aic, 1)))
    }
  }  
  end_time <- Sys.time()
  return(list(z = apply(z, 1, which.max),
              z_soft = z,
              beta0 = beta0, 
              beta = beta,
              low = beta_new$low,
              high = beta_new$high,
              theta = new_theta,
              loglik = new_ll,
              model = beta_new$model,
              time = as.numeric(difftime(end_time, start_time, units = 's'))))
}


#' Function to perform an exact EM algorithm on network data
#' @param network List object containing similarity matrix dis, 
#' @param K Number of clusters
#' @param offset Boolean indicating whether to use random-effect terms
#' @param verbose Boolean indicating the verbosity of the function
#' 
#' @return Soft clustering estimate z, beta0 and beta estimates, 
#' log-likelihood, GLM output, and computing time
#' @export
calfsbm_em_exact <- function(network, K, offset = FALSE, verbose = TRUE){
  start_time <- Sys.time()
  initial_ll <- -Inf; gain <- Inf
  ## SET UP PARAMETERS
  n <- nrow(network$A)
  log_prob_mat <- matrix(0, nrow = n, ncol = K)
  ## INITIALIZE Z AND BETA
  z <- mclust::Mclust(network$X, K)$z
  group <- gen_factor_soft(z, A = network$A, 
                           S_ij = network$dis, offset = offset)
  initial_beta <- update_beta_exact(K, group$cluster, offset = offset)
  beta0 <- initial_beta$beta0; beta <- initial_beta$beta
  if(verbose){
    print('Parameters Set!')
  }
  ## BEGIN EM ALGORITHM
  while (gain > 0.000001) {
    ## E-STEP: NODE-LEVEL PROBABILITIES
    for (i in 1:n){
      loglik <- matrix(0, K, K)
      for (j in 1:K){
        ## Average over all possible configurations
        for (k in 1:K){
          eta_ij <- beta0 + beta[j, k] * network$dis[i, -i]
          ## Calculate probabilities
          fit <- 1 / (1 + exp(-eta_ij))
          ## Calculate log-likelihood
          loglik[j, k] <- #z[i, j] * 
            sum(log(ifelse(network$A[i, -i] == 1, fit, 1 - fit)) * z[-i, k])
        }
        log_prob_mat[i, j] <- sum(loglik[j, ])
      }
    }
    for (i in 1:n){
      ## Get relative probabilities of each row
      wts <- log_prob_mat[i, ]
      wts <- wts - max(wts) #  Safety check for extremely low likelihood
      z[i, ] <- exp(wts) / sum(exp(wts))
    }
    if(verbose){
      print('E STEP COMPLETED!')
    }
    ## M-STEP: RUN LOGISTIC REGRESSION AND UPDATE BETAS/THETAS
    group <- gen_factor_soft(z, A = network$A, 
                             S_ij = network$dis, offset = offset)
    beta_new <- update_beta_exact(K, group$cluster, 
                                  init = c(beta0, mat2vec(beta)), offset = offset)
    beta0 <- beta_new$beta0
    beta <- beta_new$beta
    ## Update expected log-likelihood
    new_ll <- sum(log_prob_mat * z)
    gain <- new_ll - initial_ll
    initial_ll <- new_ll
    if(verbose){
      print('M STEP COMPLETED!')
      print(paste('Current Log-Likelihood:', round(new_ll, 1)))
    }
  }  
  end_time <- Sys.time()
  return(list(z = z,
              beta0 = beta0, 
              beta = beta,
              low = beta_new$low,
              high = beta_new$high,
              loglik = new_ll,
              time = as.numeric(difftime(end_time, start_time, units = 's')),
              model = beta_new$model))
}


#' Fixed point EM iteration to be used for squarem
#' @param beta description
#' @param network List object containing the following attributes:
#' dis: Similarity matrix
#' A: Adjacency matrix
#' @param z Clustering vector
#' @return Vector of updated betas
calfsbm_q <- function(beta, network, z){
  ## E-STEP: NODE-LEVEL PROBABILITIES
  beta0 <- beta[1]
  beta <- vec2mat(beta[-1])
  initial_z <- apply(z, 1, which.max)
  updated_z <- initial_z
  log_prob_mat <- matrix(0, n, K)
  for (i in 1:n){
    for (j in 1:K){
      eta_ij <- beta0 + 
        beta[j, updated_z[-i]] * network$dis[i, -i]
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
  #beta_new <- update_beta_em(K, group$cluster, offset = offset, conf = conf)
  beta_new <- update_beta_em(K, group$cluster)
  return(c(beta_new$beta0, mat2vec(beta_new$beta)))
}


#' Exact fixed-point EM iteration to be used for squarem
#' @param beta Vector containing initial beta parameter estimates
#' @param network List object containing the following attributes:
#' dis: Similarity matrix
#' A: Adjacency matrix
#' @param z Soft clustering matrix, with each row corresponding to a node,
#' and each column corresponding to a cluster
#' 
#' @return Vector of updated betas
calfsbm_q_exact <- function(beta, network, z){
  ## E-STEP: NODE-LEVEL PROBABILITIES
  beta0 <- beta[1]
  beta <- vec2mat(beta[-1])
  n <- nrow(z); K <- nrow(beta)
  log_prob_mat <- matrix(0, n, K)
  for (i in 1:n){
    loglik <- matrix(0, K, K)
    for (j in 1:K){
      ## Average over all possible configurations
      for (k in 1:K){
        eta_ij <- beta0 + beta[j, k] * network$dis[i, -i]
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
  beta_new <- update_beta_exact(K, group$cluster, init = c(beta0, mat2vec(beta)), 
                                offset = offset)
  beta0 <- beta_new$beta0
  beta <- beta_new$beta
  #print(beta0); print(mat2vec(beta))
  return(c(beta_new$beta0, mat2vec(beta_new$beta)))
}


#' Function to estimate covariance matrix of EM parameters using the method
#' described in Louis (1982)
#' @param z Matrix of clustering probabilities
#' @param beta0 Intercept
#' @param beta Matrix of beta parameters for each pair of clusters
#' @param group Matrix used in logistic regression fit
em_covar <- function(z, S, A, group){
  z_est <- apply(z, 1, which.max)
  mod_mat <- stats::model.matrix(~ 0 + as.factor(cl), group$cluster) * group$cluster$x
  mod_mat2 <- as.data.frame(mod_mat)
  mod_mat2$y <- group$cluster$y
  logit_fit <- glm(y ~ ., family = 'binomial', data = mod_mat2)
  beta0 <- logit_fit$coefficients[1]
  beta <- vec2mat(logit_fit$coefficients[-1])
  ## Observed Fisher Information
  fisher_info <- unname(solve(vcov(logit_fit)))
  ## Constants
  n <- nrow(z); K <- ncol(z)
  n_param <- 1 + (K^2 + K) / 2
  yhat <- unname(logit_fit$fitted.values)
  ## Observed score
  score <- matrix(0, nrow(mod_mat), n_param)
  mod_mat <- cbind(1, mod_mat)
  for (col in 1:n_param){
    score[, col] <- (mod_mat2$y - yhat) * mod_mat[, col]
  }
  ## Expected score
  e_score <- matrix(0, nrow(mod_mat), n_param)
  e_score[, 1] <- score[, 1] # Intercept (universal, no expectation needed)
  
  cl <- which(upper.tri(matrix(1:K^2, K, K), TRUE))
  for (col in 2:n_param){ # Probability of each node multiplied by score
    e_score[, col] <- 
      z[group$grid[, 1], ((cl[col - 1] - 1) %% K) + 1] * 
      z[group$grid[, 2], ((cl[col - 1] - 1) %/% K) + 1] *
      (mod_mat2$y - yhat) * group$cluster$x
    print(c(cl[col - 1], ((cl[col - 1] - 1) %% K) + 1, ((cl[col - 1] - 1) %/% K) + 1))
  }
  ## Covariance of score (S - E[S])'(S - E[S])
  cov_score <- t(score - e_score) %*% (score - e_score)
  #print(cov_score)
  print(solve(fisher_info + cov_score))
  return(list(score = score, 
              e_score = e_score,
              fisher_info = fisher_info,
              cov_score = cov_score,
              cov = solve(fisher_info + cov_score)))
}


#' Function to perform forward stepwise variable selection
#' @param network List object containing adjacency, similarity, and covariates
#' @param K Number of clusters
#' @param cv Whether to use cross-validation (default = FALSE)
#' 
#' @return Vector of covariate indices which correspond to the optimal model
forward_stepwise_em <- function(network, K, cv = FALSE) {
  current_aic <- Inf
  remaining_predictors <- 1:ncol(network$X)  # List of all predictor names
  selected_predictors <- c(1)          # List to store selected predictors
  while(length(remaining_predictors) > 0) {
    aic_values <- numeric(length(remaining_predictors))
    print(remaining_predictors)
    for (i in 1:length(remaining_predictors)) {
      predictor <- remaining_predictors[i]
      S_ij <- dist(network$X[, c(selected_predictors[-1], predictor)], upper = TRUE, diag = TRUE)
      S_ij <- as.matrix(Matrix::forceSymmetric(as.matrix(S_ij)))
      diag(S_ij) <- 0
      network$dis <- S_ij
      if (cv){
        aic_values[i] <- calfsbm_cv(network, K)
        #+ 2 * length(selected_predictors)
      } else {
        aic_values[i] <- calfsbm_em(network, K, verbose = FALSE)$aic + 
          2 * length(selected_predictors)
      }
    }
    print(aic_values)
    best_index <- which.min(aic_values)
    if (aic_values[best_index] < current_aic) {
      selected_predictors <- c(selected_predictors, remaining_predictors[best_index])
      remaining_predictors <- remaining_predictors[-best_index]
      current_aic <- aic_values[best_index]  # Update the current AIC
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
#' @param directed Boolean indicating whether the network should be directed
#' @param offset Boolean indicating whether to include node-level 
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
#' @param directed Boolean indicating whether network is directed
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
#' @param beta Fitted beta values from the EM algorithm
#' @param z Fitted matrix of z probabilities
#' @param obs_fisher Observed Fisher information from the EM fit
#' @param increment Small increment to calculate numerical Jacobian 
#' @return Covariance matrix of betas
#' @note Internal helper
sem_richardson <- function(network, K, beta, z, obs_fisher, increment = 0.0001){
  m <- length(beta)
  DM <- matrix(0, m, m)
  ## ITERATE OVER BETA_T VECTOR
  for (j in 1:m){
    u_j <- ifelse(1:m == j, increment, 0)
    ## FULL ITERATION OF EM FROM MODIFIED BETA VALUES
    beta_h <- calfsbm_q_exact(beta + u_j, network, z)
    beta_2h <- calfsbm_q_exact(beta + 2 * u_j, network, z)
    beta_nh <- calfsbm_q_exact(beta - u_j, network, z)
    beta_n2h <- calfsbm_q_exact(beta - 2 * u_j, network, z)
    DM[j, ] <- (beta_n2h - 8 * beta_nh + 8 * beta_h - beta_2h) / 
      (12 * increment)
  }
  V <- solve(obs_fisher) %*% solve(diag(m) - DM)
  return(V)
}


#' Dr. Yan's stochastic derivative code
#' @param g Fixed-point function
#' @param x Parameters
#' @param n Noise
#' @param B Number of iterations
#' @param s Noise parameter
#' @param network Network
#' @param z_soft Soft clustering matrix
#' 
#' @return Sampled Jacobian matrix
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
#' @param directed Boolean indicating whether network is directed
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
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param offset Boolean indicating whether to use offset term
#' @param conf Confidence level of beta estimates (default = 0.95)
#' 
#' @return beta0 and beta, with beta as a matrix
#' 
#' @note Function \code{update_beta_em} is a helper in the 
#' in the \code{calfsbm_em} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
update_beta_em <- function(K, group, init = NULL, directed = FALSE, 
                           cluster_intercept = FALSE, conf = 0.95){
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
    logit_fit <- glm(y ~ 0 + . , family = 'binomial', data = mod_mat2)
  }
  sampled_beta <- logit_fit$coefficients
  ## Inverse Fisher Information for uncertainty quantification
  vars <- diag(vcov(logit_fit))
  fac <- -qnorm((1 - conf) / 2)
  ## Intercept term indexing
  if (cluster_intercept) {interc <- 1:(K * (K + 1) / 2); } 
  else {interc <- 1}
  
  beta0_low <- sampled_beta[interc] - fac * sqrt(vars[interc])
  beta0_high <- sampled_beta[interc] + fac * sqrt(vars[interc])
  ## Convert beta_kl to matrix form
  if (!directed) {
    beta_mat <- vec2mat(sampled_beta[-interc])
    low <- vec2mat(sampled_beta[-interc] - fac * sqrt(vars[-interc]))
    high <- vec2mat(sampled_beta[-interc] + fac * sqrt(vars[-interc]))
  } else {
    beta_mat <- matrix(sampled_beta[-interc], K, K)
    low <- matrix(sampled_beta[-interc] - fac * vars[-interc], K, K)
    high <- matrix(sampled_beta[-interc] + fac * vars[-interc], K, K)
  }
  if(cluster_intercept){
    intercept <- vec2mat(sampled_beta[interc])
    beta0_low <- vec2mat(sampled_beta[interc] - fac * vars[-interc])
    beta0_high <- vec2mat(sampled_beta[interc] + fac * vars[-interc])
  } 
  else {intercept <- sampled_beta[1]}
  return (list(beta0 = intercept, 
               beta = beta_mat, 
               aic = logit_fit$aic,
               low = list(beta0 = beta0_low, beta = low),
               high = list(beta0 = beta0_high, beta = high),
               model = logit_fit))
}


#' Update Beta in exact EM form
#' 
#' Helper function to update beta according to adjacency and node membership
#' @param K A positive integer indicating the true number of clusters
#' @param group Model matrix to be fitted on
#' @param init Starting values for logistic regression
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @param offset Boolean indicating whether to use offset term
#' @param conf Confidence level of beta estimates (default = 0.95)
#' 
#' @return beta0 and beta, with beta as a matrix
#' 
#' @note Function \code{update_beta_exact} is a helper in the 
#' in the \code{calfsbm_em_exact} function, deriving an estimate for beta using 
#' the initial clustering configuration as input
update_beta_exact <- function(K, group, init = NULL, directed = FALSE, 
                              offset = FALSE, conf = 0.95){
  mod_mat <- stats::model.matrix(~ 0 + as.factor(cl), group) * group$x 
  mod_mat2 <- as.data.frame(mod_mat)
  mod_mat2$y <- group$y
  mod_mat2$wts <- group$wts
  ## Maximize likelihood through logistic regression
  logit_fit <- suppressWarnings(stats::glm(y ~ . - wts, family = 'binomial', 
                                           data = mod_mat2, weights = wts, start = init))
  sampled_beta <- logit_fit$coefficients
  ## Inverse Fisher Information for uncertainty quantification
  vars <- diag(vcov(logit_fit))
  fac <- -qnorm((1 - conf) / 2)
  ## Convert beta_kl to matrix form
  beta_mat <- low <- high <- matrix(0, K, K)
  if (!directed) {
    beta_mat <- vec2mat(sampled_beta[-1])
    low <- vec2mat(sampled_beta[-1] - fac * sqrt(vars[-1]))
    high <- vec2mat(sampled_beta[-1] + fac * sqrt(vars[-1]))
  } else {
    beta_mat <- matrix(sampled_beta[-1], K, K)
    low <- matrix(sampled_beta[-1] - fac * sqrt(vars[-1]), K, K)
    high <- matrix(sampled_beta[-1] + fac * sqrt(vars[-1]), K, K)
  }
  return (list(beta0 = sampled_beta[1], 
               beta = beta_mat, 
               aic = logit_fit$aic,
               low = low,
               high = high,
               model = logit_fit))
}

