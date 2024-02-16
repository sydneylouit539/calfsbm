############## CODE CONTAINING FUNCTIONS FOR OUR SBM ###########################
library(bcdc)
library(boa)
library(igraph)
library(label.switching)
library(latentnet)
library(MASS)
library(mclust)
library(nett)
library(network)
library(nimble)
library(raster)

#' Visualize Link Probability vs. Actual Matrix
#' 
#' Using the adjacency matrix, returns a plot of the adjacency matrix, and 
#' a raster of the true probabilities
#' @param observed Observed adjacency matrix
#' @param prob Matrix of intrinsic probabilities
#' @param z_tru Vector of the true node membership
#' @return Graph of the connectivity and graph of true probabilities
#' @examples links <- gen_az(n_nodes = 50, K = 2, m = 2, prob = c(0.5, 0.5),
#'beta0 = 1, beta = diag(2) - 3, sigma = 0.3, spat = 0.5)
#' plot_conn(links$A, links$lp, links$z)
#' @export
plot_conn <- function(observed, prob, z_tru){
    n <- length(z_tru)
    par(mar = c(1, 1, 1, 1))
    par(mfrow = c(1, 2))
    B <- observed[order(z_tru), order(z_tru)]
    xx <- raster::rasterFromXYZ(cbind(expand.grid(1:n, n:1), c(B)))
    plot(xx, legend=FALSE, axes=FALSE)
    ## Plot true probabilities
    B <- prob[order(z_tru), order(z_tru)]
    xx <- raster::rasterFromXYZ(cbind(expand.grid(1:n, n:1), c(B)))
    plot(xx, axes=FALSE)
}


#' Helper function for CALF-SBM Initialization
#' 
#' Get the X matrix for logistic regression by taking clusters, links, and dist.
#' @param initial_z Initial node membership vector
#' @param A Known adjacency matrix
#' @param S_ij Distance matrix
#' @param directed Boolean 
#' @param offset Logistic regression offset terms (currently not used)
#' @return List of X matrix and grid of (x, y) indices
#' @note
#' Internal helper function
#' @export
gen_factor <- function(initial_z, A, S_ij, directed = FALSE, offset = FALSE){
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
    if (offset){
        ## Offset terms theta_i and theta_j to be used in logistic model
        cluster$offset <- degrees[eg[wuta, 1]] + degrees[eg[wuta, 2]]
    }
    return (list(cluster = cluster, grid = eg[wuta,]))
}


#' Helper Function to Make Adjacency Matrix Symmetric
#' 
#' Makes a non-symmetric matrix symmetric using the lower diagonal
#' @param mat The matrix to be symmetrized
#' @return Symmetric matrix with the same upper diagonal as the original
#' @note Internal helper function
#' @export
makesymmetric <- function(mat){
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    return (mat)
}


#' (Obsolete) Helper Function to Generate Covariates
#' 
#' Generate matrix of covariates X according to parameters
#' @param x_dim Dimension of X
#' @param n_nodes Number of nodes
#' @param sigma correlation (or correlation matrix if input is xdim by xdim)
#' @return Matrix of covariates
#' @note Deprecated
#' @export
gen_x <- function(x_dim, n_nodes, sigma){
    ## If sigma is given as a matrix
    if (length(sigma) == x_dim^2){
        X <- mvrnorm(n_nodes, mu = rep(0, x_dim), Sigma = sigma)
    } else if (length(sigma) == 1){
        X <- mvrnorm(n_nodes, mu = rep(0, x_dim), 
                     Sigma = (1 - sigma) * diag(x_dim) + sigma)
    } else {
        X <- mvrnorm(n_nodes, mu = rep(0, x_dim), Sigma = diag(x_dim))
    }
    return(X)
    z <- rcat(n_nodes, rep(K^-1, K))
    X <- matrix(0, nrow = n_nodes, ncol = x_dim)
    ## Centers of new clusters
    areas <- mvrnorm(K, mu = rep(0, x_dim), 
                     Sigma = spat * (1 - sigma) * diag(x_dim) + sigma)
    for (i in 1:K){
      zz <- which(z == i)
      X[zz, ] <- sweep(mvrnorm(length(zz), mu = rep(0, x_dim), 
                               Sigma = diag(x_dim)), 2, areas[i, ], '+')
    }
}


#' CALF-SBM Network Generation Function
#'
#' Generates network using probabilities, betas, and signal strength
#' @param n_nodes Number of nodes in the network
#' @param K A positive integer indicating the true number of clusters 
#' @param n_covar A positive integer indicating the number of node-specific
#' covariates to generate.
#' @param prob Vector of weights of being in each group. This should be
#' a numeric vector of length K, and can be weights or probabilities
#' @param beta0 Intercept term, lower values correspond to higher sparsity
#' @param beta K-by-K matrix of within and between-cluster coefficients. 
#' Diagonal values correspond to within-cluster effects
#' @param sigma Non-negative number for the standard deviation of the random 
#' effects offset term, set to 0 if not using offset
#' @param spat Non-negative number to represent the ratio of between-cluster 
#' variance to within-cluster variance of the covariates. 0 indicates no signal
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from 
#' an undirected network
#' @return List of adjacency, membership, link probability, and distance
#' @examples generate_calfsbm_network(100, 2, 2, c(0.5, 0.5), 1, diag(2) - 3,
#'  0.3, 1.5)
#' @export
generate_calfsbm_network <- function(n_nodes, K, n_covar, prob, beta0, beta, 
                                     sigma, spat, directed = FALSE){
    z_tru <- sample(1:K, n_nodes, replace = TRUE, prob = prob)
    ## Spatial correlation
    X <- matrix(0, nrow = n_nodes, ncol = n_covar)
    ## Centers of new clusters
    angles <- 2 * pi * (1:K) / K
    initial_mids <- cbind(sin(angles), cos(angles), -sin(angles), -cos(angles))
    ## Initialize midpoints to have variance 1
    mids <- sqrt(spat * 2) * initial_mids[, 1 + 1:n_covar %% K]
    for (i in 1:K){
      cli <- which(z_tru == i)
      X[cli, ] <- sweep(mvrnorm(length(cli), mu = rep(0, n_covar), 
                Sigma = diag(n_covar)), 2, mids[i, ], '+')
    }
    ## Degree
    theta <- rnorm(n_nodes, 0, sigma)
    ## Generate probabilities, and link probability matrix
    ## Using Euclidean distance, calculate true link probabilities
    ZZ <- dist(X, upper = TRUE, diag = TRUE)
    S_ij <- matrix(0, n_nodes, n_nodes)
    S_ij[upper.tri(S_ij)] <- ZZ; S_ij <- makesymmetric(S_ij); diag(S_ij) <- 0
    eta <- beta0 + outer(theta, theta, '+') + beta[z_tru, z_tru] * S_ij
    diag(eta) <- -Inf
    link_prob_tru <- 1 / (1 + exp(-eta))
    ## Finally, generate adjacency matrix and make symmetric if undirected
    A <- matrix(rbinom(n_nodes^2, 1, link_prob_tru), n_nodes, n_nodes)
    if (!directed){ A <- makesymmetric(A) }
    return (list(A = A, z = z_tru, X = X, lp = link_prob_tru, 
                 dis = S_ij, theta = theta, beta0 = beta0, beta = beta))
}


#' Latentnet Network Generator
#'
#' Function to generate network according to Handcock's RSS latentnet paper 
#' @param n_nodes Number of nodes
#' @param K Number of clusters
#' @param x_dim Number of covariates
#' @param betavec Vector of intercept and slope terms
#' @param noise noise/signal ratio (default=1)
#' @return Network as list object with characteristics
#' @export
generate_latentnet_network <- function(n_nodes, K, x_dim, betavec, noise = 1, 
                                       directed = FALSE){
    z_tru <- sample(1:K, n_nodes, replace = TRUE)
    centroids <- matrix(rnorm(x_dim * K), K, x_dim)
    X <- centroids[z_tru, ] + rnorm(x_dim * n_nodes, sd = sqrt(noise))
    ## Generate probabilities, and link probability matrix
    ## Using Euclidean distance, calculate true link probabilities
    ZZ <- dist(X, upper = TRUE, diag = TRUE)
    S_ij <- matrix(0, n_nodes, n_nodes)
    S_ij[upper.tri(S_ij)] <- ZZ; S_ij <- makesymmetric(S_ij); diag(S_ij) <- 0
    eta <- betavec[1] + betavec[2] * S_ij; diag(eta) <- -Inf
    link_prob_tru <- exp(eta) / (1 + exp(eta))
    ## Finally, generate adjacency matrix and make symmetric if undirected
    A <- matrix(rbinom(n_nodes^2, 1, link_prob_tru), n_nodes, n_nodes)
    if (!directed){ A <- makesymmetric(A) }
    return (list(A = A, z = z_tru, X = X, lp = link_prob_tru, dis = S_ij))
}


#' Latentnet Cluster Finder
#' 
#' Function to find the optimal number of groups using latentnet package
#' @param p Number of clusters. The function will iterate from 1 to p
#' @param AA Network object (from network package) containing adjacency and 
#' covariate information
#' @param sample_size Number of non-rejected MCMC samples
#' @param burnin Number of initial MCMC samples to be discarded
#' @param by Number of samples to skip after burn-in period
#' @return A dataframe with the AIC and BIC values for each K value
#' @examples net <- network::network(m <- matrix(rbinom(25, 1, 0.4), 5, 5));
#' find_K_optimal(6, net, 10000, 2000, 5)
#' @export
find_K_optimal <- function(p, AA, sample_size, burnin, by){
    ## Set a dataframe for BIC and number of clusters
    best_estimate <- data.frame(K = 1:p, BIC = rep(Inf, p), AIC = rep(Inf, p))
    n <- AA$gal$n
    z_est <- rep(1, p)
    for (i in 1:p){
        samp.fit <- latentnet::ergmm(AA ~ euclidean(d = 2, G = i), 
                                     control = control.ergmm(
            sample.size = sample_size, burnin = burnin, interval = by), 
            verbose = TRUE)
        fit <- summary(samp.fit)$bic$overall
        if (fit < min(best_estimate$BIC)){
            z_est <- samp.fit$mkl$Z.K
        }
        best_estimate$BIC[i] <- fit
        n_params <- (i^2 + i) / 2
        best_estimate$AIC[i] <- fit + 2 * n_params - n_params * log(n * (n - 1) / 2)
        print(c(i, fit))
    }
    return (list(K = which.min(best_estimate$BIC), data = best_estimate, z = z_est))
}


#' (Obsolete) Function to Simulate K Frequencies
#'
#' Function to expand K optimization into a simulation study
#' @param nsim The total number of simulations
#' @param n_vals Vector of the total node numbers
#' @param K True number of clusters, can also be vector form
#' @param p The maximum number of clusters to test for
#' @param x_dim Number of covariates in X
#' @param beta0 Intercept term
#' @param beta Effect of clusters on link probability
#' @note Deprecated
#' @return Matrix of true K vs. estimated K
sim_study_K_finder <- function(nsim, n_vals, K, p, x_dim, beta0, beta, 
                               sample_size = 2000, burnin = 4000, by = 4){
    derived_k <- rep(0, nsim)
    if (length(n_vals) == 1){n_vals = rep(n_vals, nsim)}
    if (length(K) == 1){K = rep(K, nsim)}
    for (i in 1:nsim){
        links <- gen_az(n_vals[i], K[i], x_dim, 
                        prob = rep(1/K[i], K[i]), beta0, beta[1:K[i], 1:K[i]])
        AA <- network::network(links$A, as.data.frame(links$X))
        sim_i <- find_K_optimal(p, AA, sample_size, burnin, by)
        derived_k[i] <- which.min(sim_i$BIC)
    }
    return(table(data.frame(n_tru = n_vals, k_tru = K, k_est = derived_k)))
}


#' Update Beta in Gibbs Sampler
#' 
#' Helper function to update beta according to adjacency and node membership
#' @param K Number of clusters
#' @param group Model matrix to be fitted on
#' @param directed Boolean indicating whether network is directed
#' @param offset Boolean indicating whether to use offset term
#' @return beta0 and beta, with beta as a matrix
#' @export
update_beta <- function(K, group, directed = FALSE, offset = FALSE){
    mod_mat <- model.matrix(~ 0 + as.factor(cl), group) * group$x 
    mod_mat2 <- as.data.frame(mod_mat)
    mod_mat2$y <- group$y
    ## Fit Bayesian logistic regression to the data
    if (!offset){
      logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
                        prior.mean = 0, prior.scale = 1)
    } else {
      logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
                        prior.mean = 0, prior.scale = 1, offset = group$offset)
    }
    predicted_beta <- arm::sim(logit_fit, n.sims = 2) 
    sampled_beta <- coef(predicted_beta)[1, ]
    #sampled_beta <- coef(logit_fit)
    ## Convert beta_kl to matrix form
    beta_mat <- matrix(0, K, K)
    if (!directed) {
        beta_mat[upper.tri(beta_mat, diag = TRUE)] <- sampled_beta[-1]
        beta_mat <- makesymmetric(beta_mat)
    } else {
        beta_mat <- matrix(sampled_beta[-1], K, K)
    }
    return (list(beta0 = sampled_beta[1], beta = beta_mat, aic = logit_fit$aic))
}


#' (Obsolete) Helper for Gibbs Sampler
#'
#' Update node membership using likelihood
#' @param z Initial values of z
#' @param X Matrix of covariates
#' @param beta0 Draw from posterior distribution
#' @param beta Draw from posterior distribution
#' @param S_ij Similarity or distance matrix
#' @param A Observed adjacency matrix
#' @param K Number of clusters
#' @param directed Boolean indicating whether network is directed
#' @return Updated node membership vector
#' @export
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

#' Gelman diagnostic
#' 
#' Takes matrix and converts it to a useful format
#' @param results A matrix, where each column represents a chain
#' @param as_matrix Boolean indicating whether matrix or array input is used
#' @param n_chains Number of MCMC chains
#' @return Gelman-Rubin diagnostic (ideal value = 1) 
#' @note Function has been superseded by the return_gelman option in 
#' \code{calf_sbm_nimble} function
#' @export 
gelman <- function(results, as_matrix = TRUE, n_chains = 3) {
  n <- nrow(results)
  if (as_matrix) {
    results <- array(results, dim = c(dim(results) / c(n_chains, 1), n_chains))
    W <- rowMeans(apply(results, c(2, 3), var))
    B <- apply(apply(results, c(2, 3), mean), 1, var)
    n <- dim(results)[1]
    V_hat <- ((n - 1) * W + B) / n
    R <- sqrt(V_hat / W)
    return(mean(R))
  } else {
    W <- mean(apply(results, 2, var)) # Within-chain variance
    B <- n * var(colMeans(results)) # Between-chain variance
    n <- nrow(results)
    V_hat <- ((n - 1) * W + B) / n
    R <- sqrt(V_hat / W) # Gelman-Rubin statistic
    return(R)
  }
}


#' (Obsolete) Check MCMC for convergence
#'
#' Checks MCMC by determining if node membership is stable
#' @param z Node membership from last five iterations
#' @param group Output from \code{gen_factor}
#' @param niter Total iterations initially planned
#' @param iter The current iteration of the Gibbs sampler
#' @return Simulated coefficients
#' @note Deprecated
#' @export
convergence_procedure <- function(z, group, niter, i, K){
    mod_mat <- model.matrix(~ 0 + as.factor(cl), group$cluster) * 
        group$cluster$x
    mod_mat2 <- as.data.frame(mod_mat)
    mod_mat2$y <- group$cluster$y
    ## Fit Bayesian logistic regression to the data
    logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
                             prior.mean = 0, prior.scale = 10)
    predicted_beta <- coef(arm::sim(logit_fit, n.sims = niter - i + 1))
    ## Fill in the history for the rest of the iterations
    beta_history <- array(0, dim = c(niter - i + 1, K, K))
    all_z <- matrix(0, niter - i + 1, ncol(z))
    for (j in i:niter){
        beta_mat <- matrix(0, K, K)
        if (!directed){
            beta_mat[upper.tri(beta_mat, diag = TRUE)] <- 
              predicted_beta[j - i + 1, -1]
            beta_mat <- makesymmetric(beta_mat)
        } else {
            beta_mat <- predicted_beta[j - i + 1, -1]
        }
    beta_history[j - i + 1, , ] <- beta_mat
    all_z[j - i + 1, ] <- z[1, ]
    }
    return(list(preds = predicted_beta, ints = predicted_beta[, 1], 
              bics = logit_fit$aic - 2 * K^2 + K^2 * log(n * (n - 1) /2), 
              beta_his = beta_history, z_vals = all_z))
}


#' (Obsolete) Gibbs sampler
#' 
#' Gibbs sampler with fixed number of clusters
#' @param K Number of clusters (should be optimized already)
#' @param alpha Vector of prior probabilities of being in each group
#' @param beta0 Prior intercept for logistic regression
#' @param beta Prior parameters for logistic regression, in matrix form
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @param directed logical (default = FALSE); if FALSE, the network is undirected
#' @return List containing the estimated beta, node membership, and history
#' @note Deprecated
#' @export
gibbs_sampler_fixed_k <- function(K, alpha, beta0, beta, niter, A, S_ij, 
                                  directed = FALSE){
    n <- nrow(A)
    converged <- FALSE
    ## Initialize beta
    initial_beta <- beta
    ## Initialize node membership
#    initial_z <- sample(1:K, n, replace = TRUE, prob = alpha)
    initial_z <- cluster::pam(links$X, K)$clustering
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


#' (Obsolete) Unknown K Gibbs sampler
#'
#' Gibbs sampler with CRP as initial state
#' @param dp Positive number for CRP concentration parameter (recommended 1)
#' @param X Matrix of node-specific covariates, each row corresponding to a node
#' @param niter Total number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @return List containing the estimated beta, node membership, and history
#' @note Deprecated for now, but may be used in future
#' @export
gibbs_sampler_unknown_k <- function(dp, niter, A, S_ij, directed){
    n <- nrow(A)
    converged <- FALSE
    ## Initialize node membership using Chinese Restaurant Process
    initial_z <- rCRP(n = 1, conc = dp, n)
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
#' @param p Number of clusters to test (from 1 to p)
#' @param alpha Vector of prior probabilities of being in each group
#' @param beta0 Prior intercept for logistic regression
#' @param beta Prior parameters for logistic regression, in matrix form
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @return List containing the estimated beta, node membership, and history
#' @note Deprecated
#' @export
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
#' @return List of beta, z, and history of K
#' @examples \code{links <- generate_calfsbm_network(n_nodes = 50, K = 2, m = 2, 
#' prob = c(0.5, 0.5), beta0 = 1, beta = diag(2) - 3, sigma = 0.3, spat = 0.5);
#' mfm_sbm(links$z, links$A, 0.8, links$dis, 1000)}
#' @export
mfm_sbm <- function(z, A, conc, S_ij, niter = 100){
    n <- nrow(A)
    ## Initialize node membership using Chinese Restaurant Process
    z <- nimble::rCRP(n = 1, conc, n)
    K <- length(unique(z))
    print(K)
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
#' @return Updated node membership
#' @note Internal helper
#' @export
update_z_from_q <- function(z, Q, A, conc){
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
        v <- dpois(K + 1, 1)/dpois(K, 1) * conc
        probs[i, K + 1] <- v * prod(new_clus)
        probs[i, ] <- probs[i, ]/sum(probs[i, ])
    }
    z <- apply(probs, 1, rcat, n = 1)
    print(table(z))
    return(serialize(z))
}

#' Condense Adjacency into SBM with Known Node Membership
#' 
#' Function to calculate the implied stochastic block probabilities
#' @param A The observed adjacency matrix
#' @param z The true node membership (values assumed to be in [1, K])
#' @return A K x K matrix with the observed density of each block
#' @examples \code{set.seed(123)
#' links <- gen_az(n_nodes = 50, K = 2, m = 2, prob = c(0.5, 0.5),
#'beta0 = 1, beta = diag(2) - 3, sigma = 0.3, spat = 0.5)
#' find_sbm(links$A, links$z)}
#' @export
find_sbm <- function(A, z){
  K <- max(z)
  sbm <- matrix(0, K, K)
  for (i in 1:K){
    for (j in 1:K){
      sbm[i, j] <- mean(A[z == i, z == j])
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
#' @return Entire MCMC samples data reorganized according to constraint
#' @examples
#' set.seed(123)
#' links <- generate_calfsbm_network(100, 2, 2, c(0.5, 0.5), diag(2) - 3, 
#' 1, 0.3, 1.5)
#' X <- calf_sbm_nimble(links, 1000, 500, 2, 3, 2)
#' post_label_mcmc_samples(X$mcmcSamples, 2, 1)
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
#' @param directed logical; if \code{FALSE} (default), the MCMC output is from an 
#' undirected network
#' @param lab logical; if \code{TRUE} (default), labels the clusters from 
#' smallest to largest mean values of within-cluster effects beta_{ii}
#' @return Reorganized dataframe according to apparent beta_{ii} ordering
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
#' @param links List of elements of the network, 
#' requires adjacency matrix A, matrix of covariates X, and distance matrix dis
#' @param nsim Total number of MCMC iterations per chain
#' @param burnin Number of iterations in each chain to be discarded
#' @param thin Post-burnin thinning parameter
#' @param nchain Number of MCMC chains to run
#' @param K Number of clusters
#' @param offset logical (default = \code{TRUE}); where \code{TRUE} 
#' indicates to use offset terms \theta in the \code{NIMBLE} model
#' @param beta_scale Prior standard deviation of beta terms
#' @param return_gelman logical (default = \code{FALSE}); if \code{TRUE}, 
#' returns the Gelman-Rubin diagnostic for all \beta terms
#' @return List of beta, z, and history of K
#' @export
calf_sbm_nimble <- function(links, nsim, burnin, thin, nchain, K, 
                            offset = TRUE, beta_scale = 10, 
                            return_gelman = FALSE){
  ## Inits
  const <- list(n = nrow(links$A), K = K)
  data <- list(A = links$A, x = links$dis)
  directed <- !isSymmetric(links$A)
  inits <- list(beta0 = rnorm(1, 0, 5)
                , beta = rnorm(const$K^2, 0, 5)
                , z = cluster::pam(links$X, const$K)$clustering
                , gamma = matrix(1, const$n, const$K)
  )
  if (offset){
    if (!directed){
      inits$theta <- log(rowSums(links$A) * const$n / sum(links$A) + 0.0001)
    } else {
      inits$theta_in <- log(rowSums(links$A) * const$n / sum(links$A) + 0.0001)
      inits$theta_out <- log(colSums(links$A) * const$n / sum(links$A) + 0.0001)
    }
  }
  ## Initialize betas
  group <- gen_factor(inits$z, links$A, links$dis)
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
    param_names <- colnames(mcmcSamples$chain1)
    gelman.diag <- boa::boa.chain.gandr(
      mcmcSamples, 
      list(mcmcSamples$chain1 - Inf, mcmcSamples$chain1 + Inf), 
      alpha = 0.05, pnames = param_names[grep('beta', param_names)])
  }
  mcmcSamples <- rbind(mcmcSamples$chain1, mcmcSamples$chain2, mcmcSamples$chain3)
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
#' @return Serialized node membership
#' @examples \code{serialize(c(2, 2, 3, 3))}
#' @export
serialize <- function(z){
    z <- as.factor(z)
    levels(z) <- 1:length(levels(z))
    return(as.numeric(levels(z))[z])
}
