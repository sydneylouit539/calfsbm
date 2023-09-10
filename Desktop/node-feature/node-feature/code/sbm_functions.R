############## CODE CONTAINING FUNCTIONS FOR SBM.R #############################
library(bcdc)
library(igraph)
library(latentnet)
library(MASS)
library(nett)
library(network)
library(nimble)
library(raster)

## Visualize link probability vs. actual matrix
#' @param observed Observed adjacency matrix
#' @param prob Matrix of intrinsic probabilities
#' @param z_tru Vector of the true node membership
#' @return Graph of the connectivity and graph of true probabilities
plot_conn <- function(observed, prob, z_tru){
    n <- length(z_tru)
    par(mar = c(1, 1, 1, 1))
    par(mfrow = c(1, 2))
    B <- observed[order(z_tru), order(z_tru)]
    xx <- rasterFromXYZ(cbind(expand.grid(1:n, n:1), c(B)))
    plot(xx, legend=FALSE, axes=FALSE)
    ## Plot true probabilities
    B <- prob[order(z_tru), order(z_tru)]
    xx <- rasterFromXYZ(cbind(expand.grid(1:n, n:1), c(B)))
    plot(xx, axes=FALSE)
}


## Get the X matrix for logistic regression by taking clusters, links, and dist.
#' @param initial_z Initial node membership vector
#' @param A Known adjacency matrix
#' @param S_ij Distance matrix
#' @param directed Boolean 
#' @param offset Logistic regression offset terms (currently not used)
#' @return List of X matrix and grid of (x, y) indices
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


## Make a non-symmetric matrix symmetric through lower diagonal
#' @param mat The matrix to be symmetrized
#' @return Symmetric matric with the same upper diagonal as the original
makesymmetric <- function(mat){
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    return (mat)
}


## Generate X according to parameters
#' @param x_dim Dimension of X
#' @param n_nodes Number of nodes
#' @param sigma correlation (or correlation matrix if input is xdim by xdim)
#' @return Matrix of covariates
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
}


## Generate A and z from K and probabilities
#' @param n_nodes Number of nodes
#' @param K True number of groups
#' @param x_dim The number of covariates in X
#' @param prob Vector of probability of being in each group
#' @param beta0 Intercept term, lower values indicate higher sparsity
#' @param beta Matrix of logistic regression coefficients
#' @param sigma Correlation matrix of the X matrix of covariates
#' @param spat Degree of spatial heterogeneity of z
#' @param directed Boolean indicating whether the network should be directed
#' @param offset Boolean indicating whether the degree terms theta should be used
#' @return list of adjacency, membership, link probability, and distance
gen_az <- function(n_nodes, K, x_dim, prob, beta0, beta, sigma, 
                   spat, directed = FALSE, offset = FALSE){
    #z_tru <- sample(1:K, n_nodes, replace = TRUE, prob = prob)
    #X <- mvrnorm(n_nodes, mu = rep(0, x_dim), Sigma = 0.5 * diag(x_dim) + 0.5)
    X <- gen_x(x_dim, n_nodes, sigma)
    ## Add spatial correlation between groups
    mid <- mvrnorm(K, rep(0, x_dim), spat * diag(x_dim))
    z_tru <- rep(1, n_nodes)
    for (i in 1:n_nodes){
        z_tru[i] <- rcat(1, exp(-colSums((t(mid) - X[i, ])^2) / 2))
    }
    ## Degree
    theta <- rnorm(n_nodes, 0, 0.2)
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
                 dis = S_ij, theta = theta, beta0 = beta0, beta = beta_tru))
}


## Function to generate network according to Handcock's RSS latentnet paper 
#' @param n_nodes Number of nodes
#' @param K Number of clusters
#' @param x_dim Number of covariates
#' @param betavec Vector of intercept and slope terms
#' @param noise noise/signal ratio (default=1)
gen_az_ln <- function(n_nodes, K, x_dim, betavec, noise = 1, directed = FALSE){
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


## Function to find the optimal number of groups using latentnet package
#' @param p Number of clusters. The function will iterate from 1 to p
#' @param AA Network object containing adjacency and covariate information
#' @param sample_size Number of non-rejected MCMC samples
#' @param burnin Number of MCMC samples to allow the algorithm to converge
#' @param by Number of samples to skip after burn-in
#' @return A dataframe with the AIC and BIC values for each K value
find_K_optimal <- function(p, AA, sample_size, burnin, by){
    ## Set a dataframe for BIC and number of clusters
    best_estimate <- data.frame(K = 1:p, BIC = rep(Inf, p), AIC = rep(Inf, p))
    n <- AA$gal$n
    z_est <- rep(1, p)
    for (i in 1:p){
        samp.fit <- ergmm(AA ~ euclidean(d = 2, G = i), control = control.ergmm(
            sample.size = sample_size, burnin = burnin, interval = by), 
            verbose = TRUE)
        fit <- summary(samp.fit)$bic$overall
        if (fit < min(best_estimate$BIC)){
            z_est <- samp.fit$mkl$Z.K
        }
        best_estimate$BIC[i] <- fit
        best_estimate$AIC[i] <- fit + 2 * (i^2 + i)/2 - (i^2 + i)/2 * log(n * (n - 1) / 2)
        print(c(i, fit))
    }
    return (list(K = which.min(best_estimate$BIC), data = best_estimate, z = z_est))
}


## Function to expand K optimization into a simulation study
#' @param nsim The total number of simulations
#' @param n_vals Vector of the total node numbers
#' @param K True number of clusters, can also be vector form
#' @param p The maximum number of clusters to test for
#' @param x_dim Number of covariates in X
#' @param beta0 Intercept term
#' @param beta Effect of clusters on link probability
#' @return Matrix of true K vs. estimated K
sim_study_K_finder <- function(nsim, n_vals, K, p, x_dim, beta0, beta, 
                               sample_size = 2000, burnin = 4000, by = 4){
    derived_k <- rep(0, nsim)
    if (length(n_vals) == 1){n_vals = rep(n_vals, nsim)}
    if (length(K) == 1){K = rep(K, nsim)}
    for (i in 1:nsim){
        links <- gen_az(n_vals[i], K[i], x_dim, 
                        prob = rep(1/K[i], K[i]), beta0, beta[1:K[i], 1:K[i]])
        AA <- network(links$A, as.data.frame(links$X))
        sim_i <- find_K_optimal(p, AA, sample_size, burnin, by)
        derived_k[i] <- which.min(sim_i$BIC)
    }
    return(table(data.frame(n_tru = n_vals, k_tru = K, k_est = derived_k)))
}


## Update beta
#' @param K Number of clusters
#' @param group Model matrix to be fitted on
#' @param directed Boolean indicating whether network is directed
#' @return beta0 and beta, with beta as a matrix
update_beta <- function(K, group, directed = FALSE, offset = FALSE){
    mod_mat <- model.matrix(~ 0 + as.factor(cl), group) * group$x 
    mod_mat2 <- as.data.frame(mod_mat)
    mod_mat2$y <- group$y
    ## Fit Bayesian logistic regression to the data
    if (!offset){
      logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
                        prior.mean = 0, prior.scale = 10)
    } else {
      logit_fit <- arm::bayesglm(y ~ ., family = 'binomial', data = mod_mat2,
                        prior.mean = 0, prior.scale = 10, offset = group$offset)
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


## Update node membership using likelihood
#' @param z Initial values of z
#' @param X Matrix of covariates
#' @param beta0 Draw from posterior distribution
#' @param beta Draw from posterior distribution
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


## Checks MCMC for convergence
#' @param z Node membership from last five iterations
#' @param group 
#' @param niter Total iterations initially planned
#' @param iter The current iteration of the Gibbs sampler
#' @return Simulated coefficients
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


## Gibbs sampler with fixed number of clusters
#' @param K Number of clusters (should be optimized already)
#' @param alpha Vector of prior probabilities of being in each group
#' @param beta0 Prior intercept for logistic regression
#' @param beta Prior parameters for logistic regression, in matrix form
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @param directed Boolean indicating whether the network is directed
#' @return List containing the estimated beta, node membership, and history
gibbs_sampler_fixed_k <- function(K, alpha, beta0, beta, niter, A, S_ij, directed = FALSE){
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
        #if (i > 1){if (runif(1) < exp(bic_history[i - 1] - bic_history[i])) {i <- i + 1}}
        #if (i == 1){i <- i + 1}
        #print(c(i, bic_history[i - 1]))
        if (i %% 100 == 0){print(i)} #  Check to see if still running smoothly
    }
    return(list(beta = initial_beta, z = all_z[niter, ], prob_history = prob_history,#z_history = all_z, 
           bic = bic_history, int_history = intercepts, beta_history = beta_history, 
           converged = i))
}


## Gibbs sample with CRP as initial state
#' @param dp Dirichlet process parameter (recommended 1)
#' @param X Matrix of covariates
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @return List containing the estimated beta, node membership, and history
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


## Function that uses BIC to find optimal number of clusters
#' @param p Number of clusters to test (from 1 to p)
#' @param alpha Vector of prior probabilities of being in each group
#' @param beta0 Prior intercept for logistic regression
#' @param beta Prior parameters for logistic regression, in matrix form
#' @param niter Number of iterations to run
#' @param A Adjacency matrix
#' @param S_ij Distance matrix
#' @return List containing the estimated beta, node membership, and history
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
    return (list(K = which.min(best_estimate$BIC), data = best_estimate, z = z_est))
}


## Implementation of Algorithm 1 in MFM-SBM for the model
#' @param z node membership vector
#' @param A Known connectivity matrix
#' @param conc Concentration parameter of CRP
#' @param S_ij Distance
#' @param niter Number of iterations
#' @return List of beta, z, and history of K
mfm_sbm <- function(z, A, conc, S_ij, niter = 100){
    n <- nrow(A)
    ## Initialize node membership using Chinese Restaurant Process
    z <- rCRP(n = 1, conc, n)
    K <- length(unique(z))
    print(K)
    group <- gen_factor(z, A, S_ij, directed) 
    k_hist <- rep(0, niter)
    Q <- matrix(0.5, n, n)
    for (iter in 1:niter){
        ## Update beta and Q conditional on z
        beta <- update_beta(K, group$cluster, directed)
        for (i in 1:n){
          Q[i, ] <- expit(beta$beta0 + beta$beta[z[i], z[1:n]] * S_ij[i, ])
        }
        diag(Q) <- 0
        #for(i in 1:K){
        #  for (j in 1:K){
        #    c <- which(z == i)
        #    d <- which(z == j)
        #    Q[i, j] <- mean(expit(beta$beta0 + beta$beta[z[c], z[d]] * S_ij[c, d]))
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


## Helper function for MFM-SBM Algorithm 1
#' @param z Node membership
#' @param Q SBM probability matrix
#' @param A Known adjacency matrix
#' @param conc Concentration parameter
#' @return Updated node membership
update_z_from_q <- function(z, Q, A, conc){ #120ms
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
            new_clus[c] <- beta(a + a_star, b + length(which(z == c)) - a_star)/beta(a, b)
        }
        v <- dpois(K + 1, 1)/dpois(K, 1) * conc
        probs[i, K + 1] <- v * prod(new_clus)
        probs[i, ] <- probs[i, ]/sum(probs[i, ])
    }
    z <- apply(probs, 1, rcat, n = 1)
    print(table(z))
    return(serialize(z))
}


## Number clusters from 1-K in case a number was deleted
#' @param z Vector of node membership to be renumbered from 1-K
#' @return Serialized node membership
serialize <- function(z){
    z <- as.factor(z)
    levels(z) <- 1:length(levels(z))
    return(as.numeric(levels(z))[z])
}

