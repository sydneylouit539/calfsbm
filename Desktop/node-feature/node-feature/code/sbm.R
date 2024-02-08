library(bcdc)
library(igraph)
library(latentnet)
library(MASS)
library(nett)
library(network)
library(raster)


source('sbm_functions.R')

####################### SET MODEL PARAMETERS ###################################
## Set number of nodes (n) and true number of clusters (K)
K <- 3
n <- 100
## Number of covariates (m) and number of iterations in Gibbs sampler
m <- 2
niter <- 100
## Directed?
directed <- FALSE
## Set beta_ij
beta_0 <- 2
if(directed){beta_tru <- diag(rnorm(K, mean = 1.5, sd = 0.2)) + 
    rnorm(K^2, -3, sd = 0.2)
} else {
    beta_tru <- diag(K) * 2 - 4
}

############################# TEST CASES #######################################
## Balanced case
prob <- rep(1/K, K)
links <- gen_az(n, K, m, prob, beta_0, beta_tru, directed = directed)
degrees <- rowSums(links$A)

## Gibbs sampler, with predicted K and z from earlier
K_guess <- 3
updated_beta_z <- gibbs_sampler_fixed_k(K_guess, rep(1, K_guess), beta0 = 0, 
  beta = matrix(0, K_guess, K_guess), niter = niter, links$A, links$dis, directed)

## MCMC trace plots
#par(mfrow = c(K_guess, K_guess))
plot(updated_beta_z$beta_history[83:niter,1,1], type='l')
cor(updated_beta_z$beta_history[83:(niter - 1), 1, 1], 
    updated_beta_z$beta_history[84:niter, 1, 1])
updated_beta_z <- gibbs_sampler_unknown_k(2, niter = 100, links$A, links$dis, directed)
#optimized_beta_z <- find_k_best_bic(4, rep(1, 4), beta0 = 0, 
#                        beta = matrix(0, 4, 4), niter = 100, links$A, links$dis)

nett::compute_mutual_info(updated_beta_z$z, links$z)
table(updated_beta_z$z, links$z)
#sbm_full_sim(niter, links$A, links$z, links$X, links$lp, plot = TRUE)




################################################################################
## Slightly unbalanced case
#prob <- (1:K) / K
#links <- gen_az(n, K, m, prob, beta_0, beta_tru)
#AA <- network::network(links$A, as.data.frame(links$X))
## Use latentnet package to find optimal cluster number by BIC
#head_start <- find_K_optimal(5, AA, size, burn, by)
## Gibbs sampler, with predicted K and z from earlier
#gibbs_sampler_fixed_k(head_start$K, prob, beta, links$X, niter, links$A)
#sbm_full_sim(niter, links$A, links$z, links$X, links$lp, plot = TRUE)
#mcmc.diagnostics(samp.fit)

## Heavily unbalanced case
#prob <- (1:K)^2 / K^2
#links <- gen_az(n, K, m, prob, beta_0, beta_tru)
#AA <- network::network(links$A, as.data.frame(links$X))
## Use latentnet package to find optimal cluster number by BIC
#head_start <- find_K_optimal(5, AA, size, burn, by)
## Gibbs sampler, with predicted K and z from earlier

#gibbs_sampler_fixed_k(head_start$K, prob, beta, links$X, niter, links$A)

#sbm_full_sim(niter, links$A, links$z, links$X, links$lp, plot = TRUE)

#################### K OPTIMIZATION SIMULATION STUDY ###########################
#ns <- rep(c(100), c(10))
#K_reals <- rep(2, c(10))
#beta_real <- diag(5)-3

#sim_study_K_finder(10, ns, K_reals, 5, x_dim = 3, beta0 = 2, beta = beta_real)



