## Run Stan Model
library(calfsbm)
library(cmdstanr)
library(rstan)

## Set seed and parameters
set.seed(1)
K <- 2
m <- 2
n <- 400
beta0 <- 1
beta_tru <- 1.5 * diag(K) - 3

## Generate network
links <- calfsbm::sim_calfsbm(n_nodes = n, 
                              K = K, 
                              n_covar = m, 
                              prob = rep(1, K), # Unbalanced cluster sizes
                              beta0 = beta0, 
                              beta = beta_tru,
                              sigma = 0.2, 
                              spat = 2,
                              directed = FALSE)
## Necessary inputs for Stan model
links$N <- n
links$K <- K
links$S_ij <- links$dis
links$a <- links$b <- 1
links$gamma <- 1
## Set model path for Stan
#model_path <- file.path(cmdstan_path(), '../../Desktop/calfsbm/misc/calfsbm_model.stan')
#model_path <- file.path(cmdstan_path(), '../../Desktop/community-detection/calfsbm-vb/code/test.stan')
links$theta <- log(rowSums(links$A) + 0.0001)
model_path <- file.path(cmdstan_path(), '../../Desktop/community-detection/calfsbm-vb/code/fixed_z.stan')
#model_path <- file.path(cmdstan_path(), '../../Desktop/community-detection/calfsbm-vb/code/no_theta.stan')
mod <- cmdstan_model(model_path)

## Initialization for soft z
initial_z <- cluster::pam(links$X, links$K)$clustering
c_1 <- ifelse(initial_z == 1, 0.8, 0.2)

## Variational fit
fit_vb <- mod$variational(data = links, seed = 123, draws = 1000,
                          algorithm = 'fullrank',
                          #init = list(
                            #list(z = cbind(c_1, 1 - c_1)),
                          #  list(z = cbind(c_1, 1 - c_1))
                          #  )
                          )

initial_beta <- initialize_beta_from_z(initial_z[1:subset, ], 
                                       links$dis[1:subset, 1:subset], 
                                       links$A[1:subset, 1:subset])
## Pathfinder fit
fit_pf <- mod$pathfinder(data = links, seed = 123,
                         init = list(list(z = cbind(c_1, 1 - c_1),
                                          beta = initial_beta$beta),
                                     list(z = cbind(c_1, 1 - c_1)),
                                     list(z = cbind(c_1, 1 - c_1)),
                                     list(z = cbind(c_1, 1 - c_1))),
                         )

summ <- fit_vb$summary()

## Check clustering assignments
cor(summ$mean[2 + 1:links$N], links$z)
## Check thetas against true values
cor(links$theta, summ$mean[2 + links$K * links$N + 1:links$N])
plot(links$theta, summ$mean[2 + links$K * links$N + 1:links$N],
     xlab = 'True Theta', ylab = 'Predicted Theta')
abline(0, 1)
## Check betas
tail(summ, links$K^2 + 2)




