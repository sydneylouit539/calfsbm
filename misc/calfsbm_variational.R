## Run Stan Model
library(calfsbm)
library(cmdstanr)
library(rstan)

links <- calfsbm::sim_calfsbm(n_nodes = 200, 
                              K = 2, 
                              n_covar = 2, 
                              prob = c(0.5, 0.5), 
                              beta0 = 1, 
                              beta = diag(2) - 3,
                              sigma = 0.3, 
                              spat = 2)
links$N <- length(links$theta)
links$K <- nrow(links$beta)
links$S_ij <- links$dis
links$a <- links$b <- 1


model_path <- file.path(cmdstan_path(), '../../Desktop/calfsbm/misc/calfsbm_model.stan')
#fit1 <- stan(
#  file = 'calfsbm_model.stan', iter = 100, verbose = TRUE,
#  data = links
#)

mod <- cmdstan_model(model_path)

fit_vb <- mod$variational(data = links, seed = 123, draws = 20000)

summ <- fit_vb$summary()

## Check thetas against true values
cor(summ$mean[3:(2 + links$N)], links$theta)
## Check betas
tail(summ)

