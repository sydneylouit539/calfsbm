## Run Stan Model
library(calfsbm)
library(cmdstanr)
library(rstan)

## Set seed and generate network
set.seed(123)
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
links$alpha <- c(1, 1)
links$gamma <- 1
## Set model path for Stan
#model_path <- file.path(cmdstan_path(), '../../Desktop/calfsbm/misc/calfsbm_model.stan')
model_path <- file.path(cmdstan_path(), '../../Desktop/calfsbm/misc/test.stan')
#fit1 <- stan(
#  file = 'calfsbm_model.stan', iter = 100, verbose = TRUE,
#  data = links
#)

mod <- cmdstan_model(model_path)

fit_vb <- mod$variational(data = links, seed = 123, draws = 1000) 
#                          init = list(list(beta = 0*diag(2))))

summ <- fit_vb$summary()

## Check thetas against true values
cor(summ$mean[2 + 1:links$N], links$theta)
## Check betas
tail(summ, links$K^2 + 2)





stan_code <- "
data{
  int<lower=0> N;
  int<lower=0> K;
}

parameters {
  array[N] simplex[K] alpha;
}

model {
  for (i in 1:N){
    alpha[i] ~ dirichlet(rep_vector(1, K));
  }
}
"
links <- list (
  N = 100,
  K = 3,
  y = t(rmultinom(100, 1, rep(1, 3)))
)

model_path <- file.path(cmdstan_path(), '../../Desktop/calfsbm/misc/dirmult_example.stan')

mod <- cmdstan_model(model_path)












