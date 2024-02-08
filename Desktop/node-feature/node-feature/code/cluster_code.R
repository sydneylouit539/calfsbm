## Block of code to be copy/pasted into the cluster
library(nimble)
library(MASS)
## Batch identifier for .csv writing purposes
iden <- commandArgs(trailingOnly = TRUE)
source('code/sbm_functions.R')
## Set constants
K <- 2
n <- 100
m <- 2
beta_0 <- 1
spat <- 1
directed <- FALSE
beta_tru <- diag(K) * 1.5 - 3
prob <- rep(1/K, K)
## Include offset terms theta_i and theta_j?
offset <- TRUE
## Generate data
links <- gen_az(n, K, m, prob, beta_0, beta_tru, 0.6 * diag(m) + 0.4, 
                spat, directed, offset)
## MCMC parameters
nsim <- 50000
burnin <- 10000
thin <- 10
n_chain <- 1
## NIMBLE inits
const <- list(n = n, K = K)
data <- list(A = links$A, x = links$dis)

inits <- list(beta0 = rnorm(1, 0, 5)
              , beta = rnorm(const$K^2, 0, 5)
              , z = cluster::pam(links$X, K)$clustering
              , alpha = matrix(1/const$K, const$n, const$K)
)
## Initialize betas
group <- gen_factor(inits$z, links$A, links$dis)
initial_beta <- update_beta(const$K, group$cluster)
inits$beta0 <- initial_beta$beta0
inits$beta <- c(initial_beta$beta)


## NIMBLE code
monitors <- c('z', 'beta', 'beta0')

if (offset) {
  code <- nimbleCode({
    ## Priors for parameter matrix
    beta0 ~ dnorm(mean = 0, sd = 3)
    for (a in 1:K^2){
      beta[a] ~ dnorm(mean = 0, sd = 3)
    }
    ## Priors for offset
    for (i in 1:n){
      theta[i] ~ dnorm(mean = 0, sd = 1)
    }
    ## Node membership
    for (i in 1:n){
      z[i] ~ dcat(alpha[i, 1:K])
    }
    ## Adjacency matrix from fitted values
    for (i in 1:n){
      for (j in (i+1):n){
        A[i, j] ~ dbin(expit(beta0 +
            theta[i] + theta[j] + beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
      }
    }
  })
  ## Initialize theta
  inits$theta <- rnorm(const$n, 0, 0.2)
  monitors <- c(monitors, 'theta')
} else {
  code <- nimbleCode({
    ## Priors for parameter matrix
    beta0 ~ dnorm(mean = 0, sd = 5)
    for (a in 1:K^2){
      beta[a] ~ dnorm(mean = 0, sd = 5)
    }
    ## Node membership
    for (i in 1:n){
      z[i] ~ dcat(alpha[i, 1:K])
    }
    ## Adjacency matrix from fitted values
    for (i in 1:n){
      for (j in (i+1):n){
        A[i, j] ~ dbin(expit(beta0 + beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
      }
    }
  })
  
}


## Compile model
model <- nimbleModel(code, const, data, inits, check = FALSE)
cmodel <- compileNimble(model)
## Compile MCMC sampler
modelConf <- configureMCMC(model, monitors)
modelMCMC <- buildMCMC(modelConf)
cmodelMCMC <- compileNimble(modelMCMC, project = model) #1 min
## ~ 30s per 1000 MCMC
sim_results <- data.frame(
  intercept.low = rep(0, 1),
  intercept.high = rep(0, 1),
  within.low = rep(0, 1),
  within.high = rep(0, 1),
  between.low = rep(0, 1),
  between.high = rep(0, 1),
  NMI = rep(0, 1)
)
#for (i in 1:10){
cmodelMCMC$run(nsim)
mcmcSamples <- as.matrix(cmodelMCMC$mvSamples)[seq(burnin, nsim, by = thin), ]

write.csv(mcmcSamples, paste0('all_data_', iden, '_', K, '_groups.csv'), row.names = FALSE)
## Add MCMC diagnostics
interc <- const$K^2 + 101
interc_vals <- quantile(mcmcSamples[, interc], c(0.025, 0.5, 0.975))
diag_vals <- quantile(mcmcSamples[, 1], c(0.025, 0.5, 0.975))
offdiag_vals <- quantile(mcmcSamples[, 2], c(0.025, 0.5, 0.975))
interc_sd <- sd(mcmcSamples[, interc])
diag_sd <- sd(mcmcSamples[, 1])
offdiag_sd <- sd(mcmcSamples[, 2])

nmi <- nett::compute_mutual_info(mcmcSamples[nrow(mcmcSamples), interc + 1:100], links$z)
print(nmi)
table(mcmcSamples[nrow(mcmcSamples), interc + 1:100], links$z)
#tab <- apply(mcmcSamples[1:nrow(mcmcSamples), interc + 1:100], 1, table)
#if (class(tab)==c('matrix','array')){print(1)}
#K.vals <- unlist(lapply(apply(mcmcSamples[1:nrow(mcmcSamples), interc + 100 + 1:100], 1, table), length))
sim_results[1, ] <- c(interc_vals, interc_sd, diag_vals, diag_sd, offdiag_vals, offdiag_sd, nmi)
#}
#write.csv(sim_results, paste0('sim_results', iden, '.csv'), row.names = FALSE)

