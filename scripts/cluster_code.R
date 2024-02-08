## Block of code to be copy/pasted into the cluster
library(nimble)
library(MASS)
## Batch identifier for .csv writing purposes
iden <- commandArgs(trailingOnly = TRUE)
## Generate data
source('code/sbm_functions.R')
K <- 4
n <- 200
m <- 2
beta_0 <- 1
directed <- FALSE
offset <- TRUE
lab <- TRUE
## Generate true betas in a label-switching scenario
if (!lab){
  beta_tru <- diag(K) * 1.5 - 3
} else {
  beta_tru <- diag(seq(1, 2, length = K)) - 3
}
prob <- rep(1/K, K)
set.seed(iden)
links <- gen_az(n, K, m, prob, beta_0, beta_tru, 0.7 * diag(m) + 0.3, 5, directed, offset)
## MCMC parameters
nsim <- 15000
burnin <- 5000
thin <- 10
n_chain <- 3
## NIMBLE inits
const <- list(n = n, K = K)
data <- list(A = links$A, x = links$dis)

inits <- list(beta0 = rnorm(1, 0, 5)
              , beta = rnorm(const$K^2, 0, 5)
              , z = cluster::pam(links$X, const$K)$clustering
              , alpha = matrix(1/const$K, const$n, const$K)
              , theta = log(rowSums(links$A) * const$n / sum(links$A) + 0.0001)
)

## Initialize betas
group <- gen_factor(inits$z, links$A, links$dis)
initial_beta <- update_beta(const$K, group$cluster)
inits$beta0 <- initial_beta$beta0
inits$beta <- c(initial_beta$beta)

## NIMBLE code
monitors <- c('z', 'beta', 'beta0')
if (offset){monitors <- c(monitors, 'sigma', 'theta')}

code <- nimbleCode({
  ## Priors for parameter matrix
  beta0 ~ dnorm(mean = 0, sd = 3)
  for (a in 1:K^2){
    beta[a] ~ dnorm(mean = 0, sd = 3)
  }
  ## Priors for offset
  if (offset) {
    for (i in 1:n){
      theta[i] ~ dnorm(mean = 0, var = sigma)
    }
    sigma ~ dexp(1)
  }
  ## Node membership
  for (i in 1:n){
    z[i] ~ dcat(alpha[i, 1:K])
  }
  ## Adjacency matrix from fitted values
  for (i in 1:n){
    for (j in (i+1):n){
      if (!directed){
        A[i, j] ~ dbin(expit(beta0 +
                               theta[i] + theta[j] +
                               beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
      } else {
        A[i, j] ~ dbin(expit(beta0 +
                               theta[i] + theta[j] +
                               beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
      }
    }
  }
})


## Compile model
model <- nimbleModel(code, const, data, inits, check = FALSE)
cmodel <- compileNimble(model)
## Compile MCMC sampler
modelConf <- configureMCMC(model, monitors)
modelMCMC <- buildMCMC(modelConf)
cmodelMCMC <- compileNimble(modelMCMC, project = model) # 1 min

## Multiple chains runner
mcmcSamples <- nimble::runMCMC(cmodelMCMC, niter = nsim, nburnin = burnin,
                               thin = thin, nchains = n_chain)
mcmcSamples <- rbind(mcmcSamples$chain1, mcmcSamples$chain2, mcmcSamples$chain3)#[, 1:(const$K^2 + 2)]
## Post-process MCMC results for label-switching
#mcmcSamples[, 1:const$K^2] <- post_label_mcmc(mcmcSamples[, 1:const$K^2], const$K, directed)
mcmcSamples <- post_label_mcmc_samples(mcmcSamples, const$K, directed)

#cmodelMCMC$run(nsim)
#mcmcSamples <- as.matrix(cmodelMCMC$mvSamples)[seq(burnin, nsim, by = thin), ]

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

