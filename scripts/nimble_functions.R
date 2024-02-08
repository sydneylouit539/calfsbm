########## CONTAINS FUNCTIONS COMPATIBLE WITH RNIMBLE ##########################
library(nimble)
############## BASE FUNCTIONS ###########
## Test function for nimblecode
Unique <- nimbleFunction(
    run = function(vec = double(1)){
        uni <- numeric(length(vec))
        vals <- 1
        for (i in 1:length(vec)){
            if (any(uni == vec[i]) == FALSE){
                uni[vals] <- vec[i]
                vals <- vals + 1
            }
        }
        return(uni[1:vals])
        returnType(double(1))
    }
)


update_z_from_beta_nimble <- nimbleFunction(
    run = function(z = double(1), beta0 = double(0), beta = double(2), 
                 X_ij = double(2), A = double(2), K = double(0)){
    n <- length(z)
    log_prob_mat <- matrix(0, nrow = n, ncol = K)
    updated_z <- rep(0, n)
    for (i in 1:n){
        for (j in 1:K){
            eta_ij <- beta0 + beta[j, z[-i]] * X_ij[i, -i]
            result_vec <- A[i, -i]
            print(result_vec)
            ## Calculate probabilities
            fit <- exp(eta_ij) / (1 + exp(eta_ij))
            ## Calculate log-likelihood
            loglik <- 0
            for (k in 1:(n - 1)){
                if (result_vec[k] == 1) {loglik <- loglik - log(fit[k])} 
                else {loglik <- loglik - log(1 - fit[k])}
            }
            log_prob_mat[i, j] <- loglik
            print(loglik)
        }
        ## Get relative probabilities of each row
        wts <- log_prob_mat[i, ]
        wts <- wts - max(wts) #  Safety check for extremely low likelihood
        log_prob_mat[i, ] <- exp(wts)/ sum(exp(wts))
        #updated_z[i] <- sample(1:K, 1, prob=log_prob_mat[i, ])
        updated_z[i] <- which(rmulti(1, 1, prob = log_prob_mat[i, ]) == 1)
    }
    ## Make cluster assignments from 1 to K in case a cluster was deleted
    #serial_z <- updated_z; cl <- 1
    #for (k in Cunique(updated_z)){
    #    updated_z[serial_z == k] <- cl
    #    cl = cl + 1
    #}
    return(updated_z)
    returnType(double(1))
    }
)

compileNimble(update_z_from_beta_nimble, showCompilerOutput = TRUE)
################################################################################
## Block of NIMBLE code (Assumes fixed number of clusters)
source('sbm_functions.R')
## K = number of clusters, n = number of nodes, m=number of covariates
K <- 4
n <- 2000
m <- 2
beta_0 <- 1
directed <- FALSE
offset <- TRUE
lab <- FALSE
#beta_tru <- diag(K) * 1 - 2
beta_tru <- diag(seq(1.4, 2, length = K)) - 3
prob <- rep(1/K, K)
#links <- gen_az(n, K, m, prob, beta_0, beta_tru, directed = directed)
## Optional: Set seed for replicability
set.seed(36)
links <- gen_az(n, K, m, prob, beta_0, beta_tru, 0.99 * diag(m) + 0.01, 1, directed)
## MCMC parameters
nsim <- 1000
burnin <- 500
thin <- 5
## NIMBLE inits
const <- list(n = n, K = K)#, alpha = rep(1, K))
data <- list(A = links$A, x = links$dis)

inits <- list(beta0 = rnorm(1, 0, 5)
            , beta = rnorm(const$K^2, 0, 5)
            , z = cluster::pam(links$X, K)$clustering
            , gamma = matrix(1, const$n, const$K)
#            , theta = log(rowSums(links$A) * const$n / sum(links$A) + 0.0001)
)

## Initialize betas
group <- gen_factor(inits$z, links$A, links$dis)
initial_beta <- update_beta(const$K, group$cluster)
inits$beta0 <- initial_beta$beta0
inits$beta <- c(initial_beta$beta)

## NIMBLE code
#monitors <- c('z', 'beta', 'beta0', 'alpha')
#if(offset){monitors <- c(monitors, 'sigma', 'theta')}
#if(directed){monitors <- c(monitors, 'sigma_in', 'sigma_out', 'theta_in', 'theta_out')}
monitors <- c('beta', 'beta0', 'z', 'sigma')
              #'sigma_in', 'sigma_out')
code <- nimbleCode({
  ## Priors for parameter matrix
  beta0 ~ dnorm(mean = -3, sd = 3)
  for (a in 1:K^2){
    beta[a] ~ dnorm(mean = -3, sd = 3)
  }
  ## Priors for offset
  if (offset){
    if(!directed){
      for (i in 1:n){
        theta[i] ~ dnorm(mean = 0, var = sigma)
      }
    sigma ~ dinvgamma(1, 1)
    } else {
      for (i in 1:n){
        theta_in[i] ~ dnorm(mean = 0, var = sigma_in)
        theta_out[i] ~ dnorm(mean = 0, var = sigma_out)
      }
    sigma_in ~ dexp(1)
    sigma_out ~ dexp(1)
    }
  }
  ## Node membership
  for (i in 1:n){
    alpha[i, 1:K] ~ ddirch(gamma[i, 1:K])
    z[i] ~ dcat(alpha[i, 1:K])
    #z[i] ~ dcat(alpha[i, 1:K])
  }
  ## Adjacency matrix from fitted values
  for (i in 1:n){
    if (!directed){
      for (j in (i+1):n){
        ## Undirected network
        A[i, j] ~ dbin(expit(beta0 + 
                   theta[i] + theta[j] + 
                   beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
      }
    } else {
      for (j in 1:n){
        A[i, j] ~ dbin(expit(beta0 + 
                   theta_in[i] + theta_out[j] + 
                   beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
      }
    }
  }# (z[i] - 1) * K + z[j]
})
## Compile model
model <- nimbleModel(code, 
                     constants = const, 
#                     dimensions = list(),
                     data = data, 
                     inits = inits,
                     calculate = FALSE,
                     check = FALSE
)
cmodel <- compileNimble(model)
## Compile MCMC sampler
modelConf <- configureMCMC(model, monitors = monitors, enableWAIC = TRUE)
modelMCMC <- buildMCMC(modelConf)
cmodelMCMC <- compileNimble(modelMCMC, project = model) #1 min
cmodelMCMC$run(nsim)


## WAIC
cmodelMCMC$getWAIC()

mcmcSamples <- as.matrix(cmodelMCMC$mvSamples)[seq(burnin, nsim, by = thin), ]


## Diagnostics
plot(mcmcSamples[, 1], ylim = c(-2.5, -0.5), type = 'l', col = 'blue')
lines(mcmcSamples[, 4], col = 'red')
nett::compute_mutual_info(links$z, mcmcSamples[400, const$K^2 + 2 + 1:const$n])
mclust::adjustedRandIndex(links$z, mcmcSamples[400, const$K^2 + 2 + 1:const$n])
## Post-process MCMC results
mcmcSamples[, 1:const$K^2] <- post_label_mcmc(mcmcSamples[, 1:const$K^2], const$K, directed, lab)



## Add MCMC diagnostics
interc <- const$K^2 + 1 + 100
interc_vals <- quantile(mcmcSamples[, interc], c(0.025, 0.975))
diag_vals <- quantile(mcmcSamples[, 1], c(0.025, 0.975))
offdiag_vals <- quantile(mcmcSamples[, 2], c(0.025, 0.975))
nmi <- nett::compute_mutual_info(mcmcSamples[nrow(mcmcSamples), interc + 1:100], links$z)
nmi
table(mcmcSamples[nrow(mcmcSamples), interc + 1:100], links$z)
## MCMC trace plots and diagnostics

tail(mcmcSamples)


## CRP distribution implementation for improved flexibility over default
dcrp <- nimbleFunction(
  run = function(x = double(1), conc = double(0), size = double(0), 
                 log = integer(0, default = 0)){
    returnType(double(0))
    n <- length(x)
    if (n != size) {
      nimStop("dCRP: length of 'x' has to be equal to 'size'.\n")
    }
    if (conc <= 0 | is.na(conc)) {
      return(NaN)
    }
    if (any_na(x)) 
      return(NaN)
    ldens <- 0
    if (n > 1) {
      for (i in 2:n) {
        counts <- sum(x[i] == x[1:(i - 1)])
        if (counts > 0) {
          ldens <- ldens + log(counts/(i - 1 + conc))
        }
        else {
          ldens <- ldens + log(conc/(i - 1 + conc))
        }
      }
    }
    if (log) 
      return(ldens)
    else return(exp(ldens))
  }
)

rcrp <- nimbleFunction(
  run = function(n = double(0), conc = double(0), size = double(0)){
    returnType(double(1))
      if (n != 1) {
        nimStop("rCRP only handles n = 1 at the moment.\n")
      }
      if (conc <= 0 | is.na(conc)) {
        return(nimRep(NaN, size))
      }
      x <- nimNumeric(size)
      x[1] <- 1
      if (size > 1) {
        numComponents <- 1
        ones <- nimRep(1, size)
        for (i in 2:size) {
          if (runif(1) <= conc/(conc + i - 1)) {
            numComponents <- numComponents + 1
            x[i] <- numComponents
          }
          else {
            index <- rcat(n = 1, ones[1:(i - 1)])
            x[i] <- x[index]
          }
        }
      }
      return(x)
    }
)

