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
n <- 200
m <- 2
## beta_0 = 
beta_0 <- 1
directed <- FALSE
beta_tru <- diag(K) * 1.5 - 3
#prob <- rep(1/K, K)
#links <- gen_az(n, K, m, prob, beta_0, beta_tru, directed = directed)
links <- gen_az(n, K, m, 0, beta_0, beta_tru, 0.6 * diag(m) + 0.4, 1, directed)
## MCMC parameters
nsim <- 5000
burnin <- 2500
thin <- 5
## NIMBLE code
code <- nimbleCode({
  ## Priors for parameter matrix
  beta0 ~ dnorm(mean = 0, sd = 10)
  for (a in 1:K^2){
    beta[a] ~ dnorm(mean = 0, sd = 10)
  }
  ## Priors for offset
  for (i in 1:n){
    theta[i] ~ dnorm(mean = 0, sd = 10)
  }
  ## Node membership
  for (i in 1:n){
    z[i] ~ dcat(alpha[i, 1:K])
  }
  ## Adjacency matrix from fitted values
  for (i in 1:n){
    for (j in (i+1):n){
      A[i, j] ~ dbin(expit(beta0 + 
                             theta[i] + theta[j] + 
                             beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
    }
  }
})

set.seed(1)
const <- list(n = 100, K = 4)
const$inds <- expand.grid(1:const$n, 1:const$n)
data <- list(A = links$A, x = links$dis
             #             , inds = expand.grid(1:const$n, 1:const$n)
)

inits <- list(beta0 = rnorm(1, 0, 10)
              , beta = rnorm(const$K^2, 0, 10)
              , theta = rnorm(const$n, 0, 1)
#              , z = sample(1:const$K, const$n, replace = TRUE)
              , z = cluster::pam(links$X, const$K)$clustering
              , alpha = matrix(1/const$K, const$n, const$K)
)
table(inits$z)
## Compile model
model <- nimbleModel(code, const, data, inits, check = FALSE)
cmodel <- compileNimble(model)
monitors <- c('z', 'beta', 'beta0', 'theta')
## Compile MCMC sampler
modelConf <- configureMCMC(model, monitors)
modelMCMC <- buildMCMC(modelConf)
cmodelMCMC <- compileNimble(modelMCMC, project = model) #1 min
cmodelMCMC$run(nsim)

mcmcSamples <- as.matrix(cmodelMCMC$mvSamples)[seq(burnin, nsim, by = thin), ]
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

