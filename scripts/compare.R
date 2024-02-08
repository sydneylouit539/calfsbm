## SCRIPT TO ASSESS PERFORMANCE OF DIFFERENT MODELS ON DATASET
source('compare_base.R')
iden <- 1
set.seed(iden)


## Read in Weddell Sea dataset and take some minor cleaning steps
weddell <- read.csv('../../bcdc-master/Data/weddell_net.csv', row.names = 1)
full <- read.csv('../../bcdc-master/Data/weddell_cov.csv', sep = ';')

## Match covariates and adjacency
ind <- match(rownames(weddell), full$Species)
ind.null <- which(full[ind, "FeedingType"] == "NULL")
ind <- ind[-ind.null]

## Threshold adjacency matrix and make it symmetric
G <- as.matrix(weddell[-ind.null, -ind.null])
#A <- as(ifelse(tcrossprod(G) >= 5, 1, 0), "dgCMatrix")
#A <- ifelse(tcrossprod(G) >= 25, 1, 0)
A <- as(G, "dgCMatrix")
#diag(A) <- 0
n <- nrow(G); K <- 4
#find_sbm(A, links$z)

Cov_unscaled <- log(as.numeric(full[ind, "BodyWeight"])) # for plot
Cov <- scale(Cov_unscaled, center = TRUE, scale = TRUE)

feed <- as.factor(full$FeedingType[ind])
levels(feed) <- c(1, 1, 1, 2, 3, 4)

## Assume links is a list with A being adjacency matrix, X being covariate(s), and z being true node membership
directed <- FALSE
offset <- TRUE
n <- length(ind); K <- 4
links <- list(
  A = G,
  z = as.numeric(as.vector(feed)),
  X = Cov,
  dis = abs(outer(c(Cov), c(Cov), '-'))
)
## Initialize dataframe object to save results
final <- data.frame(
  k_means = 0,
  k_medians = 0,
  spectral = 0,
  casc = 0,
  bcdc = 0,
  chsbm = 0
)



## K-MEANS
k_means <- kmeans(links$X, K)$cluster
final$k_means <- nett::compute_mutual_info(k_means, links$z)

## K-MEDIANS (OUR INITIALIZED STATE)
k_meds <- cluster::pam(links$X, K)$clustering
final$k_medians <- nett::compute_mutual_info(k_meds, links$z)

## SPECTRAL CLUSTERING
spectre <- nett::spec_clust(links$A, K)
final$spectral <- nett::compute_mutual_info(spectre, links$z)

## CASC (COVARIATE-ASSISTED SPECTRAL CLUSTERING)

ghost <- kmeans(getCascAutoSvd(links$A, cbind(links$X, ifelse(links$z <= K / 2, 1, 0)), 
                               K, enhancedTuning = TRUE)$singVec
         , centers = K, nstart = 20)$cluster
final$casc <- nett::compute_mutual_info(ghost, links$z)

## BCDC (BAYESIAN COMMUNITY DETECTION WITH COVARIATES)
ind <- match(rownames(weddell), full$Species)
ind.null <- which(full[ind, "FeedingType"] == "NULL")
ind <- ind[-ind.null]

G <- as.matrix(weddell[-ind.null, -ind.null])
A <- as(G, 'dgCMatrix')#as(ifelse(tcrossprod(G) >= 5, 1, 0), "dgCMatrix")
diag(A) <- 0
n <- nrow(A)

Cov_unscaled <- log(as.numeric(full[ind, "BodyWeight"])) # for plot
Cov <- scale(Cov_unscaled, center = TRUE, scale = TRUE)

train_ratio <- 1 ## Use full training adjacency
M <- matrix(0, n, n)
no_pairs <- choose(n, 2)
no_observe <- as.integer(train_ratio * no_pairs)
no_latent <- no_pairs - no_observe
M[upper.tri(M)] <- sample(c(rep(1, no_observe)
                            , rep(0, no_latent)))
M[lower.tri(M)] <- t(M)[lower.tri(M)]
M <- Matrix(M, sparse = TRUE)

test_idx <- which(triu(!M,1))

# w/ covariates
mod <- new(CovarPSBM, A, M, alpha = 1, beta = 1, dp_concent = 1)
mod$set_continuous_features(Cov)
n_iter <- 2000
zmat <- mod$run_gibbs(n_iter)
final$bcdc <- nett::compute_mutual_info(zmat[, n_iter], links$z)
################### OUR METHODOLOGY ############################################

## Set constants and parameters
directed <- FALSE
offset <- TRUE
n <- length(ind); K <- 4
links <- list(
  A = G,
  z = as.numeric(as.vector(feed)),
  X = Cov,
  dis = abs(outer(c(Cov), c(Cov), '-'))
)


## MCMC parameters
nsim <- 20000
burnin <- 10000
thin <- 10
## NIMBLE inits
const <- list(n = n, K = K)
data <- list(A = links$A, x = links$dis)

inits <- list(beta0 = rnorm(1, 0, 5)
              , beta = rnorm(const$K^2, 0, 5)
              , z = cluster::pam(links$X, const$K)$clustering
              , alpha = matrix(1/const$K, const$n, const$K)
              , theta = log(rowSums(links$A) * const$n / sum(links$A) + 0.01)
)

## Initialize betas
group <- gen_factor(inits$z, links$A, links$dis)
initial_beta <- update_beta(const$K, group$cluster)
inits$beta0 <- initial_beta$beta0
inits$beta <- c(initial_beta$beta)

## Initial parameters for priors
#beta0_mean <- inits$beta0
#abs_beta0 <- abs(inits$beta0)
#beta_mean <- mean(inits$beta)
#beta_sd <- abs(beta_mean)
## NIMBLE code
monitors <- c('z', 'beta', 'beta0')
if(offset){monitors <- c(monitors, 'sigma', 'theta')}

code <- nimbleCode({
  ## Priors for parameter matrix
  beta0 ~ dnorm(mean = -3, sd = 3)
  for (a in 1:K^2){
    beta[a] ~ dnorm(mean = -3, sd = 3)
  }
  ## Priors for offset
  if (offset){
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
      ## Undirected network
      if (!directed){
        A[i, j] ~ dbin(expit(beta0 + 
                               theta[i] + theta[j] + 
                               beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
      } else {
        A[i, j] ~ dbin(expit(beta0 + 
                               theta[i] + theta[j] + 
                               beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
      }
    }# (z[i] - 1) * K + z[j]
  }
})
## Compile model
model <- nimbleModel(code, const, data, inits, check = FALSE)
cmodel <- compileNimble(model)
## Compile MCMC sampler
modelConf <- configureMCMC(model, monitors, enableWAIC = TRUE)
modelMCMC <- buildMCMC(modelConf)
cmodelMCMC <- compileNimble(modelMCMC, project = model) #1 min
cmodelMCMC$run(nsim)

mcmcSamples_wd <- as.matrix(cmodelMCMC$mvSamples)[seq(burnin, nsim, by = thin), ]
## Post-process samples using label.switching library
mcmcSamples_wd <- post_label_mcmc_samples(mcmcSamples_wd, const$K, directed)


## Diagnostics
last_estimate <- mcmcSamples_wd[nrow(mcmcSamples_wd), K^2 + 2 + 1:n]
find_sbm(links$A, last_estimate)
colMeans(mcmcSamples_wd)

#nett::compute_mutual_info(links$z, last_estimate) # 0.285074
mclust::adjustedRandIndex(links$z, last_estimate) # 0.3408694
table(links$z, last_estimate) # 0.3408694

final$lyz <- nett::compute_mutual_info(links$z, last_estimate)
write.csv(, paste0('results/weddell/diag_', iden, '.csv'), row.names = FALSE)
write.csv(final, paste0('results/weddell/', iden, '.csv'), row.names = FALSE)


















