## VISUALIZATIONS
library(data.table)
library(ggplot2)
library(maps)
library(patchwork)
library(sf)
library(xtable)

n <- c(200, 400); spat <- c(1, 3)
############### NODE ASSIGNMENT BOXPLOTS #######################################
## Working directory should be Desktop/community_detect_raw_data
## Canadian province shapefile should be in this folder
ari_base_200 <- ari_base_400 <- data.frame(matrix(0, nrow = 200, ncol = 5))
names(ari_base_200) <- names(ari_base_400) <- c('K-means', 'K-medians',
                                        'Spectral', 'CASC', 'CALF-SBM')
path <- '../Archives/fixed_k_raw_data/dynamic/'
files <- matrix(
  c(list.files(path, pattern = paste0(n[1], '_sn_', spat[1])),
    list.files(path, pattern = paste0(n[1], '_sn_', spat[2])),
    list.files(path, pattern = paste0(n[2], '_sn_', spat[1])),
    list.files(path, pattern = paste0(n[2], '_sn_', spat[2]))),
  nrow = 200, ncol = 2
)
for (h in 1:2){
  for (i in 1:200){
    file <- as.data.frame(fread(paste0(path, files[i, h])))
    #print(tr(table(c(file[1:400, 1])$truth, c(file[1:400, 6])$chsbm)))
    if(h == 1){ari_base_200[i, ] <- file[n[h] + 1, 2:6]}
    else{ari_base_400[i, ] <- file[n[h] + 1, 2:6]}
  }
}
ari <- expand.grid( rep(paste('S/N Ratio =', spat / 2), rep(100, 2)), 
                    names(ari_base_200), 
                    paste('n =', n))
#ari_200 <- expand.grid( rep(spat, rep(100, 2)), names(ari_base_200))
#ari_400 <- expand.grid( rep(spat, rep(100, 2)), names(ari_base_400))
#nmi <- expand.grid( rep(spat, rep(100, 2)), names(nmi_base))
names(ari) <- c('S/N Ratio', 'Method', 'Nodes')
#names(ari_200) <- names(ari_400) <- c('S/N Ratio', 'Method')
ari$ARI <- c(unname(unlist(c(ari_base_200, ari_base_400))))
#ari_200$ARI <- c(unname(unlist(ari_base_200)))
#ari_400$ARI <- c(unname(unlist(ari_base_400)))
ari$`S/N Ratio` <- as.factor(ari$`S/N Ratio`/2)
#ari_200$`S/N Ratio` <- as.factor(ari_200$`S/N Ratio`/2)
#ari_400$`S/N Ratio` <- as.factor(ari_400$`S/N Ratio`/2)
theme_set(theme_gray(base_size = 14))
#col_200 <- ifelse(ari_200$`S/N Ratio` == 0.5, 0.5, 1.5)
#a <- ggplot(ari_200, aes(x=`Method`, y=`ARI`,  fill = `Method`, color=`S/N Ratio`)) + 
#  geom_boxplot() + ylim(c(0, 1)) + 
#  theme(legend.position = 'none') +
#  ggtitle('n = 200')
#b <- ggplot(ari_400, aes(x=`Method`, y=`ARI`,  fill = `Method`, color=`S/N Ratio`)) + 
#  geom_boxplot() + ylim(c(0, 1)) + 
#  ggtitle('n = 400')
#gridExtra::grid.arrange(a, b, ncol = 2)
#combined <- a + b & theme(legend.position = "bottom")
#combined + plot_layout(guides = "collect")

ggplot(ari, aes(x=`Method`, y=`ARI`,  fill = `Method`)) + 
  geom_boxplot() + ylim(c(0, 1)) + 
  facet_grid(`Nodes` ~`S/N Ratio`) + 
#  ggtitle('Signal-to-Noise Ratio') +
  theme(legend.position = 'bottom')

############## K-SELECTION TABLES ##############################################
data_lib <- '../Archives/fixed_k_raw_data/dynamic/'
ln_files <- list.files(data_lib, pattern = 'latentnet')
ln_K <- rep(Inf, length(ln_files))
iter <- 1
for (i in ln_files){
  basedata <- as.data.frame(fread(paste0(ln_str, i)))
  ln_K[iter] <- basedata[which.min(basedata[, 2]), 1]
  iter <- iter + 1
}
print(table(ln_K))

## Evaluate BCDC file
bcdc_K <- rep(0, 100)
bcdc_file <- read.csv(paste0(data_lib, 'bcdc.csv'))
for (x in 1:100) {bcdc_K[x] <- length(table(unlist(bcdc_file[x, 401:800])))}
print(table(bcdc_K))

## Evaluate CHSBM files
chsbm_k <- rep(0, 100)
for (j in 1:100){
  base_files <- list.files(data_lib, pattern = glob2rx(paste0('chsbm_', j, '_*_groups*')))
  if (length(base_files) == 5){
    optimal_K <- Inf
    for(k in 1:5){
      basedata <- as.data.frame(fread(paste0(data_lib, base_files[k])))
      print(c(j, basedata[1, 1]))
      waic <- ifelse(is.na(basedata[1, 1]), Inf, basedata[1, 1])
      if (waic < optimal_K){
        optimal_K <- waic
        chsbm_k[j] <- round(sqrt(ncol(basedata) - 403))
      }
    }
  }
}
print(table(chsbm_k))
## Xtable summaries

xtable::xtable()




############ FIVE NUMBER SUMMARY FOR INITIALIZED RANDOM EFFECTS ################
n <- 400
K <- 4
m <- 2
beta0 <- 1
beta <- diag(seq(1.4, 2, length = K)) - 3
sigma <- 0.1
spat <- 3
directed <- FALSE
offset <- TRUE

predicted_means <- true_means <- predicted_vars <- true_vars <- cors <- rep(0, 100)
for (l in 1:100){
  set.seed(l)
  links <- gen_az(n, K, m, prob, beta0, beta, sigma, spat, directed, offset)
  total_variance <- var(links$X)
  print(total_variance)
  
  true_means[l] <- mean(links$theta)
  true_vars[l] <- var(links$theta)
  preds <- log(rowSums(links$A) * n / sum(links$A) + 0.01)
  predicted_means[l] <- mean(preds)
  predicted_vars[l] <- var(preds)
  cors[l] <- cor(links$theta, preds)
}
hist(cors)


######################### AIRPORT CLUSTERING ###################################
par(bg = 'grey100')
#five_k <- read.csv('real-data-results/airports_5_groups_alpha_prior.csv')
optimal_k <- read.csv('real-data-results/airports_8_groups_weak_prior.csv')
links <- readRDS('airports.RData')
n <- 456
z_est <- optimal_k$chsbm[1:n]; K_est <- length(unique(z_est))
## Label-switching (In case not done before)
z_est1 <- as.factor(z_est)
levels(z_est1)[order(diag(find_sbm(links$A, z_est)))] <- levels(z_est1)
z_est1 <- as.numeric(as.character(z_est1))
#pal <- palette.colors(palette = "R4")#[-c(5, 6)]
pal <- c('darkorange', 'red', 
         'lightgreen', 'purple', 
         'grey', 'darkblue', 
         'gold', 'lightblue')

## Plotting map with cities overlaid
canadian_provinces <- st_read('canada.shp')
plot(canadian_provinces$geometry, xlim = c(-162, -50), ylim = c(19, 65))
points(links$X[, 2], links$X[, 3], col = pal[z_est1], 
       cex = 0.1 * sqrt(colSums(links$A)), pch = 20)
monde <- maps::map('world', regions = c('Canada', 'USA:Alaska', 
                                        'USA:Hawaii', 'Mexico'), plot = FALSE)
etats <- maps::map('state', plot = FALSE)
lines(monde); lines(etats)
#maps::map('state', xlim = c(-162, -50), ylim = c(19, 65))
legend('bottomright', legend = paste('Cluster', 1:K_est), 
       col = pal[1:K_est], pch = 20, cex = 1)


########### MCMC TRACE PLOTS, GELMAN RUBIN ##############
mcmc_raw <- read.csv('../Archives/fixed_k_raw_data/dynamic/chsbm_1_4_groups.csv')
plot.ts(mcmc_raw$beta.1.[1:1000])
lines(1:1000, mcmc_raw$beta.1.[1001:2000], col = 'green')
lines(1:1000, mcmc_raw$beta.1.[2001:3000], col = 'orange')


