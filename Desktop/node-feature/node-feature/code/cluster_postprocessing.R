##### SHORT CODE BLOCK TO PARSE RAW OUTPUT FROM THE CLUSTER ####################
K <- 2
base_path <- paste0('../../Archives/fixed_k_raw_data/', K, '_groups/all_data_')

sim_mean <- matrix(0, nrow = 100, ncol = K^2 + 1)
sim_se <- matrix(0, nrow = 100, ncol = K^2 + 1)
sim_low <- matrix(0, nrow = 100, ncol = K^2 + 1)
sim_high <- matrix(0, nrow = 100, ncol = K^2 + 1)

for (i in 1:100){
  replicate <- read.csv(paste0(base_path, i, '_', K, '_groups.csv'))
  for (j in 2:(K^2 + 2)){
    sim_mean[i, j - 1] <- mean(replicate[, j])
    sim_se[i, j - 1] <- sd(replicate[, j])
    sim_low[i, j - 1] <- quantile(replicate[, j], 0.025)
    sim_high[i, j - 1] <- quantile(replicate[, j], 0.975)
  }
}

true_values <- c(1.5 * diag(K) - 3, 2)

for (i in 1:ncol(sim_low)){
  coverage <- length(which(sim_low[, i] <= true_values[i]
                           & true_values[i] <= sim_high[, i]))
  print(coverage / 100)
}
> postprocess_results(2)
[1] 0.39
[1] 0.37
[1] 0.31
[1] 0.39
[1] 0.46
> postprocess_results(3)
[1] 0.38
[1] 0.36
[1] 0.37
[1] 0.36
[1] 0.43
[1] 0.39
[1] 0.36
[1] 0.36
[1] 0.38
[1] 0.45








