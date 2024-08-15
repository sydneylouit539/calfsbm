//
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// Function
//functions {
//  int which_max(vector x){
//    real max_value = max(x);
//    int the_max = 0;
//    for (i in 1:size(x))
//      if (x[i] == max_value)
//        the_max = i;
//    return the_max;
//  }
//}


data {
  int<lower=0> N;
  int<lower=0> K;
//  vector[N] y;
//  matrix[N, N] A;
  array[N, N] int<lower=0,upper=1> A;
  matrix[N, N] S_ij;
//  matrix[N, K] alpha;
  int<lower=0> a;
  int<lower=0> b;
//  array[N] int<lower=1,upper=K> z;
}

// The parameters accepted by the model. 
parameters {
  vector[N] theta;
  real beta0;
  matrix[K, K] beta;
  //vector[N] z;
  //int z[N];
  real<lower=0> sigma;
//  matrix[K, N] soft_z; // log unnormalized clusters
}

//transformed parameters {
//  int z[N]; // Modes of the clustering
//  for (n in 1:N)
//    for (k in 1:K)
//      soft_z[n, k] = -(beta0 + beta[k, z[-n]] * S_ij[n, -n]);
      //eta_ij <- beta0 + beta[k, z[-n]] * S_ij[n, -n]
      // Calculate probabilities
      //fit <- 1 / (1 + exp(-eta_ij))
      // Calculate log-likelihood
      //loglik <- sum(log(ifelse(A[n, -n] == 1, fit, 1 - fit)));
//}

model {
  // Declare parameters
  matrix[K, N] soft_z;
  
  // Node heterogeneity
  sigma ~ inv_gamma(a, b);
  for (i in 1:N)
    theta[i] ~ normal(0, sigma);
  
  // Betas
  for (i in 1:K){
    for (j in 1:K){
      beta[i, j] ~ normal(0, 10);
    }
  }

  for (ii in 1:N) {
    // Node assignment
    z[ii] ~ categorical(softmax(col(soft_z, ii)));
    for (k in 1:K) {
//      soft_z[k, ii] = -(beta0 + beta[k, z[-ii]] * S_ij[ii, -ii]);
//      soft_z[k, ii] = -(beta0 + beta[k, z] * to_vector(S_ij[ii]));
// print(beta0 + to_vector(beta[k, z]) .* to_vector(S_ij[ii]));// This should be some kind of vector, not an integer
// print(inv_logit(beta0 + beta[k, z] * to_vector(S_ij[ii])));
        soft_z[k, ii] = sum(log(abs(to_vector(A[ii]) - 
              inv_logit(beta0 + to_vector(beta[k, z]) .* to_vector(S_ij[ii])))));
    }
//    for (jj in 1:N) {
      // Link probability
//      A[ii, jj] ~ bernoulli_logit(beta0 + 
//      theta[ii] + theta[jj] + 
//      beta[z[ii], z[jj]] * S_ij[ii, jj]);
//    }
  }
//  print(z);

}

// Working Soft K-Means example

//data {
//  int<lower=0> N;        // number of data points
//  int<lower=1> D;        // number of dimensions
//  int<lower=1> K;        // number of clusters
//  array[N] vector[D] y;  // observations
//}


//parameters {
//  array[K] vector[D] mu; // cluster means
//}

//transformed parameters {
//  array[N, K] real<upper=0> soft_z; // log unnormalized clusters
//  for (n in 1:N) {
//    for (k in 1:K) {
//      soft_z[n, k] = -log(K)
//                     - 0.5 * dot_self(mu[k] - y[n]);
//    }
//  }
//}
//model {
  // prior
//  for (k in 1:K) {
//    mu[k] ~ std_normal();
//  }

  // likelihood
//  for (n in 1:N) {
//    target += log_sum_exp(soft_z[n]);
//  }
//}
