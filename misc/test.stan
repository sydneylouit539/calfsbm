data {
  int<lower=0> N;
  int<lower=0> K;
  real<lower=0> gamma;
  array[N, N] int<lower=0,upper=1> A;
  matrix[N, N] S_ij;
//  matrix[N, K] alpha;
  int<lower=0> a;
  int<lower=0> b;
}

// The parameters accepted by the model. 
parameters {
  array[N] simplex[K] z;
  vector[N] theta;
  real beta0;
  matrix[K, K] beta;
  real<lower=0> sigma;
//  matrix[N, K] z;
}

transformed parameters {
//  int z[N]; // Modes of the clustering
  matrix[N, N] eta;
  //matrix[N, K] z;
//  eta = beta0 + (z * beta * z');
  for (i in 1:N){
//    z[i, ] = row_simplex[K];
    for (j in 1:N){
      eta[i, j] = beta0 + theta[i] + theta[j] + 
        (to_row_vector(z[i]) * (beta * to_vector(z[j]))) * S_ij[i, j];
    }
  }
  
//      eta[i, j] = theta[i] + theta[j] + eta[i, j] * S_ij[i, j];
}

model {
  // Priors
    // Node heterogeneity
  sigma ~ inv_gamma(a, b);
  
  for (i in 1:N){
      // Degree heterogeneity
    theta[i] ~ normal(0, sigma);
      // Clustering assignments
    z[i, ] ~ dirichlet(rep_row_vector(gamma, K));
//    z[i, ] ~ dirichlet(alpha);
  }
    // Betas
  beta0 ~ normal(0, 10);
  for (i in 1:K){
    for (j in 1:K){
      beta[i, j] ~ normal(0, 10);
    }
  }
  // Likelihood
  for (i in 1:N){
    for (j in 1:N){
      // beta0 + theta[i] + theta[j] + (to_row_vector(z[i]) * (beta * to_vector(z[j]))) * S_ij[i, j]
      A[i, j] ~ bernoulli_logit(eta[i, j]);
    }
  }
//  print(z);
}
