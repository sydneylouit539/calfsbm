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
}


model {
  // Priors
    // Node heterogeneity
  sigma ~ inv_gamma(a, b);
    // Betas
  beta0 ~ normal(0, 10);
  for (i in 1:K){
    for (j in i:K){
      beta[i, j] ~ normal(0, 10);
    }
  }
  // Likelihood
  for (i in 1:N){
      // Degree heterogeneity
    theta[i] ~ normal(0, sigma);
      // Clustering assignments
    z[i, ] ~ dirichlet(rep_row_vector(gamma, K));
    for (j in i:N){
      A[i, j] ~ bernoulli_logit(beta0 + theta[i] + theta[j] + 
        (to_row_vector(z[i, ]) * (beta * (z[j, ]))) * S_ij[i, j]);
    }
  }
  // print(to_row_vector(z[10, ]) * (beta * (z[10, ])));
}
