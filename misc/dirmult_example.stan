//
data {
  int K; // num dims
  int N; // num observations
  array[N, K] int y; // multinomial observations
}
parameters {
  real<lower=0> alpha;
  array[N] simplex[K] theta;
}
model {
  alpha ~ lognormal(0, 1); // or some other prior on alpha
  for (n in 1:N) {
    theta[n] ~ dirichlet(rep_row_vector(alpha, K));
    y[n] ~ multinomial(to_vector(theta[n]));
  }
}


