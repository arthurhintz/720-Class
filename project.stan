data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> J;
  array[N] vector[2] coord;
  array[N] int<lower=1,upper=K> trt;
  array[N] int<lower=1,upper=J> blocks;
  vector[N] y;   // don`t used
}
  
parameters {
  real mu_0;
  vector[K] T;
  vector[J] B;
  real<lower=0> range;
  real<lower=0> psill;
  //real<lower=0> nugget;
  vector[N] eta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  matrix[N, N] L_K;
  {
    matrix[N, N] Kcov = gp_exp_quad_cov(coord, psill, range);
    for (n in 1:N)
      Kcov[n, n] += 1e-9;           
    L_K = cholesky_decompose(Kcov);
  }
  for (i in 1:N)
    mu[i] = mu_0 + T[trt[i]] + B[blocks[i]] + eta[i];
}

model {

  mu_0  ~ normal(0, 10);
  
  for (k in 1:3) T[k] ~ normal(1, 10);
  for (k in 4:7) T[k] ~ normal(-1, 10);
  for (k in 8:12) T[k] ~ normal(5, 10);
  for (k in 13:K) T[k] ~ normal(-5, 10);
  
  B ~ normal(0, 5);
  range ~ gamma(27, 1.1);
  psill  ~ gamma(9, 0.6);
  sigma  ~ gamma(50, 5);
  eta   ~ multi_normal_cholesky(rep_vector(0, N), L_K);
  y ~ normal(mu, sigma);
  
  //entropy() in R
  
}

generated quantities {
  
  real mu_0;
  vector[K] T;
  vector[J] B;
  real<lower=0> range;
  real<lower=0> psill;
  real<lower=0> sigma;
  vector[N] eta;
  vector[N] mu;
  vector[N] y_rep;

  mu_0 = normal_rng(0, 10);
  
  for (k in 1:3) T[k] = normal_rng(1, 10);
  for (k in 4:7) T[k] = normal_rng(-1, 10);
  for (k in 8:12) T[k] = normal_rng(5, 10);
  for (k in 13:K) T[k] = normal_rng(-5, 10);
  
  for (j in 1:J) B[j] = normal_rng(0, 5);

  range = gamma_rng(27, 1.1);
  psill  = gamma_rng(9, 0.6);
  sigma  = gamma_rng(50, 5);

  {
    matrix[N, N] Kcov = gp_exp_quad_cov(coord, psill, range);
    matrix[N, N] L_K;
    for (n in 1:N) Kcov[n, n] += 1e-9;
    L_K = cholesky_decompose(Kcov);
    eta = multi_normal_cholesky_rng(rep_vector(0, N), L_K);
  }

  for (i in 1:N)
    mu[i] = mu_0 + T[trt[i]] + B[blocks[i]] + eta[i];

  for (n in 1:N)
    y_rep[n] = normal_rng(mu[n], sigma);
}

