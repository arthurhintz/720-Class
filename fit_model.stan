data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> J;
  array[N] vector[2] coord;
  array[N] int<lower=1,upper=K> trt;
  array[N] int<lower=1,upper=J> blocks;
  vector[N] y;
}
  
parameters {
  real mu_0;
  vector[K] T;
  vector[J] B;
  real<lower=0> range;
  real<lower=0> psill;
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
  
  //Intercept
  mu_0  ~ normal(0, 10);
  
  
  // Treatmetns
  T[1] ~ normal(-4, 2);
  for (k in 2:5) T[k] ~ normal(-2, 2);
  for (k in 6:9) T[k] ~ normal(0, 2);
  for (k in 10:12) T[k] ~ normal(2, 2);
  for (k in 13:K) T[k] ~ normal(4, 2);
  
  
  // Blocks
 if (J == 1) {
    B[1] ~ normal(0, 0.0001);
  } else {
    B[1] ~ normal(0, 3);
    B[2] ~ normal(4, 3);
    B[3] ~ normal(8, 3);
    for (j in 4:J) B[j] ~ normal(12, 6);
  }
  
  range ~ gamma(27, 1.1);
  psill  ~ gamma(9, 0.6);
  sigma  ~ gamma(50, 5);
  eta   ~ multi_normal_cholesky(rep_vector(0, N), L_K);
  
  y ~ normal(mu, sigma);
  
}



