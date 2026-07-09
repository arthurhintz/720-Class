data {
  int<lower=1> N;
  int<lower=1> K;
  //int<lower=1> J;
  array[N] vector[2] coord;
  //array[N] int<lower=1, upper=J> block_id;
}
  
generated quantities {
  
  real mu_0;
  //vector[J] B;
  real<lower=0> range;
  real<lower=0> psill;
  real<lower=0> sigma;
  vector[N] eta;
  vector[N] mu;
  vector[N] ec_rep;
  
  
  //Intercept
  mu_0 = normal_rng(0, 10);
  
  //Blocks
  // B[1] = normal_rng(0, 4);
  // B[2] = normal_rng(4, 4);
  // B[3] = normal_rng(8, 4);
  // B[4] = normal_rng(12, 4);
  
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

  for (n in 1:N)
    mu[n] = mu_0 + eta[n]; // B[block_id[n]]

  for (n in 1:N)
    ec_rep[n] = normal_rng(mu[n], sigma);
}

