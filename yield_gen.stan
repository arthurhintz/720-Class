data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> J;
  array[N] vector[2] coord;
  array[N] int<lower=1,upper=K> trt;
  array[N] int<lower=1,upper=J> blocks;
  vector[N] ec;
}
  
generated quantities {
  
  real a;
  vector[K] T;
  vector[J] B;
  real<lower=0> sigma;
  vector[N] eta;
  vector[N] mu;
  vector[N] y_rep;
  
  a = normal_rng(1, 5);
  
  // Treatmetns
  T[1] = normal_rng(-4, 2);
  for (k in 2:5) T[k] = normal_rng(-2, 2);
  for (k in 6:9) T[k] = normal_rng(0, 2);
  for (k in 10:12) T[k] = normal_rng(2, 2);
  for (k in 13:K) T[k] = normal_rng(4, 2);
  
  
  // Blocks
 if (J == 1) {
    B[1] = normal_rng(0, 0.0001); 
  } else {
    B[1] = normal_rng(0, 3);
    B[2] = normal_rng(4, 3);
    B[3] = normal_rng(8, 3);
    for (j in 4:J) B[j] = normal_rng(12, 6);
  }
  sigma  = gamma_rng(50, 10);


  for (i in 1:N)
    mu[i] = T[trt[i]] + B[blocks[i]] + ec[i] * a;
    

  for (n in 1:N)
    y_rep[n] = normal_rng(mu[n], sigma);
}



