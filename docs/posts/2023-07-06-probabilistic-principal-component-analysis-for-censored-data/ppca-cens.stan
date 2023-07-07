// modified from: https://github.com/nbip/PPCA-Stan/blob/master/PPCA.ipynb

data {
  int<lower=1> N; // num datapoints
  int<lower=1> D; // num data-space dimensions
  int<lower=1> K; // num latent space dimensions
  array[D, N] real y;
  // data for censored imputation:
  int<lower=0> Ncens_y1; // number of censored, y1
  int<lower=0> Ncens_y2; // number of censored, y2
  array[Ncens_y1] int<lower=1> Jcens_y1; // positions of censored, y1
  array[Ncens_y2] int<lower=1> Jcens_y2; // positions of censored, y2
  real U_y1; // left-censoring limit, y1
  real U_y2; // left-censoring limit, y2
}
transformed data {
  matrix[K, K] Sigma; // identity matrix
  vector<lower=0>[K] diag_elem;
  vector<lower=0>[K] zr_vec; // zero vector
  for (k in 1 : K) {
    zr_vec[k] = 0;
  }
  for (k in 1 : K) {
    diag_elem[k] = 1;
  }
  Sigma = diag_matrix(diag_elem);
}
parameters {
  matrix[D, K] A; // transformation matrix / PC loadings
  array[N] vector[K] x; // latent variables
  real<lower=0> sigma; // noise variance
  // censored value parameters:
  array[Ncens_y1] real<upper=U_y1> Ycens_y1; // estimated censored, y1
  array[Ncens_y2] real<upper=U_y2> Ycens_y2; // estimated censored, y2
}
transformed parameters {
  // combine observed with estimated censored:
  array[D, N] real yl = y;
  yl[1, Jcens_y1] = Ycens_y1;
  yl[2, Jcens_y2] = Ycens_y2;
}
model {
  for (i in 1 : N) {
    x[i] ~ multi_normal(zr_vec, Sigma);
  } // zero-mean, identity matrix
  for (i in 1 : N) {
    for (d in 1 : D) {
      //y[d,i] ~ normal(dot_product(row(A, d), x[i]), sigma);
      target += normal_lpdf(yl[d, i] | dot_product(row(A, d), x[i]), sigma);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1 : N) {
    for (d in 1 : D) {
      log_lik[n] = normal_lpdf(yl[d, n] | dot_product(row(A, d), x[n]), sigma);
    }
  }
}


