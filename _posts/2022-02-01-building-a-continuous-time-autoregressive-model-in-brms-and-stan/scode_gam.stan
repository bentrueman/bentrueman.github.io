// generated with brms 2.16.2
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector[N] s;  // CAR(1) exponent
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline s(x)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data needed for ARMA correlations
  int<lower=0> Kar;  // AR order
  int<lower=0> Kma;  // MA order
  // number of lags per observation
  int<lower=0> J_lag[N];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int max_lag = max(Kar, Kma);
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline s(x)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  vector<lower=0>[Kar] ar;  // autoregressive coefficients, modified for CAR(1)
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
}
model {
  // likelihood including constants
  if (!prior_only) {
    // matrix storing lagged residuals
    matrix[N, max_lag] Err = rep_matrix(0, N, max_lag);
    vector[N] err;  // actual residuals
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N) + Xs * bs + Zs_1_1 * s_1_1;
    // include ARMA terms
    for (n in 1:N) {
      err[n] = Y[n] - mu[n];
      for (i in 1:J_lag[n]) {
        Err[n + 1, i] = err[n + 1 - i];
      }
      mu[n] += Err[n, 1:Kar] * pow(ar, s[n]); // CAR(1)
    }
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 0.1, 2.5);
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
  target += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}

