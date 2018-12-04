data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p]  X;
  int<lower=0, upper=1> scale_flag;   // use standardized or not
}

transformed data {
  vector[p] means;
  vector[p] L_sigma;

  if (scale_flag) {
    means = rep_vector(0, p);
    L_sigma = rep_vector(1, p);
  } else {
    for (i in 1:p) {
      means[i] = mean(X[,i]);
      L_sigma[i] = sd(X[,i]);
    }
  }
}

parameters {
  real<lower=0, upper=1> alpha;
  cholesky_factor_corr[p] r_chol;
}

model {
  matrix[p, p] sigma_chol;
  matrix[p, p] r;
  real omega_sq;
  real omega;
  real a0;

  r_chol ~ lkj_corr_cholesky(.5);  // less prob at identity matrix
  // r_chol ~ lkj_corr_cholesky(1);  // agnostic

  sigma_chol = diag_pre_multiply(L_sigma, r_chol);

  for (n in 1:N)
    X[n] ~ multi_normal_cholesky(means, sigma_chol);
    
  r = tcrossprod(r_chol);    
  
  // note that p^2/(p-1)^2 in some depictions is already in the calculation here
  omega_sq = 2 * p^2 / ((p - 1)^2 * sum(r)^3) * 
    (sum(r) * trace(r .* r) + sum(r) * trace(r)^2 - 2 * trace(r) * sum(r .* r)) ;  
  omega = sqrt(omega_sq); 
  
  a0 = (1 - trace(r)/sum(r)) * 1.0 * p/(p - 1);  // 1.0 to avoid integer division
  
  alpha ~ normal(a0, omega/sqrt(N));   // following Zyl, Yuan
}

generated quantities {
  vector[p] ev;
  corr_matrix[p] r_std;
  cov_matrix[p] r_cov;
  matrix[p,p] r;
  real omega;
  real omega_sq;
  real alpha2;
  real theta;
  
  if (scale_flag) {
    r_std = tcrossprod(r_chol);
    r = r_std;
  } else {
    r_std = tcrossprod(r_chol);
    r_cov = quad_form_diag(r_std, L_sigma);
    r = r_cov;
  }
  
  omega_sq = 2 * p^2 / ((p - 1)^2 * sum(r)^3) * 
    (sum(r) * trace(r .* r) + sum(r) * trace(r)^2 - 2 * trace(r) * sum(r .* r)) ;  
  omega = sqrt(omega_sq); 

  ev = eigenvalues_sym(r);

  alpha2 = (1 - trace(r)/sum(r)) * 1.0 * p/(p - 1);  
  theta = (1 - 1/max(ev)) * 1.0 * p/(p - 1);   // does not rely on normal, Armor 1973
  
}
