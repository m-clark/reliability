// to do: optimize for standardized, general cleanup

data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N, p]  X;
  int<lower=0, upper=1> scale_flag;   // use standardized or not
}

transformed data {
  vector[p] means;
  vector[p] sds;
  vector[p] L_sigma;
  vector[N] sum_score;

  if (scale_flag) {
    means = rep_vector(0, p);
    for (i in 1:p) 
      sds[i] = sd(X[,i]);
    L_sigma = rep_vector(1, p);
  } else {
    for (i in 1:p) {
      means[i] = mean(X[,i]);
      L_sigma[i] = sd(X[,i]);
    }
  }
  
  for (n in 1:N) sum_score[n] = sum(X[n]);
  
}

parameters {
  vector<lower=0, upper=1>[p] lambda;  // loadings
  vector<lower=0, upper=1>[p] u_var;   // uniquenesses
  vector[scale_flag ? 0 : p] b;        // intercepts, unless scale_flag
  vector[N] Z;                         // the latent variable
}

transformed parameters {
  real<lower=0, upper=1> omega;
  real<lower=0, upper=1> ave;
  vector<lower=0, upper=1>[p] lambda_sq;

  for (i in 1:p) 
    lambda_sq[i] = lambda[i]^2;
  omega = sum(lambda)^2 / (sum(lambda)^2 + sum(u_var));
  ave = mean(lambda_sq);
  
}

model {
  if (!scale_flag)
    b ~ normal(0, sd(sum_score));
  Z ~ normal(0, 1);
  
  if (!scale_flag) {
    for (i in 1:p)
      X[, i] ~ normal(b[i] + Z*lambda[i], sqrt(u_var[i]));    
  } 
  else {
    for (i in 1:p)
      X[, i] ~ normal(Z*lambda[i], sqrt(u_var[i]));
  }


}

generated quantities {
 real<lower=-1, upper=1> rho;   // correlation of sum score and latent
 
 rho = dot_product(sum_score - mean(sum_score), Z - mean(Z)) / 
  (sd(sum_score)*sd(Z)*(N-1));
  
}
