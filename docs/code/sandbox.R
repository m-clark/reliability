# Basic congeneric ----------------------------------------------------------


library(tidyverse)
set.seed(1234)
N = 1000
N_items = 10
lambda = rep(.8, N_items)
r_congeneric = psych::sim.congeneric(loads = lambda)
# visibly::corr_heat(r_congeneric, pal = 'acton', dir = 1)
d_congeneric = 
  mvtnorm::rmvnorm(N, sigma = r_congeneric) %>% 
  as_data_frame() %>% 
  rename_all(function(x) str_replace(x, 'V', 'item_'))

d_congeneric_std = scale(d_congeneric)

glimpse(d_congeneric)
fa_congeneric = psych::fa(d_congeneric)
fa_congeneric

d_congeneric_long = d_congeneric %>% 
  rownames_to_column('subject') %>% 
  gather(key=item, value=value, -subject)

d_congeneric_long

d_congeneric_long_std = d_congeneric_std %>% 
  as_data_frame() %>% 
  rownames_to_column('subject') %>% 
  gather(key=item, value=value, -subject)

library(lme4)
mixed_congeneric = lmer(value ~ (1|subject), d_congeneric_long)
# summary(mixed_congeneric)

sjstats::icc(mixed_congeneric)

# g theory 
g_mod = gtheory::gstudy(value ~ (1|subject), data=d_congeneric_long)

# generalizability results for a single item
gtheory::dstudy(g_mod, colname.objects = 'subject')$generalizability  # same as ICC


# generalizability results for 10 items; compare to alpha; note also average r vs. ICC
gtheory::dstudy(g_mod, 
                colname.objects = 'subject', 
                colname.scores = 'value', 
                data=data.frame(d_congeneric_long))$generalizability
psych::alpha(d_congeneric)$total  # compare alpha and average r to previous; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4138888/ alpha = p*ICC/(1 + (p-1)*ICC)
psych::omega(d_congeneric, nfactors = 1)$omega.tot  # omega total

mean(fa_congeneric$loadings^2)  # compare to icc
sum(fa_congeneric$loadings)^2/var(rowSums(d_congeneric))  # From McDonald, compare to omega asymptotic?
psych::omega(d_congeneric, nfactors = 1)$omega.lim

1-sum(fa_congeneric$uniquenesses)/var(rowSums(d_congeneric)) # compare to omega tot



# Unequal items with at least one poor loading ----------------------------------------------------------

set.seed(1234)
N = 1000
lambda = c(rep(.8, 5), .7,.7,.6,.2,.2)
N_items = length(lambda)

r_congeneric = psych::sim.congeneric(loads = lambda)
# visibly::corr_heat(r_congeneric, pal = 'acton', dir = 1)
d_congeneric = 
  mvtnorm::rmvnorm(N, sigma = r_congeneric) %>% 
  as_data_frame() %>% 
  rename_all(function(x) str_replace(x, 'V', 'item_'))


# glimpse(d_congeneric)


fa_congeneric = psych::fa(d_congeneric)
fa_congeneric


# Yuan approach

library(tidyverse)
set.seed(1234)
N = 1000
N_items = 10
r_congeneric = psych::sim.congeneric(loads = rep(.8, N_items))
d_congeneric = 
  mvtnorm::rmvnorm(N, sigma = r_congeneric) %>% 
  as_data_frame() %>% 
  rename_all(function(x) str_replace(x, 'V', 'item_'))


# ke-hai 2003 using cor
cor_observed = cor(d_congeneric)
s = cor_observed[lower.tri(cor_observed, diag = T)]
init_mat = diag(N_items)
init_mat[lower.tri(init_mat, diag = T)]

a = c(init_mat[lower.tri(init_mat, diag = T)])
b = a
b[b==0] = 2
g = crossprod(a, s) / crossprod(b, s)
alpha = (1 - g) * N_items / (N_items - 1)
alpha
psych::alpha(d_congeneric)$total$std



# Bayesian approach ----------------

# This posterior predictive approach just made sense to me as a simple and
# straightforward approach, i.e. estimate corr matrix and base alpha on that,
# but Padilla & Zhang 2011 suggested.  The other is based on an asymptotic
# normal as in Zyl (2000) (also Yuan 2002, 2003)


alpha_stan = "

data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p]  X;  // standardized
}

transformed data {
  vector[p] means;
  vector[p] L_sigma;

  means = rep_vector(0, p);
  L_sigma = rep_vector(1, p);
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
  corr_matrix[p] r;
  vector[p] ev;
  real omega;
  real omega_sq;
  real alpha2;
  real theta;
  
  r = tcrossprod(r_chol);
  
  omega_sq = 2 * p^2 / ((p - 1)^2 * sum(r)^3) * 
    (sum(r) * trace(r .* r) + sum(r) * trace(r)^2 - 2 * trace(r) * sum(r .* r)) ;  
  omega = sqrt(omega_sq); 

  ev = eigenvalues_sym(r);

  alpha2 = (1 - trace(r)/sum(r)) * 1.0 * p/(p - 1);  
  theta = (1 - 1/max(ev)) * 1.0 * p/(p - 1);   // does not rely on normal, Armor 1973
}
"

library(rstan)


data_list = list(X=d_congeneric_std, N=nrow(d_congeneric_std), p=ncol(d_congeneric_std))
alpha_bayes = stan(file = 'code/alpha.stan', data = data_list, cores = 4, thin = 4)
print(alpha_bayes, par=c('alpha', 'alpha2', 'omega', 'omega_sq'), digits=4)
# launch_shinystan(alpha_bayes)

r_est = matrix(get_posterior_mean(alpha_bayes, par='r')[,5], N_items, N_items)

cronbach = get_posterior_mean(alpha_bayes, par='alpha')[,5]
cronbach
psych_alpha = psych::alpha(d_congeneric_std)$total$std
psych_alpha
psych::alpha(d_congeneric_std, n.iter = 1000)$boot.ci
psych::alpha.ci(psych_alpha, n.obs=nrow(d_congeneric_std), digits = 4)

r_resid = (cor(d_congeneric_std) - r_est)^2
sqrt(mean(r_resid[lower.tri(r_resid)]))

library(tidybayes)
alpha_draws = spread_draws(alpha_bayes, alpha, alpha2, theta)

alpha_draws %>% 
  gather(key=alpha, value=value, -.chain, -.iteration, -.draw) %>% 
  qplot(x=value, geom = 'density', color = alpha, data=.)
alpha_draws %>% mean_qi()

alpha_draws %>% 
  ggplot(aes(y = '', x = alpha)) +
  geom_halfeyeh(.width = c(.95, .8))


# bayes alpha icc ---------------------------------------------------------

# alpha as icc/generalizability statistic from g-theory

library(rstanarm)
mixed_congeneric_bayes = stan_lmer(value ~ (1|subject), d_congeneric_long_std, cores=4, thin=4)
VarCorr(mixed_congeneric_bayes)

# note that Sigma is var intercept, while sigma is residual sd
var_comp_draws = spread_draws(mixed_congeneric_bayes, Sigma[subject:(Intercept),(Intercept)], sigma) %>% 
  mutate(sigma = sigma^2,
         alpha = Sigma/(Sigma + sigma/N_items))

var_comp_draws %>% mean_qi(alpha)
print(alpha_bayes, par='alpha', digits=4)

alpha_boot = psych::alpha(d_congeneric_std, n.iter = 1000)$boot %>% 
  as_data_frame() %>% 
  rename(alpha = std.alpha) %>% 
  select(alpha) %>% 
  mutate(.chain = 1, 
         .iter = 1:1000,
         .draw = .iter,
         Group = 'alpha_boot')


alpha_draws_for_plot = alpha_draws %>% 
  mutate(Group = 'alpha_bayes_r') %>% 
  select(-alpha2, -theta) 
alpha2_draws_for_plot = alpha_draws %>% 
  mutate(Group = 'alpha2_bayes_r') %>% 
  select(-alpha) %>% 
  rename(alpha = alpha2)

# interesting, is pp_alpha bias or regularization? With only four items, all had same mean, but pp_alpha was still sharper
dist_dat = var_comp_draws %>% 
  ungroup() %>% 
  select(.chain, .iteration, .draw, alpha) %>% 
  mutate(Group = 'alpha_bayes_icc') %>% 
  bind_rows(alpha_draws_for_plot, alpha2_draws_for_plot) %>% 
  bind_rows(alpha_boot)

dist_means = dist_dat %>% tidyext::num_by(alpha, Group, digits=5)

ggplot(dist_dat, aes(x=alpha, color=Group, fill=Group)) +
  geom_point(x = psych::alpha(d_congeneric)$total$std, y=0, size=10, color='gray50', alpha=.25) +
  geom_density(alpha = .25) +
  geom_point(aes(x = Mean, y = 0), alpha=.5, size=3, data = dist_means) +
  scico::scale_color_scico_d() +
  scico::scale_fill_scico_d() +
  visibly::theme_trueMinimal()

# bayes normal approx vs. bootstrap (see Yuan 2003)
qqplot(alpha_draws_for_plot$alpha, alpha_boot$alpha)
abline(coef(lm(sort(alpha_boot$alpha) ~ sort(alpha_draws_for_plot$alpha))))
