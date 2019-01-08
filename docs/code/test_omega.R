# test omega
neuroticism_no_na = na.omit(scale(neuroticism))

data_list = list(X=scale(neuroticism_no_na), 
                 N=nrow(neuroticism_no_na), 
                 p=ncol(neuroticism_no_na),
                 scale_flag = 0)

library(rstan)

omega_bayes = stan(file = 'code/omega.stan', 
                   data = data_list, 
                   cores = 4,
                   # chains = 1,
                   # iter = 10 #,
                   thin = 4,
                   control = list(adapt_delta = .99)
                   )

print(omega_bayes, digits=3, pars = c('omega', 'rho', 'ave', 'lambda'))
launch_shinystan(omega_bayes)

fa(scale(neuroticism_no_na), fm = 'ml')
omega(scale(neuroticism_no_na), digits = 2)
MBESS::ci.reliability(scale(neuroticism_no_na), interval.type = 'ml')

# for ave
library(lavaan)
neuro_one_fac = "
 neuro =~ N1 + N2 + N3 + N4 + N5
"
semTools::reliability(cfa(neuro_one_fac, data = neuroticism_no_na))

