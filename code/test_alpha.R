# test alpha
neuroticism_no_na = na.omit(scale(neuroticism))

data_list = list(X=scale(neuroticism_no_na), 
                 N=nrow(neuroticism_no_na), 
                 p=ncol(neuroticism_no_na),
                 scale_flag = 0)

library(rstan)

alpha_bayes = stan(file = 'code/alpha.stan', 
                   data = data_list, 
                   cores = 4,
                   # chains = 1,
                   # iter = 10 #,
                   thin = 4,
                   control = list(adapt_delta = .99)
)

print(alpha_bayes, digits=3, pars = c('alpha', 'alpha2'))
launch_shinystan(alpha_bayes)

library(tidybayes)
alpha_draws = spread_draws(alpha_bayes, alpha2)

head(alpha_draws)

mean_qi(alpha_draws)

library(gganimate)

n_hops = 100

alpha_animate = alpha_draws %>% 
  ggplot(aes(x = alpha2, y= 'alpha', group=.draw)) +
  # ggplot(aes(x = alpha2, y= 'alpha')) +
  geom_point(x = .813, size=10, color='#ff5500', alpha=.25) +
  geom_point(alpha=.05, size=3, position = position_jitter()) +
  # transition_states(.draw, transition_length = 1, state_length = 1) +
  transition_reveal(.draw) + # group .draw to keep dots persistent
  theme_trueMinimal()
    

# animate(alpha_animate, nframes = n_hops*2)
animate(alpha_animate, nframes = n_hops*1)
anim_save('img/alpha_unc.gif')



alpha_animate = alpha_draws %>% 
  mutate(alpha = 'alpha') %>% 
  ggplot(aes(y = alpha2, x= alpha, group=.draw)) +
  # ggplot(aes(x = alpha2, y= 'alpha')) +
  geom_point(y = .813, size=10, color='#ff5500', alpha=.25) +
  geom_sina(aes(), bins=100) +
  # transition_states(.draw, transition_length = 1, state_length = 1) +
  transition_reveal(.draw) + # group .draw to keep dots persistent
  theme_trueMinimal()


# animate(alpha_animate, nframes = n_hops*2)
animate(alpha_animate, nframes = n_hops*1)
anim_save('img/alpha_unc.gif')



alpha(scale(neuroticism_no_na), digits = 2)
MBESS::ci.reliability(scale(neuroticism_no_na), type='alpha', interval.type = 'ml')

# for ave
library(lavaan)
neuro_one_fac = "
 neuro =~ N1 + N2 + N3 + N4 + N5
"
semTools::reliability(cfa(neuro_one_fac, data = neuroticism_no_na))

