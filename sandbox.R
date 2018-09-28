library(tidyverse)
set.seed(1234)
N = 1000
r_congeneric = psych::sim.congeneric(loads = rep(.8, 10))
# visibly::corr_heat(r_congeneric, pal = 'acton', dir = 1)
d_congeneric = 
  mvtnorm::rmvnorm(N, sigma = r_congeneric) %>% 
  as_data_frame() %>% 
  rename_all(function(x) str_replace(x, 'V', 'item_'))


glimpse(d_congeneric)
fa_congeneric = psych::fa(d_congeneric)
fa_congeneric

d_congeneric_long = d_congeneric %>% 
  rownames_to_column('subject') %>% 
  gather(key=item, value=value, -subject)

d_congeneric_long

library(lme4)
mixed_congeneric = lmer(value ~ (1|subject), d_congeneric_long)
# summary(mixed_congeneric)

sjstats::icc(mixed_congeneric)

# g theory 
g_mod = gtheory::gstudy(value ~ (1|subject), data=d_congeneric_long)

# generalizability results for a single item
gtheory::dstudy(g_mod, colname.objects = 'subject')$generalizability  # same as ICC


# generalizability results for 10 items; compare to alpha
gtheory::dstudy(g_mod, 
                colname.objects = 'subject', 
                colname.scores = 'value', 
                data=data.frame(d_congeneric_long))$generalizability
psych::alpha(d_congeneric)$total  # compare alpha and average r to previous
psych::omega(d_congeneric, nfactors = 1)$omega.tot  # omega total

mean(fa_congeneric$loadings^2)  # compare to icc
sum(fa_congeneric$loadings)^2/var(rowSums(d_congeneric))  # From McDonald, compare to omega asymptotic?
psych::omega(d_congeneric, nfactors = 1)$omega.lim

1-sum(fa_congeneric$uniquenesses)/var(rowSums(d_congeneric)) # compare to omega tot
