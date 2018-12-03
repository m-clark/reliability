library(tidyverse)
library(psych)

Ratings = data_frame(
  J1 = 1:6,
  J2 = J1,
  J3 = 6:11,
  J4 = seq(2, 12, 2),
  J5 = c(3,1,5,2,6,4),
  J6 = c(6,2,10,4,12,8)
)

n_judges = n_subjects = 6

Ratings_long = Ratings %>% 
  gather(key = 'Judge', value = Rating) %>% 
  mutate(Subject = sapply(1:n_subjects, paste0, 'S'))

cor(Ratings) # check

ICC_result = ICC(Ratings)   # by default uses lmer
print(ICC_result, all=T)

library(lme4)

# Judge = items, id = Subject in ICC_result
mod = lmer(Rating ~  (1|Subject) + (1|Judge), Ratings_long)
summary(mod)

# extract varcomp
vc = VarCorr(mod) %>% 
  data.frame() %>% 
  pull(vcov)

vc[1]/sum(vc)                                     # basic ICC, aka ICC2 single_random
vc[1]/sum(vc[1], vc[3])                           # single fixed, ICC2k
vc[1]/sum(vc[1], vc[2]/n_judges, vc[3]/n_judges)  # avg random, ICC2k 
vc[1]/sum(vc[1], vc[3]/n_judges)                  # avg fixed, ICC3k

# subject icc is ICC2 single_random_raters
lmer_icc2 = sjstats::icc(mod)
lmer_icc2

# ICC 3 = fixed raters
lmer_icc3 = sjstats::icc(lmer(Rating ~  (1|Subject) + Judge, Ratings_long))
lmer_icc3


library(gtheory)
g_result = gstudy(Rating ~  (1|Subject) + (1|Judge), data=Ratings_long)
g_result

# d-study based on single observation (raw lme4 output)
d_result_1 = dstudy(g_result, 
                    data=Ratings_long, 
                    colname.objects = 'Subject')

# dependability = ICC2 single_random_raters
# generalizability = ICC3 single_fixed_raters
d_result_1  


# d-study based on average observation (over the 6)
d_result_avg = dstudy(g_result, 
                      data=data.frame(Ratings_long), 
                      colname.objects = 'Subject',
                      colname.scores = 'Rating')

# dependability = ICC2k average_random_raters
# generalizability = ICC3k average_fixed_raters = raw_alpha
d_result_avg

alpha_result = psych::alpha(Ratings)$total
alpha_result

data_frame(
  Type = ICC_result$results$type[-c(1,4)],
  Name = rownames(ICC_result$results)[-c(1,4)],
  Gtheory = rep(c('Dependability', 'Generalizability'), 2),
  ICC = ICC_result$results$ICC[-c(1,4)],
  `Mixed model` = c(lmer_icc2[1], lmer_icc3, NA, NA),
  dstudy_1 =  c(d_result_1$dependability, d_result_1$generalizability, NA, NA),
  dstudy_avg =  c(NA, NA, d_result_avg$dependability, d_result_avg$generalizability),
  alpha = c(rep(NA, 3), alpha_result$raw_alpha)
) %>% 
  kableExtra::kable(digits = 3)


