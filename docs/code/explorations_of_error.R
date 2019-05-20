# reliability of specific effects?
library(tidyverse)
set.seed(1234)
group_N = 1:50
group = rep(group_N, group_N)
x = rnorm(sum(group_N))
re = rnorm(50, sd=.25)

y = .5*x + re[group] + rnorm(sum(group_N), sd=.25)

df = data.frame(x, y, group)

library(lme4)

mod = lmer(y ~ x + (1|group), df)
summary(mod)

re_est = data.frame(ranef(mod, condVar=T)) %>% 
  mutate(rel_sem = sigma(mod)^2/1:50,                          # relative standard error of measurement
         denom = (c(VarCorr(mod)[[1]]) + sigma(mod)^2)/1:50)   # style of denominator in d-study

qplot(group_N, condsd^2, data=re_est) +
  geom_line(aes(y=rel_sem))  +
  geom_line(aes(y=denom), color = 'darkred') 

# 56 in lme4 article (ignoring weights)
# Lambdat = transpose of relative covariance matrix
# Zt = transpse of RE matrix
V = solve(
  mod@pp$Lambdat %*% mod@pp$Zt %*% t(mod@pp$Zt) %*% t(mod@pp$Lambdat) + diag(50)
)

# 58
cond_sd = diag(sigma(mod)^2 * mod@pp$Lambdat %*% V %*% t(mod@pp$Lambdat))

all.equal(cond_sd, re_est$condsd^2)  

condvar_dat = 
  data.frame(ranef(mod, condVar=T)) %>% 
  mutate(
    condvar = condsd^2,
    V_ = diag(V),
    LVL = diag(crossprod(mod@pp$Lambdat) %*% V)
  )

condvar_dat %>% mutate(grp_N = group_N) %>% select_if(is.numeric) %>% cor


var_grp = sigma(mod)^2 * crossprod(mod@pp$Lambdat)  # note that latter is diagonal(ngrps)

var_grp[1]
VarCorr(mod) %>% data.frame()


re = V %*% mod@pp$Lambdat %*% mod@pp$Zt %*% (mod@resp$y - predict(mod, re.form=NA))
data.frame(re[,1], ranef(mod)[[1]])
