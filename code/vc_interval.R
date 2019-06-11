vc_interval <- function(model, digits=2) {
  model %>%
    confint() %>% 
    data.frame() %>%
    rownames_to_column(var = 'Component') %>% 
    filter(grepl(Component, pattern = 'sig')) %>% 
    mutate(Component = c('Person', 'Residual')) %>% 
    mutate_if(is.numeric, rnd, digits = digits) %>%
    rename(`2.5%` = X2.5.., `97.5%` = X97.5..)
}