extract_vc <- function(model) {
  data.frame(lme4::VarCorr(model)) %>% 
    mutate_if(function(x) all(is.na(x)), function(x) x = NULL)
}