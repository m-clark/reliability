# From NAJAFABADI

Alpha <- function(all_data) {
  library(psych)
  n = dim(all_data)[1]
  p = dim(all_data)[2]
  D = c()
  Q = c()
  B = c()
  B2 = c()
  
  if (p <= 2) {
    print("There is no solution for dimension less than 3")
  } else {
    # Cronbach’s alpha calculation
    Alpha_Cronbach <- function(data) {
      k = dim(data)[2]
      
      alpha = k / (k - 1) * (1 - tr(cov(data, use = "na.or.complete")) /
                               sum(cov(data, use = "na.or.complete")))
      return(alpha)
    }
    
    # Ordinal’s theta calculation
    theta_ordinal = function(data) {
      k = dim(data)[2]
      
      theta = k / (k - 1) * (1 - 1 / max(eigen(cor(data, use = "na.or.complete"))$value))
      return(theta)
    }
    
    # Bayes estimator
    Bayes <- function(data) {
      n = dim(data)[1]
      k = dim(data)[2]
      
      r = k / (k - 1) * (1 - tr(mixed.cor(data)$rho) / sum(mixed.cor(data)$rho))   # removed deprecated polycor=T
      theta = theta_ordinal (data)
      a = (1 - theta)^2 / 100
      b = (1 - theta) / 100
      
      f1 = function(x) {
        x * (1 - x)^(min(30, a + (n - 3) / 2)) *
          exp(-b * (1 - x)) * 
          (1 + (1 - x) / (k - r * (k - 1)))^(-(n - 1) * k / 2)
      }
      
      f2 = function(x) {
        (1 - x)^(min(30, a + (n - 3) / 2)) * 
          exp(-b * (1 - x)) * 
          (1 + (1 - x) / (k - r * (k - 1)))^(-(n - 1) * k / 2)
      }
      
      bayes = integrate(f1, lower = -Inf, upper = 1)$value / 
        (0.e-20 + integrate(f2, lower = -Inf, upper = 1)$value)
      
      return(bayes)
    }
    
    Bayes_LINEX <- function(data) {
      n = dim(data)[1]
      k = dim(data)[2]
      
      r = k / (k - 1) * (1 - tr(mixed.cor(data)$rho) / sum(mixed.cor(data)$rho))   # removed deprecated polycor=T
      theta = theta_ordinal (data)
      a = (1 - theta)^2 / 100
      b = (1 - theta) / 100
      c = -10
      
      f1 = function(x) {
        exp(-c * x) * (1 - x)^(min(30, a + (n - 3) / 2)) * 
          exp(-b * (1 - x)) *
          (1 + (1 - x) / (k - r * (k - 1)))^(-(n - 1) * k / 2)
      }
      
      f2 = function(x) {
        (1 - x)^(min(30, a + (n - 3) / 2)) * 
          exp(-b * (1 - x)) * 
          (1 + (1 - x) / (k - r * (k - 1))) ^ (-(n - 1) * k / 2)
      }
      
      bayes = -log(integrate(f1, lower = -Inf, upper = 1)$value / 
                     (0.e-20 + integrate(f2, lower = -Inf, upper = 1)$value)) / c
      
      return(bayes)
    }
    
    for (j in 1:p) {
      D <- c(D, Alpha_Cronbach(all_data[-j]))
      Q <- c(Q, theta_ordinal(all_data[-j]))
      B <- c(B, Bayes(all_data[-j]))
      B2 <- c(B2, Bayes_LINEX (all_data[-j]))
    }
    
    list(
      "Alpha if Item Deleted" = data.frame(
        "." = D,
        row.names = paste("Alpha Without Item", 1:p, "=")
      ) ,
      "Cronbach’s Alpha for all Items =" = Alpha_Cronbach(all_data),
      "Ordinal Theta if a Item Deleted" = data.frame(
        "." = Q,
        row.names = paste("Theta Without Item", 1:p, "=")
        ) ,
      "Ordinal Theta for all Item=" = theta_ordinal(all_data),
      "Bayesian-Squared-Error Alpha if Item Deleted" = data.frame(
        "." = B,
        row.names = paste("Bayesian-Squared-Error Alpha Without Item", 1:p, "=")
      ) , 
      "Bayesian-Squared-Error Alpha for all Item=" = Bayes(all_data), 
      "Bayesian-LINEX Alpha if a Item Deleted" = data.frame(
        "." = B2,
        row.names =
          paste("Bayesian-LINEX Alpha Without Item", 1:p, "=")
      ) , 
      "Bayesian-LINEX Alpha for all Item=" = Bayes_LINEX (all_data)
    )
  }
}



test = Alpha(d_congeneric)
test[[2]]
print(alpha_bayes, par='alpha', digits=4)
