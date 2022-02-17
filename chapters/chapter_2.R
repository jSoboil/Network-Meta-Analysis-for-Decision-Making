# Libraries ---------------------------------------------------------------
pkgs <- c("R2jags", "rstan")
sapply(pkgs, require, character.only = TRUE)
options(mc.cores = parallel::detectCores())


# The Core Model ----------------------------------------------------------
# This chapter provides an overview of the foundations of the standard NMA model, 
# specifically the binomial model with a logit link function. 


# Development of the Core Model -------------------------------------------
# Skipped due to prior knowledge

## Incorporating Multi-arm trials ------------------------------------------
# Fixed effects Thrombolytics example section 2.3.2:

# Number of arms:
n_arms <- c(rep(2, 36))
# Treatment in arm [, k] of trial [i, ]:
t <- matrix(c(
 t_1 = c(rep(1, 19), rep(2, 3), rep(3, 14)),
 t_2 = c(3, rep(2, 8), 4, 5, rep(7, 11), 5, 5, 6, rep(7, 11))
 ),
 ncol = 2)
# Numerators and denominators arm [, k] of trial [i, ]:
# Arm 1 for trial i:

r <- matrix(c(
 r_1 = c(1472, 3, 12, 7, 10, 887, 5, 1455, 9, 4, 285, 11, 1, 8, 1, 4, 14, 9, 42,
         2, 13, 2, 13, 356, 522, 3, 10, 40, 5, 5, 2, 19, 59, 5, 16, 8),
 r_2 = c(652, 3, 7, 4, 5, 929, 2, 1418, 6, 6, 270, 2, 3, 5, 1, 0, 7, 3, 29, 3, 5,
         2, 7, 757, 523, 1, 3, 32, 5, 3, 3, 20, 52, 2, 12, 6)), ncol = 2)
n <- matrix(data = c(
 n_1 = c(20251, 65, 159, 85, 135, 10396, 63, 13780, 130, 107, 3004, 149, 50, 58,
         53, 45, 99, 41, 421, 44, 200, 56, 155, 4921, 8488, 55, 94, 573, 75, 69,
         61, 419, 782, 81, 226, 66),
 n_2 = c(10396, 64, 157, 86, 135, 10372, 59, 13746, 123, 109, 3006, 152, 50, 54,
         47, 42, 101, 46, 429, 46, 195, 47, 169, 10138, 8461, 55, 95, 565, 75,
         71, 62, 421, 790, 81, 225, 71)), ncol = 2)
n_s <- 36
n_t <- length(unique(t[[2]]))

# data list for Stan:
data_list <- list(n_arms = n_arms, n_s = n_s, n_t = n_t, t = t, r = r, n = n)
data_list

# jags model:
string_FE_Bi_logit <- "# Binomial Likelihood, logit link
# Fixed effect model
model {                                # *** PROGRAMME STARTS
  for(i in 1:n_s){                      # Loop through studies
    mu[i] ~ dnorm(0,.0001)              # vague priors for all trial baselines
      for (k in 1:n_arms[i])  {          # Loop through arms
        r[i,k] ~ dbin(p[i, k],n[i, k])    # binomial likelihood
        logit(p[i, k]) <- mu[i] + d[t[i, k]] - d[t[i, 1]] # model for linear predictor
      }
    }
    
 d[1] <- 0                          # treatment effect zero for reference treatment
 for (k in 1:n_t) {
   d[k] ~ dnorm(0, 0.0001)              # vague priors for treatment effects
 }
}                                       # *** PROGRAMME ENDS

"

writeLines(text = string_FE_Bi_logit, con = "jags/jags_FE_Bi_logit.txt")

# JAGS --------------------------------------------------------------------
jags_FE_Bi_logit <- jags(model.file = "jags/jags_FE_Bi_logit.txt",
                         data = data_list, parameters.to.save = "d")

# Stan --------------------------------------------------------------------
stan(file = "stan/Ch2_FE_Bi_logit.stan", data = data_list)
