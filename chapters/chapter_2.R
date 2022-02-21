# Libraries ---------------------------------------------------------------
pkgs <- c("R2jags", "rstan", "bayesplot")
sapply(pkgs, require, character.only = TRUE)
options(mc.cores = parallel::detectCores())


# The Core Model ----------------------------------------------------------
# This chapter provides an overview of the foundations of the standard NMA model, 
# specifically the binomial model with a logit link function. 

# Development of the Core Model -------------------------------------------
## Pairwise Meta-Analysis --------------------------------------------------
### Fixed Effects -----------------------------------------------------------
# N studies:
n_s <- 11
# N trial arms:
n_arms <- 2
# Events:
r <- matrix(c(
 3, 10, 40, 5, 5, 2, 19, 59, 5, 16, 8,
 1,  3, 32, 5, 3, 3, 20, 52, 2, 12, 6),
 ncol = 2)
# Observations:
n <- matrix(c(
 55, 94, 573, 75, 69, 61, 419, 782, 81, 226, 66,
 55, 95, 565, 75, 71, 62, 421, 790, 81, 225, 71), 
 ncol = 2)
# Data list:
data_list <- list(n_s = n_s, n_arms = n_arms, r = r, n = n)

#### Jags model --------------------------------------------------------------
# model
Ch2_FE_Bi_logit_pair <- "# Binomial likelihood, logit link
# Pairwise meta-analysis
# Fixed effects model
model {                    # *** PROGRAMME STARTS
 for (i in 1:n_s) {        # Loop through studies
  mu[i] ~ dnorm(0, 0.0001) # vague prior for trial baselines
  for (k in 1:n_arms) {
   r[i, k] ~ dbinom(p[i, k], n[i, k]) # binomial likelihood
   logit(p[i, k]) <- mu[i] + d[k]     # model for linear predictor
  }
 }
 d[1] <- 0                 # treatment effect is zero for reference treatment
 d[2:n_arms] ~ dnorm(0, 0.0001)   # vague prior for treatment effect
 OR <- exp(d[2])
 prob_harm <- step(d[2])
 
}                          # *** PROGRAMME ENDS

"
writeLines(text = Ch2_FE_Bi_logit_pair, con = "jags/Ch2_FE_Bi_logit_pairs.txt")

jags_Model <- jags(model.file = "jags/Ch2_FE_Bi_logit_pairs.txt", 
                   data = data_list, parameters.to.save = c(
                    "d", "OR", "prob_harm"))
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_FE_Bi_logit_pair <- stan(file = "stan/Ch2_FE_Bi_logit_pair.stan",
                              data = data_list, chains = 4,
                              pars = c("d", "OR", "prob_harm"))
print(stan_FE_Bi_logit_pair)
sims_stan_FE_Bi_logit_pair <- extract(stan_FE_Bi_logit_pair)
hist(sims_stan_FE_Bi_logit_pair$d[, 2])

# Inspection:
mcmc_dens(stan_FE_Bi_logit_pair, pars = "d[2]")
mcmc_trace(stan_FE_Bi_logit_pair, pars = "d[2]")

# Success! Both JAGS and Stan are giving similar output...

### Random effects ----------------------------------------------------------
# In a random effects model, each study i provides an estimate of the 
# study-specific treatment effects \Delta_{i, 1, 2} - that is, the relative
# effect of treatment 2 compared to treatment 1 on some scale in trial i. 
# Hence, treatment effects are not assumed to be 'equal' but rather 
# exchangeable or 'similar'. Hence, the effects can be derived from a 'common'
# distribution.

#### Jags model --------------------------------------------------------------
# model
Ch2_RE_Bi_logit_pair <- "# Binomial likelihood, logit link
# Pairwise meta-analysis
# Random effects model
model {                    # *** PROGRAMME STARTS
 for (i in 1:n_s) {        # LOOP THROUGH STUDIES
 delta[i, 1] <- 0          # treatment effect 0 for control arm
  mu[i] ~ dnorm(0, 0.0001) # vague prior for trial baselines
  delta[i, 2] ~ dnorm(d[2], tau) # trial-specific LOR distributions
  
  for (k in 1:n_arms) {    # LOOP THROUGH ARMS
   r[i, k] ~ dbinom(p[i, k], n[i, k]) # binomial likelihood
   logit(p[i, k]) <- mu[i] + d[k] + delta[i, k]  # model for linear predictor
  }
 }
 d[1] <- 0                 # treatment effect is zero for reference treatment
 d[2:n_arms] ~ dnorm(0, 0.0001) # vague prior for treatment effect
 tau <- pow(sd, -2)        # between-trial precision
 sd ~ dunif(0, 2)          # prior for between-trial stdev
 OR <- exp(d[2])           # LOR to OR
 prob_harm <- step(d[2])   # OR to p
 
}                          # *** PROGRAMME ENDS

"
writeLines(text = Ch2_RE_Bi_logit_pair, con = "jags/Ch2_RE_Bi_logit_pairs.txt")

jags_Model <- jags(model.file = "jags/Ch2_RE_Bi_logit_pairs.txt", 
                   data = data_list, parameters.to.save = c(
                    "d", "OR", "prob_harm"))
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_FE_Bi_logit_pair <- stan(file = "stan/Ch2_RE_Bi_logit_pair.stan",
                              data = data_list, chains = 4,
                              pars = c("d", "OR", "prob_harm"))
print(stan_FE_Bi_logit_pair)
sims_stan_FE_Bi_logit_pair <- extract(stan_FE_Bi_logit_pair)
hist(sims_stan_FE_Bi_logit_pair$d[, 2])

# Inspection:
mcmc_dens(stan_FE_Bi_logit_pair, pars = "d[2]")
mcmc_trace(stan_FE_Bi_logit_pair, pars = "d[2]")

# Success! Both JAGS and Stan are giving similar output...

## Extension to Indirect Comparisons and Network Meta-Analysis -------------