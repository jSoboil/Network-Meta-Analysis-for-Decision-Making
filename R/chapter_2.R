# Libraries ---------------------------------------------------------------
# required packages for analysis
pkgs <- c("R2jags", "rstan", "rstanarm", "bayesplot", "dplyr", 
          "tidyr", "igraph")
# install packages if not installed
installed_packages <- pkgs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
 install.packages(pkgs[!installed_packages])
}
# load required packages
invisible(lapply(pkgs, library, character.only = TRUE))
# set number of parallel cores
options(mc.cores = parallel::detectCores())

# The Core Models ---------------------------------------------------------
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
writeLines(text = Ch2_FE_Bi_logit_pair, 
           con = "jags/02_JAGS/Ch2_FE_Bi_logit_pair.txt")

jags_Model <- jags(model.file = "jags/02_JAGS/Ch2_FE_Bi_logit_pair.txt", 
                   data = data_list, parameters.to.save = c(
                    "d", "OR", "prob_harm"))
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_FE_Bi_logit_pair <- stan(file = "stan/02_Stan/Ch2_FE_Bi_logit_pair.stan",
                              data = data_list, chains = 4,
                              pars = c("d", "OR", "prob_harm"), 
                              control = list(adapt_delta = 0.9))
print(stan_FE_Bi_logit_pair)
sims_stan_FE_Bi_logit_pair <- extract(stan_FE_Bi_logit_pair)
hist(sims_stan_FE_Bi_logit_pair$d[, 2])

# Inspection:
mcmc_dens(stan_FE_Bi_logit_pair, pars = "d[2]")
mcmc_trace(stan_FE_Bi_logit_pair, pars = "d[2]")

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
writeLines(text = Ch2_RE_Bi_logit_pair, con = "jags/02_JAGS/Ch2_RE_Bi_logit_pair.txt")

jags_Model <- jags(model.file = "jags/02_JAGS/Ch2_RE_Bi_logit_pair.txt", 
                   data = data_list, parameters.to.save = c(
                    "d", "OR", "prob_harm"))
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_RE_Bi_logit_pair <- stan(file = "stan/02_Stan/Ch2_RE_Bi_logit_pair.stan",
                              data = data_list, chains = 4,
                              pars = c("d", "OR", "prob_harm"),
                              control = list(adapt_delta = 0.9)
                              )
print(stan_RE_Bi_logit_pair)
sims_stan_RE_Bi_logit_pair <- extract(stan_FE_Bi_logit_pair)
hist(sims_stan_FE_Bi_logit_pair$d[, 2])

# Inspection:
mcmc_dens(stan_FE_Bi_logit_pair, pars = "d[2]")
mcmc_trace(stan_FE_Bi_logit_pair, pars = "d[2]")

# Success! Both JAGS and Stan are giving similar output...

## Extension to Indirect Comparisons and Network Meta-Analysis -------------
# For NMA, we must assume the exchangeability of both treatment effects over 
# both 1 versus 2 *and* 1 versus 3 trials. The theory extends readily to 
# additional treatments k = 4. 5, ..., S. In each case, we must assume
# exchangeability of the \Delta's across the **entire** set of trials. Then,
# the within-trial transitivity relation is enough to imoly the exchangeability
# of all the treatment effects \Delta_{i, xy}. The consistency equations

# d_{2, 3} = d_{1, 3} - d_{1, 2}
# d_{2, 4} = d_{1, 4} - d_{1, 4}
# .
# .
# .
# d_{s-1, s} = d_{1, S} - d_{1, S - 1}

# are also therefore implied. 

### Fixed Effects -----------------------------------------------------------
##### Data --------------------------------------------------------------------
# N studies:
n_s <- 36
# N treatments:
n_t <- 7
# N trial arms:
n_a <- c(3, rep(2, 35))
# Treatment:
t <- matrix(c(
 rep(1, 19), rep(2, 3), rep(3, 14),
 rep(3, 1), rep(2, 8), 4, 5, rep(7, 11), rep(5, 2), 6, rep(7, 11),
 4, rep(NA, 35)),
 ncol = 3)
# Events:
r <- matrix(c(
 1472, 3, 12, 7, 10, 887, 5, 1455, 9, 4, 285, 11, 1, 8, 1, 4, 14, 9, 42, 2, 13,
 2, 13, 356, 522, 3, 10, 40, 5, 5, 2, 19, 59, 5, 17, 8,
 652, 3, 7, 4, 5, 929, 2, 1418, 6, 6, 270, 2, 3, 5, 1, 0, 7, 3, 29, 3, 5, 2, 7,
 757, 523, 1, 3, 32, 5, 3, 3, 20, 52, 2, 12, 6,
 723, rep(NA, 35)),
 ncol = 3)
# Observations:
n <- matrix(c(
 20251, 65, 159, 85, 135, 10396, 63, 13780, 130, 107, 3004, 149, 50, 58, 53, 45,
 99, 41, 421, 44, 200, 56, 155, 4921, 8488, 55, 94, 573, 75, 69, 61, 419, 782, 
 81, 226, 66,
 10396, 64, 157, 86, 135, 10372, 59, 13746, 123, 109, 3006, 152, 50, 54, 47, 42,
 101, 46, 429, 46, 195, 47, 169, 10138, 8461, 55, 95, 565, 75, 71, 62, 421, 790,
 81, 225, 71,
 10374, rep(NA, 35)), 
 ncol = 3)

# Data list:
data_list <- list(n_s = n_s, n_t = n_t, n_a = n_a, t = t, r = r, n = n)

#### Jags model --------------------------------------------------------------
# model
Ch2_FE_Bi_logit <- "# Binomial likelihood, logit link
# NMA model
# Fixed effects model
model {                    # *** PROGRAMME STARTS
 for (i in 1:n_s) {        # LOOP THROUGH STUDIES
  mu[i] ~ dnorm(0, 0.0001) # vague prior for i'th trial baselines
  for (j in 1:n_a[i]) {
   r[i, j] ~ dbin(p[i, j], n[i, j]) # binomial model
   logit(p[i, j]) <- mu[i] + d[t[i, j]] - d[t[i, 1]] # y_hat
   }
 }
  d[1] <- 0
  for (k in 2:n_t) {
   d[k] ~ dnorm(0, 0.0001) # vague priors for k treatment effects
  }
}                          # *** PROGRAMME ENDS

"
writeLines(text = Ch2_FE_Bi_logit, con = "jags/02_JAGS/Ch2_FE_Bi_logit.txt")

jags_Model <- jags(model.file = "jags/02_JAGS/Ch2_FE_Bi_logit.txt", 
                   data = data_list, parameters.to.save = c(
                    "d")
                   )
print(jags_Model)

#### Stan model --------------------------------------------------------------
# Since Stan can't handle NA values, we have to index the data in a different 
# way. A simple approach is to impute an arbitrary value and then exclude the
# loop from any NA values (so that we do not loop over NAs in Stan model so 
# these do not affect inference). We use -1 since we know that this value is 
# not possible given the model:
r[is.na(r)] <- -1
n[is.na(n)] <- -1
t[is.na(t)] <- -1
# Data list:
data_list <- list(n_s = n_s, n_t = n_t, n_a = n_a, t = t, r = r, n = n, 
                  max_arms = max(n_a))

stan_FE_Bi_logit <- stan(file = "stan/02_Stan/Ch2_FE_Bi_logit.stan",
                         data = data_list, chains = 4,
                         pars = c("d", "OR")
                         )
print(stan_FE_Bi_logit, digits = 4)
# Success! Both JAGS and Stan are giving similar output...

# Posterior inspection:
plot(stan_FE_Bi_logit)
mcmc_dens(stan_FE_Bi_logit, pars = c("d[2]", "d[3]", "d[4]", "d[5]", 
                                          "d[6]", "d[7]"))
mcmc_trace(stan_FE_Bi_logit, pars = c("d[2]", "d[3]", "d[4]", "d[5]", 
                                          "d[6]", "d[7]"))
mcmc_areas(stan_FE_Bi_logit, pars = c("d[2]", "d[3]", "d[4]", "d[5]", 
                                          "d[6]", "d[7]"))
sims_FE_Bi_logit <- rstan::extract(stan_FE_Bi_logit)
# LORs that 1 > 4:
mean(sims_FE_Bi_logit$d[, 1] > sims_FE_Bi_logit$d[, 4])

### Random effects model ----------------------------------------------------
#### Jags model --------------------------------------------------------------
# model
Ch2_RE_Bi_logit <- "# Binomial likelihood, logit link
# NMA model
# Random effects model
model {                    # *** PROGRAMME STARTS
 for (i in 1:n_s) {        # LOOP THROUGH STUDIES
  w[i, 1] <- 0             # 0 adjustment for control arm
  delta[i, 1] <- 0         # 0 treatment effect for control arm
  mu[i] ~ dnorm(0, 0.0001) # vague prior for i'th trial baselines
  for (j in 1:n_a[i]) {    # LOOP THROUGH ARMS
   r[i, j] ~ dbin(p[i, j], n[i, j]) # binomial model
   logit(p[i, j]) <- mu[i] + delta[i, j] # y_hat
  }
  for (j in 2:n_a[i]) {     # LOOP THROUGH ARMS
   delta[i, j] ~ dnorm(m_d[i, j], tau_d[i, j]) # trial-specific LOR dists
   m_d[i, j] <- d[t[i, j]] - d[t[i, 1]] + s_w[i, j] # mu of LOR dists
   tau_d[i, j] <- tau * 2 * (j - 1) / j # precision for LOR dists
   w[i, j] <- (delta[i, j] - d[t[i, j]] + d[t[i, 1]]) # multi-arm adjustment
   s_w[i, j] <- sum(w[i, 1:(j - 1)]) / (j - 1)
  }
 }
  d[1] <- 0
  for (k in 2:n_t) {
   d[k] ~ dnorm(0, 0.0001) # vague priors for k treatment effects
  }
  sd ~ dunif(0, 2)         # vague prior for between-trial stdev
  tau <- pow(sd, -2)       # between-trial precision
}                          # *** PROGRAMME ENDS

"
writeLines(text = Ch2_RE_Bi_logit, con = "jags/Ch2_RE_Bi_logit.txt")

jags_Model <- jags(model.file = "jags/02_JAGS/Ch2_RE_Bi_logit.txt", 
                   data = data_list, parameters.to.save = c(
                    "d")
                   )
print(jags_Model)

#### Stan model --------------------------------------------------------------
r[is.na(r)] <- -1
n[is.na(n)] <- -1
t[is.na(t)] <- -1
# Data list:
data_list <- list(n_s = n_s, n_t = n_t, n_a = n_a, t = t, r = r, n = n, 
                  max_arms = max(n_a))

stan_RE_Bi_logit <- stan(file = "stan/02_Stan/Ch2_RE_Bi_logit.stan",
                         data = data_list, chains = 4,
                         pars = c("d")
                         )
print(stan_RE_Bi_logit, digits = 4)
# Success! Both JAGS and Stan are giving similar output...

# Posterior inspection:
plot(stan_RE_Bi_logit)
mcmc_dens(stan_RE_Bi_logit, pars = c("d[2]", "d[3]", "d[4]", "d[5]", 
                                          "d[6]", "d[7]"))
mcmc_trace(stan_RE_Bi_logit, pars = c("d[2]", "d[3]", "d[4]", "d[5]", 
                                          "d[6]", "d[7]"))
mcmc_areas(stan_RE_Bi_logit, pars = c("d[2]", "d[3]", "d[4]", "d[5]", 
                                          "d[6]", "d[7]"))
sims_RE_Bi_logit <- rstan::extract(stan_RE_Bi_logit)
# LORs that 1 > 4:
mean(sims_RE_Bi_logit$d[, 1] > sims_RE_Bi_logit$d[, 4])

# End file ----------------------------------------------------------------