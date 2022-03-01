# Libraries ---------------------------------------------------------------
pkgs <- c("R2jags", "rstan", "rstanarm", "bayesplot", "dplyr", 
          "tidyr", "igraph")
sapply(pkgs, require, character.only = TRUE)
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
writeLines(text = Ch2_FE_Bi_logit_pair, con = "jags/Ch2_FE_Bi_logit_pair.txt")

jags_Model <- jags(model.file = "jags/Ch2_FE_Bi_logit_pair.txt", 
                   data = data_list, parameters.to.save = c(
                    "d", "OR", "prob_harm"))
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_FE_Bi_logit_pair <- stan(file = "stan/Ch2_FE_Bi_logit_pair.stan",
                              data = data_list, chains = 4,
                              pars = c("d", "OR", "prob_harm"), 
                              control = list(adapt_delta = 0.9)
                              )
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
writeLines(text = Ch2_RE_Bi_logit_pair, con = "jags/Ch2_RE_Bi_logit_pair.txt")

jags_Model <- jags(model.file = "jags/Ch2_RE_Bi_logit_pair.txt", 
                   data = data_list, parameters.to.save = c(
                    "d", "OR", "prob_harm"))
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_RE_Bi_logit_pair <- stan(file = "stan/Ch2_RE_Bi_logit_pair.stan",
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
# N studies:
n_s <- 11
# N treatments:
n_t <- 7
# N trial arms:
n_arms <- c(3, rep(2, 35))
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
data_list <- list(n_s = n_s, n_t = n_t, n_arms = n_arms, t = t, r = r, n = n)

# Since Stan can't handle NA values, we have to index the data in a different 
# way. This can be done in the same way to the more efficient approach for 
# parsing data to BUGS, using nested indexing. Basically, we group rows by the 
# number of arms in each study, so the data for a two arm study will comprise 
# two rows, and so on. This avoids the use of NA values across columns, 
# especially when only a few studies have more than one comparison.

# N studies:
n_s <- 11
# N treatments:
n_t <- 7
# Study number:
s_n <- c(1:36)
# N trial arms:
n_arms <- c(3, rep(2, 35))
# Treatment:
t_1 <- data_list$t[, 1]
t_2 <- data_list$t[, 2]
t_3 <- data_list$t[, 3]
# Events:
r_1 <- data_list$r[, 1]
r_2 <- data_list$r[, 2]
r_3 <- data_list$r[, 3]
# Observations:
n_1 <- data_list$n[, 1]
n_2 <- data_list$n[, 2]
n_3 <- data_list$n[, 3]
# Study number
thrombo_data <- tibble(r_1 = r_1, n_1 = n_1, r_2 = r_2, n_2 = n_2, r_3 = r_3,
                       n_3 = n_3, t_1 = t_1, t_2 = t_2, t_3 = t_3, 
                       n_arms = n_arms, s_n = s_n, t_b = t_1)


# Make into long (tidy) format
df_thrombo <- thrombo_data |>
 group_by(s_n) |>
 gather(key, var, r_1:t_3) |>
 mutate(n_arm = gsub("[rnt]+_", "\\1", key),
        key = gsub("\\_*[1-7]", "\\1", key)) |>
 spread(key, var) |>
 filter(!is.na(t)) |>
 ungroup() |>
 transmute(s = s_n,
           t = t,
           r = r,
           n = n, 
           t_b = t_b) |>
 arrange(s, t)
df_thrombo
df_thrombo <- as.list(df_thrombo)

#### Jags model --------------------------------------------------------------
# model
Ch2_FE_Bi_logit <- "# Binomial likelihood, logit link
# Pairwise meta-analysis
# Random effects model
model {                    # *** PROGRAMME STARTS
 for (j in 1:36) {         # LOOP THROUGH STUDIES
  mu[j] ~ dnorm(0, 0.0001) # vague prior for trial baselines
  }
 
 for (k in 2:7) {          # LOOP THROUGH TREATMENTS
  d[k] ~ dnorm(0, 0.0001)  # prior on treatment effects
  OR[k] <- exp(d[k])           # LOR to OR
  prob_harm[k] <- step(d[k])   # OR to p
 }
 
  for (i in 1:73) {        # LOOP THROUGH DATA
   r[i] ~ dbinom(p[i], n[i]) # binomial likelihood
   logit(p[i]) <- mu[s[i]] + d[t[i]] - d[t_b[i]] # model for linear predictor
  }
  d[1] <- 0                # treatment effect is zero for reference treatment
}                          # *** PROGRAMME ENDS

"
writeLines(text = Ch2_FE_Bi_logit, con = "jags/Ch2_FE_Bi_logit.txt")

jags_Model <- jags(model.file = "jags/Ch2_FE_Bi_logit.txt", 
                   data = df_thrombo, parameters.to.save = c(
                    "d", "OR", "prob_harm")
                   )
print(jags_Model)

#### Stan model --------------------------------------------------------------
stan_FE_Bi_logit <- stan(file = "stan/Ch2_FE_Bi_logit.stan",
                         data = df_thrombo, chains = 4,
                         pars = c("d", "OR", "prob_harm"),
                         control = list(adapt_delta = 0.9)
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