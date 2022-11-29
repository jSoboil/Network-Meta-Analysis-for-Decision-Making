// Random Effects
// Pair-wise meta-analysis
// Binomial model, logit link
data {
 int <lower = 0> n_s;              // n studies
 int <lower = 0> n_arms;           // n trial arms
 int <lower = 0> r[n_s, n_arms];   // n events
 int <lower = 0> n[n_s, n_arms];   // n obs
}
parameters {
 real mu[n_s];                       // i'th trial baseline parameter
 real d_draws;                       // parameter to transform to d
 real delta_draws[n_s, n_arms];      // parameter to transform to delta
 real <lower = 0, upper = 2> tau;    // between-trial stdev parameter
}
transformed parameters {
 real OR;                           // OR parameter
 real prob_harm;                    // transform OR to p parameter
 real d[n_arms];                    // ave. treatment effect (ATE)
 real delta[n_s, n_arms];           // trial-specific LORs parameter
 d[1] = 0;                          // 0 effect for reference treatment
 d[2] = d_draws;                    // add sample to [i + 1] ATE index
 for (i in 1:n_s) {                 // LOOP THROUGH STUDIES
  delta[i, 1] = 0;                  // i'th ATE == 0 for control arm
  delta[i, 2] = delta_draws[i, 2];  // i'th trial-specific LOR distributions 
 }                                  // END LOOP
 OR = exp(d[2]);                    // OR of treatment effect
 prob_harm = step(d[2]);            // probability of treatment effect
}
model {
 for (i in 1:n_s) {                    // LOOP THROUGH STUDIES
  mu[i] ~ normal(0, 100);              // i'th trial baseline effect
  delta_draws[i] ~ normal(d[2], tau);  // samples for trial-specific LORs
  for (k in 1:2) {                     // LOOP THROUGH ARMS
  // Likelihood
   r[i, k] ~ binomial_logit(n[i, k], mu[i] + d[k] + delta[i, k]);
   }                                   // END LOOP
 }                                     // END LOOP
  // Priors
  d_draws ~ normal(0, 100);            // ATE prior
}
// END FILE
