// Random Effects
// Pair-wise meta-analysis
// Binomial model, logit link
data {
 int <lower = 0> s[73];           // n studies
 int <lower = 0> r[73];           // n events
 int <lower = 0> n[73];           // n obs
 int <lower = 0> t[73];           // trt indicator
 int <lower = 0> t_b[73];         // trial baseline arm
}
parameters {
 real mu[36];                     // study baseline parameter (36 studies)
 real d_draws[7];                 // parameter to transform to d (7 treatments)
}
transformed parameters {
 real OR[7];                      // OR parameter (7 treatments)
 real prob_harm[7];               // transform OR to p parameter
 real d[7];                       // treatment effect (7 treatments)
 d[1] = 0;                        // 0 effect for reference treatment (1)
 OR[1] = 0;                       // 0 effect for reference treatment (1)
 prob_harm[1] = 0;                // 0 effect for reference treatment (1)
 for (k in 2:7) {                 // LOOP THROUGH k TREATMENTS
  d[k] = d_draws[k];              // trial-specific LORs (studies 2-7)
  OR[k] = exp(d[k]);              // OR of treatment effect (studies 2-7)
  prob_harm[k] = step(d[k]);      // probability of trt effect (studies 2-7)
 }                                // END LOOP
}
model {
 // Priors
 for (k in 2:7) {                 // LOOP THROUGH k TREATMENTS
  d_draws[k] ~ normal(0, 100);    // treatment effect (studies 2-7)
 }                                // END LOOP
 for (j in 1:36) {                // LOOP THROUGH j STUDIES
  mu[j] ~ normal(0, 100);         // baseline effect (for j study)
 }                                // END LOOP
 for (i in 1:73) {                // LOOP THROUGH i DATA FOR EACH ARM
  // Likelihood
   r[i] ~ binomial_logit(n[i], mu[s[i]] + d[t[i]] - d[t_b[i]]);
   }                              // END LOOP
}
// END FILE
