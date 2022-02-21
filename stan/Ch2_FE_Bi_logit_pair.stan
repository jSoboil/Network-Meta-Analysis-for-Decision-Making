// Fixed Effects
// Pair-wise meta-analysis
// Binomial model, logit link
data {
 int <lower = 0> n_s;                   // n studies
 int <lower = 0> n_arms;                // n trial arms
 int <lower = 0> r[n_s, n_arms];        // n events
 int <lower = 0> n[n_s, n_arms];        // n obs
}
parameters {
 real mu[n_s];                                 // ave. i'th trial baseline
 real <lower = 0, upper = 10> sigmasq_d;       // stdev of ATE, implicit uniform
 real d_sample;   
}
transformed parameters {
 real d[n_arms];                   // ave. treatment effect
 real OR[n_arms];
 real <lower = 0, upper = 1> prob_harm;
 d[1] = 0;                         // 0 effect for reference treatment
 d[2] = d_sample;                  // add sample to [i + 1] ATE index
 OR = exp(d);
 prob_harm = step(d[2]);
 
}
model {
 // Priors
 d_sample ~ normal(0, sqrt(1.0E4));             // on ATE
 for (i in 1:n_s) {
  mu[i] ~ normal(0, sqrt(1.0E4));               // on i'th trial baseline
  for (k in 1:2) {
   // Likelihood
   r[i, k] ~ binomial_logit(n[i, k], mu[i] + d[k]);
   }
  }
}
// END FILE
