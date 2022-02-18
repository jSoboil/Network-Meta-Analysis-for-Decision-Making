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
 real mu[n_s];                                 // ave. event baseline
 real <lower = 0, upper = 10> sigmasq_d;       // stdev of d, implicit uniform
 real d_sample;   
}
transformed parameters {
 real <lower = 0> sigma_d;            // parameter for d stdev
 real d[2];                      // ave. treatment effect
 d[1] = 0;
 d[2] = d_sample;
                                      // treatment
 sigma_d = sqrt(sigmasq_d);           // var to stdev
}
model {
 for (i in 1:n_s) {
  // Prior on trial baseline for Pr of event
  mu[i] ~ normal(0, sqrt(1.0E4));
  for (k in 1:2) {
   // Likelihood
   r[i, k] ~ binomial_logit(n[i, k], mu[i] + d[k]);
   }
  }
  
  d_sample ~ normal(0, sigma_d);
}
// END FILE
