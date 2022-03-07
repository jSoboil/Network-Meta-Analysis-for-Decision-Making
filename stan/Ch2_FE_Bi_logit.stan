// Fixed Effects
// NMA
// Binomial model, logit link
data {
 int n_s;                        // n studies
 int n_t;                        // n treatments
 int n_a[n_s];                   // n arms
 int max_arms;                   // indexing variable for max arms
 int r[n_s, max_arms];           // r events [studies, max arms]
 int n[n_s, max_arms];           // r observations [studies, max arms]
 int t[n_s, max_arms];           // t treatments [studies, max arms]
}
parameters {
 real mu[n_s];                   // declare baseline variable
 real d_samples[n_t - 1];        // declare d_sample variable
}
transformed parameters {
 real delta[n_s, max_arms];     // declare delta variable
 real p[n_s, max_arms];         // declare p variable
 real d_1;                      // d_{1} variable for reference treatment
 real d[n_t];                   // d_{n_t} variable for comparison treatment
 d_1 = 0;                       // treatment effect 0 for reference treatment
 d[1] = d_1;                    
 d[2:n_t] = d_samples;          // index [i + 1] treatment effects
 for (i in 1:n_s) {             // LOOP THROUGH STUDIES
   for (j in 1:n_a[i]) {        // LOOP THROUGH ARMS
     delta[i, j] = d[t[i, j]] - d[t[i, 1]];    // effect variance set to 0
     p[i, j] = inv_logit(mu[i] + delta[i, j]); // y_hat
   }
 }
}
model {
 // Priors
 d_samples ~ normal(0, 100000); // vague prior for k'th treatment effect
 for (i in 1:n_s) {              // LOOP THROUGH STUDIES
  mu[i] ~ normal(0, 100000);        // vague prior for i'th trial baseline
   for (j in 1:n_a[i]) {         // LOOP THROUGH ARMS
    // Likelihood
     r[i, j] ~ binomial(n[i, j], p[i, j]);
   }
 }
}
generated quantities {
 real OR[n_t, n_t];
 real p_t[n_t];
 for (i in 1:(n_t - 1)) {
   OR[i, i] = 1;
   for (j in (i + 1):n_t) {
     OR[i, j] = exp(d[i] - d[j]);
     OR[j, i] = 1 / OR[i, j];
     }
   }
  OR[n_t, n_t] = 1;
 }
// END FILE
