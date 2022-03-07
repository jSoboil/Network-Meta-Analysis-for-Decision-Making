// Random Effects
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
 real mu[n_s];                   // declare baseline parameter
 real d_samples[n_t - 1];        // declare d_sample parameter
 real delta_samples[n_s, max_arms - 1];    // declare delta_sample parameter
 real <lower = 0> stdev;         // between trial standard deviation
}
transformed parameters {
 real p[n_s, max_arms];         // declare p of success parameter
 real w[n_s, max_arms];         // multi-arm adjustment parameter
 real s_w[n_s, max_arms - 1];   // cumulative adjustment parameter
 real m_d[n_s, max_arms - 1];   // mean of LOR distributions parameter
 real tau_d[n_s, max_arms - 1]; // precision for LOR distributions parameter
 real sd_d[n_s, max_arms - 1];  // stdev for LOR distributions parameter
 real <lower = 0> tau;          // between trial precision
 real d_1;                      // d_{1} variable for reference treatment
 real d[n_t];                   // d_{n_t} variable for comparison treatment
 real delta_1[n_s];
 real delta[n_s, max_arms];     // declare delta variable
 tau = pow(stdev, -2);
 d_1 = 0;                       // treatment effect 0 for reference treatment
 d[1] = d_1;                    
 d[2:n_t] = d_samples;          // index [i + 1] treatment effects
 delta_1 = rep_array(0, n_s);
 delta[, 1] = delta_1;
 delta[, 2:max_arms] = delta_samples;
 for (i in 1:n_s) {             // LOOP THROUGH STUDIES
  w[i, 1] = 0;                  // adjustment for multi-arms = 0 for control
   for (j in 1:n_a[i]) {        // LOOP THROUGH ARMS
     delta[i, j] = d[t[i, j]] - d[t[i, 1]];    // effect variance set to 0
     p[i, j] = inv_logit(mu[i] + delta[i, j]); // y_hat
   }
   for (j in 2:n_a[i]) {
    w[i, j] = (delta[i, j] - d[t[i, j]] + d[t[i, 1]]); // multi-arm adjustment
    s_w[i, j - 1] = sum(w[i, 1:(j - 1)]) / (j - 1); // cumulative adjustment
    m_d[i, j - 1] = d[t[i, j]] - d[t[i, 1]] + s_w[i, j - 1]; // mean of LORs with multi-arm correction
    tau_d[i, j - 1] = tau * 2 * (j - 1) / j; // precision of LORs with multi-arm trial correction
    sd_d[i, j - 1] = pow(tau_d[i, j - 1], -0.5);
   }
 }
}
model {
 // Priors
 stdev ~ uniform(0, 2);
 d_samples ~ normal(0, 100000); // prior for k'th treatment effect
 for (i in 1:n_s) {             // LOOP THROUGH STUDIES
  mu[i] ~ normal(0, 100000);    // prior for i'th trial baseline
   for (j in 1:n_a[i]) {        // LOOP THROUGH ARMS
    // Likelihood
     r[i, j] ~ binomial(n[i, j], p[i, j]);
   }
   for (j in 1:(n_a[i] - 1)) {
     delta_samples[i, j] ~ normal(m_d[i, j], sd_d[i, j]);
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
