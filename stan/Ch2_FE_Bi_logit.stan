data {
 int <lower = 0> n_s;         // number of studies
 int <lower = 0> n_t;         // number of treatments
 int <lower = 0> n_arms;      // number of arms per study
 real t[n_s, n_arms];                 // i x k treatments
 real r[n_s, n_arms];                 // i x k events
 real n[n_s, n_arms];                 // i x k obs
}
parameters {
 vector[n_s] mu;
 real <lower = 0> sigmasq_delta;  // std of delta
}
transformed parameters {
 // Transformed parameters
 real <lower = 0> sigma_delta; // transformed variance to std
 real d[n_t];
 real OR[n_arms];
 sigma_delta = sqrt(sigmasq_delta);
 d[1] = 0;
}
model {
 for (i in 1:n_s) {
  mu[i] ~ normal(0, sqrt(1.0E4));
  for (k in 1:n_arms) {
   r[n_s, n_arms] ~ binomial_logit(n[n_s, n_arms], );
  }
 }
}
// End file
