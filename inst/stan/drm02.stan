data {
  int N;
  int y[N];
  int n[N];
  vector[N] x;
  real pri_b0_mu;
  real pri_b0_s;
  real pri_b50_mu;
  real pri_b50_s;  
  real pri_bmax_mu;
  real pri_bmax_s; 
  real<lower=0,upper=1> pr_med;
  int prior_only;
}
transformed data{
  real eta_med = logit(pr_med);
}
parameters {
  real b0;
  real<lower=0> b50; 
  real bmax; 
}
transformed parameters{
  real b2;
  vector[N] eta;
  b2 = bmax - b0;
  eta = b0 + b2 * x .* inv(x + b50);
}
model {
  target += normal_lpdf(b0 | pri_b0_mu, pri_b0_s);
  target += normal_lpdf(b50 | pri_b50_mu, pri_b50_s);
  target += normal_lpdf(bmax | pri_bmax_mu, pri_bmax_s);
  if(!prior_only){ target += binomial_logit_lpmf(y | n, eta); }  
}
generated quantities{
  int yrep[N];
  real med;
  vector[N] p;
  p = inv_logit(eta);
  med = b50 * (eta_med - b0) .* inv(b2 - eta_med + b0);
  yrep = binomial_rng(n, inv_logit(eta));
}


