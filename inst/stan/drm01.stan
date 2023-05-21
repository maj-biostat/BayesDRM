data {
  int N;
  int y[N];
  int n[N];
  vector[N] x;
  real pri_p0_a;
  real pri_p0_b;
  real pri_b50_mu;
  real pri_b50_s;  
  real<lower=0,upper=1> pr_med;
  int prior_only;
}
parameters {
  real<lower=0,upper=1> p0;
  real<lower=0> b50;  
}
transformed parameters{
  vector<lower=0,upper=1>[N] p;
  real pemax = 1-p0;
  p = p0 + pemax * x .* inv(x + b50);
  
}
model {
  target += beta_lpdf(p0 | pri_p0_a, pri_p0_b);
  target += normal_lpdf(b50 | pri_b50_mu, pri_b50_s);
  if(!prior_only){ target += binomial_lpmf(y | n, p); }  
}
generated quantities{
  int yrep[N];
  real med;
  vector[N] log_lik;
  med = b50 * (pr_med - p0) * inv(1 - p0);
  yrep = binomial_rng(n, p);
  for(i in 1:N) { log_lik[i] = binomial_lpmf(y[i] | n[i], p[i]); }
}


