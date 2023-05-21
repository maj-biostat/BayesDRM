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
  real pri_bh_mu;
  real pri_bh_s; 
  // need to constrain hill otherwise is usually bimodal
  int<lower = 0, upper = 1> pos_hill;
  int prior_only;
}
transformed data{
}
parameters {
  real b0;
  real<lower=0> b50; 
  real bmax; 
  vector<lower=0>[pos_hill ? 1 : 0] bh1;
  vector<upper=0>[!pos_hill ? 1 : 0] bh2;
}
transformed parameters{
  real b2;
  vector[N] eta;
  real bh;
  if(pos_hill){ bh = bh1[1]; } else { bh = bh2[1]; };
  b2 = bmax - b0;
  for(i in 1:N){
    real bt = pow(x[i], bh);
    eta[i] = b0 + b2 * bt * inv(bt + pow(b50, bh));
  }
}
model {
  target += normal_lpdf(b0 | pri_b0_mu, pri_b0_s);
  if(pos_hill) { target += normal_lpdf(bh1 | pri_bh_mu, pri_bh_s); }
  if(!pos_hill) { target += normal_lpdf(bh2 | pri_bh_mu, pri_bh_s); }
  target += normal_lpdf(b50 | pri_b50_mu, pri_b50_s);
  target += normal_lpdf(bmax | pri_bmax_mu, pri_bmax_s);
  if(!prior_only){ target += binomial_logit_lpmf(y | n, eta); }  
}
generated quantities{
  int yrep[N];
  vector[N] p;
  vector[N] log_lik;
  p = inv_logit(eta);
  yrep = binomial_rng(n, inv_logit(eta));
  for(i in 1:N) { log_lik[i] = binomial_lpmf(y[i] | n[i], p[i]); }
}


