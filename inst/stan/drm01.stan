data {
  int N;
  int y[N];
  int n[N];
  vector[N] x;
  real pri_p0_a;
  real pri_p0_b;
  real pri_b50_mu;
  real pri_b50_s;  
  int prior_only;
  real<lower=0,upper=1> pr_med;
}
parameters {
  real<lower=0,upper=1> p0;
  real<lower=0> b50;  
}
transformed parameters{
  vector<lower=0,upper=1>[N] p;
  real pemax = 1-p0;
  for(i in 1:N){
    p[i] = p0 + pemax * x[i] * inv(x[i] + b50);
  }
}
model {
  target += beta_lpdf(p0 | pri_p0_a, pri_p0_b);
  target += normal_lpdf(b50 | pri_b50_mu, pri_b50_s);
  if(!prior_only){
    target += binomial_lpmf(y | n, p);
  }  
}
generated quantities{
  real med;
  med = b50 * (pr_med - p0) * inv(1 - p0);
}


