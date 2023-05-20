# BayesDRM

Bayesian implementation of common dose response models using stan.

```
git clone https://github.com/maj-biostat/BayesDRM.git
cd BayesDRM
R CMD INSTALL .
```

Simulate data under this model.

```
emax_2p <- function(x, p0, b50){
  pemax <- 1 - p0
  p0 + pemax * x / (x + b50)
}

emax_2p_med <- function(pr_med, p0, b50){
  b50 * (pr_med - p0) / (1-p0)
}

p0 <- 0.15
b50 <- 2
pr_med <- 0.8
x <- c(0, 5, 10, 15, 20)
p <- emax_2p(x, p0, b50)
med <- emax_2p_med(pr_med, p0, b50)
p
K <- length(p)
# 30 obs per dose
n <- 30
y <- rbinom(K, n, p)

plot(x, p, ylim = c(0, 1), type = "l")
points(x, y/n)


y <- rbinom(K, n, p)
ld <- list(N = K, y = y, n = rep(n, K), x = x,
           pri_p0_a = 1, pri_p0_b = 1,
           pri_b50_mu = 10, pri_b50_s = 4,
           prior_only = F, pr_med = 0.8)
f1 <- BayesDRM::drm_emax2_bin(ld, refresh = 0)
f1
m <- as.matrix(f1, pars = c("p0", "b50", "p", "med", "yrep"))
```









