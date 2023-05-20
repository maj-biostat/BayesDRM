# BayesDRM

Bayesian implementation of common dose response models using stan.

## Installation

```
git clone https://github.com/maj-biostat/BayesDRM.git
cd BayesDRM
R CMD INSTALL .
```

## 2-parameter emax - upper limit assumed to be unity 

$$
\begin{aligned}
y_i &\sim Binomial(n_i, p_i) \\
p_i &= p_0 + p_{emax} * \frac{x_i}{x_i + b50} \\
p_0 &\sim Beta(a, b) \\
p_{emax} &= 1 - p_0 \\
b50 &\sim Normal^+(\mu_{b50}, \sigma_{b50})
\end{aligned}
$$

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



## 3-parameter emax

$$
\begin{aligned}
y_i &\sim Binomial(n_i, p_i) \\
p_i &= inv_logit(\eta_i) \\
\eta_i &= b_0 + b_2 * \frac{x_i}{x_i + b50} \\
b_0 &\sim Normal(\mu_{b0}, \sigma_{b0}) \\
b_{max} &\sim Normal(\mu_{bmax}, \sigma_{bmax}) \\
b50 &\sim Normal^+(\mu_{b50}, \sigma_{b50}) \\
b_2 &= b_{max} - b_0 \\
\end{aligned}
$$

Simulate data under this model.

```
emax_3p <- function(x, b0, bmax, b50){
  b2 <- bmax - b0
  eta <- b0 + b2 * x / (x + b50)
  plogis(eta)
}

emax_3p_med <- function(pr_med, b0, bmax, b50){
  eta_med <- qlogis(pr_med)
  b2 <- bmax - b0
  b50 * (eta_med - b0) / (b2 - eta_med + b0)
}


b0 <- qlogis(0.15)
bmax <- qlogis(0.95)
b50 <- 3
pr_med <- 0.8
x <- c(0, 5, 10, 15, 20)
p <- emax_3p(x, b0, bmax, b50)
med <- emax_3p_med(pr_med, b0, bmax, b50)
p
K <- length(p)
# 30 obs per dose
n <- 30
y <- rbinom(K, n, p)

plot(x, p, ylim = c(0, 1), type = "l")
points(x, y/n)

ld <- list(N = K, y = y, n = rep(n, K), x = x,
           pri_b0_mu = -1.4, pri_b0_s = 2,
           pri_b50_mu = 4, pri_b50_s = 3, 
           pri_bmax_mu = 3, pri_bmax_s = 0.5,
           pr_med = 0.8, prior_only = F)
f1 <- BayesDRM::drm_emax3_bin(ld, refresh = 0)
f1
```









