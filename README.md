# BayesDRM

Bayesian implementation of common dose response models using stan.

## Installation

```
git clone https://github.com/maj-biostat/BayesDRM.git
cd BayesDRM
R CMD INSTALL .
```

## Examples

You will need stan installed along with `data.table`, `ggplot2`, `kableExtra` and `loo` to run some of these.


### 2-parameter emax - upper limit assumed to be unity 

Parameterised in terms of upper and lower asymptote.
Emax parameter is derived.

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

# Leave one out cross validation:

log_lik_1 <- loo::extract_log_lik(f1, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
```


### 3-parameter emax

Parameterised in terms of upper and lower asymptote.
Emax parameter is derived.

$$
\begin{aligned}
y_i &\sim Binomial(n_i, p_i) \\
p_i &= expit(\eta_i) \\
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

Smaller data set:


```
ld <- list(N = 5, 
           y = c(4, 5, 5, 10, 10),
           n = c(10, 10, 10, 10, 10), 
           x = c(0, 3, 6, 12, 20),
           pri_b0_mu = 0, pri_b0_s = 1.5,
           pri_b50_mu = 4, pri_b50_s = 3, 
           # note - a stronger prior on the upper will also force 
           # down the lower bound.
           pri_bmax_mu = 7, pri_bmax_s = 0.8,
           pr_med = 0.8, prior_only = F)
f1 <- BayesDRM::drm_emax3_bin(ld, refresh = 0)
f1

dfig <- data.table(as.matrix(f1, pars = "p"))
dfig <- melt(dfig, measure.vars = names(dfig))
dfig <- dfig[, .(
  mu = mean(value),
  q_025 = quantile(value, prob = 0.025),
  q_975 = quantile(value, prob = 0.975)
), keyby = variable]
dfig[, x := gsub("p[", "", variable, fixed = T)]
dfig[, x := as.integer(gsub("]", "", x, fixed = T))]
dfig[, x := ld$x[x]]

kableExtra::kable(dfig, format = "simple", digits = 2)

variable      mu   q_025   q_975    x
---------  -----  ------  ------  ---
p[1]        0.22    0.08    0.42    0
p[2]        0.61    0.42    0.78    3
p[3]        0.82    0.69    0.92    6
p[4]        0.94    0.89    0.98   12
p[5]        0.98    0.95    0.99   20
```


### 4-parameter emax

Parameterised in terms of upper and lower asymptote.
Emax parameter is derived.
Hill parameter is constrained by config, otherwise will usually be bimodal.

$$
\begin{aligned}
y_i &\sim Binomial(n_i, p_i) \\
p_i &= expit(\eta_i) \\
\eta_i &= b_0 + b_2 * \frac{x_i^{b_h}}{x_i^{b_h} + b50^{b_h}} \\
b_0 &\sim Normal(\mu_{b0}, \sigma_{b0}) \\
b_{max} &\sim Normal(\mu_{bmax}, \sigma_{bmax}) \\
b50 &\sim Normal^+(\mu_{b50}, \sigma_{b50}) \\
b_h &\sim Normal(\mu_{b_h}, \sigma_{b_h}) \\
b_2 &= b_{max} - b_0 \\
\end{aligned}
$$

Simulate data under this model.

```
emax_4p <- function(x, b0, bmax, bh, b50){
  b2 <- bmax - b0
  eta <- b0 + b2 * x^bh / (x^bh + b50^bh)
  plogis(eta)
}

b0 <- qlogis(0.15)
bmax <- qlogis(0.95)
b50 <- 7
bh <- 3
x <- c(0, 5, 10, 15, 20)
p <- emax_4p(x, b0, bmax, bh, b50)
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
           pri_bh_mu = 1, pri_bh_s = 2,
           # constrain hill to positive
           pos_hill = 1, prior_only = F)
f1 <- BayesDRM::drm_emax4_bin(ld, refresh = 0)
f1
```





