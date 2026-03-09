Scode0 <- '
data {
  int<lower=1> n;                 // observations
  int<lower=1> p;
  vector[p] y[n];
  real<lower=0> eta;
 }
parameters {
  cholesky_factor_corr[p] L_Omega;
  vector[p] x[n];
}
transformed parameters {
  corr_matrix[p] rho = tcrossprod(L_Omega);
}
model {
  L_Omega ~ lkj_corr_cholesky(eta);
  x ~ multi_normal_cholesky(rep_vector(0.0, p), L_Omega);
  for(j in 1:p) {
    y[j] ~ gamma(exp(x[j]), 1.0); 
  }
}
'

## compile STAN code
library(rstan)
options(mc.cores = 4L)

system.time(
    S0_cmpld <- stan_model(
        model_code = Scode0,
        model_name = "g2"
    )
)

## STEP 2: prepare the data
data(SAheart, package = "msos")

str(SAheart)

(n <- nrow(SAheart))

jj <- pmatch(c("tobacco", "alcohol"), colnames(SAheart))
sdat <- sweep(SAheart[, jj], 2, sapply(SAheart[, jj], sd), "/")

cov(sdat)
cor(log(sdat+0.001))

## STAN (initial) data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(2),
    y = as.matrix(sdat),
    eta = 2
)

Samples <- sampling(
    S0_cmpld, 
    data = Sdata0,
    iter = 3000,
    warmup = 500,
    thin = 10,
    chains = 4
)

## PLOTS
library(coda)

p=2

rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

library(ggpubr)

## CI for the graphpcor model parameters (all edges are 'significant')

ycorr <- cor(sdat)
ylcorr <- cor(sapply(sdat, function(x)
    ifelse(x<0.1,-3,log(x))), use = "pair")
ycorr
ylcorr

iil <- which(lower.tri(ycorr))
iil

## CI for the correlation pairs (observed in 'red')
## g0 not enough
mcmc_intervals(Samples, rhonams) +
    geom_segment(aes(x = x, y = ya,
                     xend = x, yend = yb, color = ctype),
                 data.frame(
                     x = c(ycorr[iil], ylcorr[iil]),
                     ya = length(iil):1-0.1,
                     yb = length(iil):1+0.5,
                     ctype = rep(c('y', 'log(y)'), each = length(iil))))
