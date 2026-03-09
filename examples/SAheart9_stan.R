library(graphpcor)

## STEP 1: data model definition

##   A data frame with 462 observations on the following 10 variables.
##  sbp systolic blood pressure
##  tobacco cumulative tobacco (kg)
##  ldl low density lipoprotein cholesterol
##  adiposity a numeric vector
##  famhist family history of heart disease,
##       a factor with levels ‘"Absent"’ and ‘"Present"’
##  typea type-A behavior
##  obesity a numeric vector
##  alcohol current alcohol consumption
##  age age at onset
##  chd response, coronary heart disease

## The chd variable is the outcome
## There are 8 continuous variables, but two are "non-Gaussian"
## and one binary (famhist) as covariates

## Start by standardization of the continuous 8 covariates by
##   performing the two following steps on the first 6 and
##   only the second step for the covariates 7 and 8
##    step 1: subtract its mean
##    step 2: divide by its standard deviation
## Model the first p-3 contiuous variables as 
##     y[i,j] = x[i,j], j = 1, ..., p-3
## Model the next two (p-2 to p-1) contiuous variables as a Gamma(a, b=1)
##     y[i,j] ~ G1(a=exp(x[i,j]), j = p-2, p-1
## Model famhist, coded as 0 (if 'Absent') or 1 (if 'Present') as
##     y[i,p] ~ bernoulli(inv_logit(alpha_2 + alpha_3 * x[i,p]))
## NOTE: E(y[i,j]) = 0 and V(y[i,j]) = 1, for j in 1, ..., p-3 
##       log(E(y[i,j])=V(y[i,j])) = x[i,j] for j in p-2 and p-1
##       E(y[i,p]) = p_{ij} = Phi(alpha+x[i,p]), Phi(.) is the normal cdf
## Model alpha ~ N(0, 1)
## Model x[i, ] ~ N(0, R), R is a correlation matrix

## STAN model base code
## supply yx as the first p-1 covariates
Scode0 <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  array[n] vector[p] yx;
}
parameters {
  vector[p] mu;
  vector<lower=0>[p] y_sigmas;
  vector<lower=0>[p] x_sigmas;
}
transformed parameters {
  matrix[p,p] rho;
}
model {
  for(j in 1:p) {
    mu[j] ~ normal(0, 10);
    y_sigmas[j] ~ exponential(1);
    x_sigmas[j] ~ exponential(1);
  }
  array[n] vector[p] x;
  x ~ multi_normal(mu, quad_form_diag(rho, x_sigmas));
  for(j in 1:p) {
    for(i in 1:n) {
      yx[i,j] ~ normal(x[i,j], y_sigmas[j]);
    }
  }
}
"

## STEP 2: update the STAN code with code for the
## graphpcor prior for 'rho'
Sgrpc <- stan_add(
    Scode0, 'graphpcor',
    lambda = 1, name = "rho")

##Sgrpc <- Scode0

Sgrpc <- 'data {
  int<lower=1> n;
  int<lower=1> p;
  array[n] vector[p] yx;
    int grpc_dim;
    int grpc_theta_dim;
    real<lower=0> grpc_lambda;
    real<lower=0> grpc_logDetIhalf;
    vector[grpc_theta_dim] grpc_theta_base;
    matrix[grpc_theta_dim, grpc_theta_dim] grpc_Ibase_half;
    vector[grpc_dim] grpc_d0;
    array[grpc_theta_dim] int grpc_ii;
    array[grpc_theta_dim] int grpc_jj;
    int grpc_nfi;
    array[grpc_nfi] int grpc_ifi;
    array[grpc_nfi] int grpc_jfi;
    }
parameters {
  vector[p] mu;
  vector<lower=0>[p] y_sigmas;
  vector<lower=0>[p] x_sigmas;
    vector[grpc_theta_dim] grpc_theta;
    }
transformed parameters {
    vector[grpc_theta_dim] grpc_xi = grpc_Ibase_half * (grpc_theta - grpc_theta_base);
    real<lower=0> grpc_radius = dot_self(grpc_xi);
}
model {
    matrix[grpc_dim,grpc_dim] grpc_L0;
    matrix[grpc_dim,grpc_dim] grpc_V0 = chol2inv(grpc_L0);
    vector[grpc_dim] grpc_s0inv = inv_sqrt(diagonal(grpc_V0));
    matrix[p,p] rho = quad_form_diag(grpc_V0, grpc_s0inv);
    for(i in 1:grpc_dim) {
      for(j in 1:grpc_dim) {
        if(i==j) {
          grpc_L0[i,i] = grpc_d0[i];
        } else {
          grpc_L0[i,j] = 0.0;
        }
      }
    }
    {
      int i,j;
      for(k in 1:grpc_theta_dim) {
        i = grpc_ii[k];
        j = grpc_jj[k];
        grpc_L0[i,j] = grpc_theta[k];
      }
      if(grpc_nfi>0) {
        for(l in 1:grpc_nfi) {
          i = grpc_ifi[l];
          j = grpc_jfi[l];
          if(j>0) {
            for(k in 1:(j-1)) {
              grpc_L0[i,j] -= grpc_L0[i,k] * grpc_L0[j,k] / grpc_L0[j,j];
            }
          }
        }
      }
    }
  for(j in 1:p) {
    mu[j] ~ normal(0, 10);
    y_sigmas[j] ~ exponential(1);
    x_sigmas[j] ~ exponential(1);
  }
  array[n] vector[p] x;
  x ~ multi_normal(mu, quad_form_diag(rho, x_sigmas));
  for(j in 1:p) {
    for(i in 1:n) {
      yx[i,j] ~ normal(x[i,j], y_sigmas[j]);
    }
  }
    grpc_radius ~ exponential(grpc_lambda);
    target += lgamma(grpc_theta_dim * 0.5) -grpc_theta_dim*0.5 * log(pi());
    target += grpc_logDetIhalf - (grpc_theta_dim-1.0)*log(grpc_radius);
  }
'

cat(Sgrpc)

## STEP 3:  compile STAN code
library(rstan)
options(mc.cores = 4L)

system.time(
    Sgpc_cmpld <- stan_model(
        model_code = Sgrpc, 
        model_name = "graphpcor"
    )
)

## STEP 4: prepare the data
data(SAheart, package = "msos")

str(SAheart)

(n <- nrow(SAheart))
(p <- ncol(SAheart)-1)

jj0 <- c(1, 3:4, 6:7, p)
jj1 <- c(2, 8)
jj2 <- c(jj0, jj1)

xdata <- data.frame(
    scale(SAheart[, jj0]), ### first p-2 continuous (standardized)
    (SAheart[, jj1]<0.01)*0.01+SAheart[, jj1], ## p-2 and p-1
    famHist = (SAheart$famhist=="Present")+0L  ## p
); xnames <- colnames(xdata)
for(j in -2:-1+p) ## standardize p-2 and p-1 so its variance = 1
    xdata[, j] <- xdata[, j]/sd(xdata[, j]) 

str(xdata)
sapply(xdata, sd)

## STAN (initial) data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(p),
    chd = as.integer(SAheart[, p + 1]),
    yx = as.matrix(SAheart[, jj2]), 
    famHist = (SAheart$famhist=="Present")+0L, 
    alphas_sigma = 10,
    betas_sigma = 10
)

str(Sdata0)

round((cc <- cor(xdata) * 100))

lcc <- chol(cc)
qc <- chol2inv(lcc)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(xnames, xnames)
round(pC*100, 2)

## define a minimum spanning tree considering the functions
## that I have implemented (>20y ago!) in the spdep package
## this MST consider the weights w_{ij} = 1-|\rho_{ij}|
nb <- lapply(1:p, function(i)
    setdiff(1:p,i)); class(nb) <- 'nb'
nbc <- lapply(1:p, function(i) 1-abs(pC[i, -i]))

library(spdep)
nbw <- nb2listw(nb, nbc, style="B")
mst <- mstree(nbw)

G0 <- matrix(0, p, p, dimnames = dimnames(pC))
for(i in 1:nrow(mst)) {
    G0[mst[i,1], mst[i,2]] <- 1
    G0[mst[i,2], mst[i,1]] <- 1
}
G0

library(graphpcor)
g0 <- graphpcor(G0)
(dg0 <- dim(g0))

p*(p-1)/2

## define a graph consideing partial correlations (abs) > 0.09
tanh(qnorm(c(0.025, 0.975)) / sqrt(n-3))

g1 <- graphpcor(abs(pC[1:6, 1:6])>0.07)
G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

if(FALSE)
    png("g01SAheart.png", width = 1800, height = 900, res = 300)
par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(g0, Rgraphviz = TRUE)
b <- plot(g1, Rgraphviz = TRUE)
if(FALSE)
    dev.off()

if(FALSE)
    system("eog g01SAheart.png &")

## base model
baseM0 <- basepcor(diag(p), iLtheta = g0)
baseM1 <- basepcor(diag(6), iLtheta = g1)

##   update STAN data
SdataM0 <- stan_add(Sdata0, baseM0, lambda = 1, name = 'rho')
SdataM1 <- stan_add(Sdata0, baseM1, lambda = 1, name = 'rho')

SdataM1 <- stan_add(
    list(n = n, p = 6L, yx = xdata[, 1:6]),
    baseM1, lambda = 1, name = 'rho'
)

## STAN sampling
## STAN sampling for each graph
Samples0 <- sampling(
    Sgpc_cmpld, 
    data = SdataM0,
    iter = 30000,
    warmup = 5000,
    thin = 10,
    chains = 4
)

Samples1 <- sampling(
    Sgpc_cmpld, 
    data = SdataM1,
    iter = 30000,
    warmup = 5000,
    thin = 10,
    chains = 4
)

## PLOTS
library(coda)

munams <- paste0("mu[", 1:p, "]")
sxnams <- paste0("sigmas[", 1:p, "]")

thnams0 <- paste0("grpc_theta[", 1:dim(g0)[2], "]")
thnams1 <- paste0("grpc_theta[", 1:dim(g1)[2], "]")

rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

library(ggpubr)

mcmc_intervals(Samples1, paste0("x_sigmas[", 1:p, "]"))

## CI for the graphpcor model parameters (all edges are 'significant')
ggarrange(mcmc_intervals(Samples0, thnams0) +
         geom_abline(slope = 0, intercept = 0),
         mcmc_intervals(Samples1, thnams0) +
         geom_abline(slope = 0, intercept = 0))


ycorr <- cc[1:6, 1:6]
iil <- which(lower.tri(ycorr))

## CI for the correlation pairs (observed in 'red')
## g0 not enough
ggarrange(
    mcmc_intervals(Samples0, rhonams) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                 data.frame(
                     x = ycorr[iil],
                     ya = length(iil):1-0.1,
                     yb = length(iil):1+0.5), color = 'red'),
    mcmc_intervals(Samples1, rhonams) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                 data.frame(
                     x = ycorr[iil],
                     ya = length(iil):1-0.1,
                     yb = length(iil):1+0.5), color = 'red')
    )

library(ggplot2)


### the mu and sigmas from both graphpcor models
ggmu0 <- mcmc_areas(
  Samples0,
  pars = munams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
ggmu1 <- mcmc_areas(
  Samples1,
  pars = munams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
ggsd0 <- mcmc_areas(
  Samples0,
  pars = sxnams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
ggsd1 <- mcmc_areas(
  Samples1,
  pars = sxnams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

ggm_add <- geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                        data.frame(
                            x = colMeans(Sdata0$y),
                            ya = ncol(Sdata0$y):1-0.1,
                            yb = ncol(Sdata0$y):1+0.5), color = 'red')
ggsd_add <-  geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                          data.frame(
                              x = apply(Sdata0$y, 2, sd),
                              ya = ncol(Sdata0$y):1-0.1,
                              yb = ncol(Sdata0$y):1+0.5), color = 'red') 


ggarrange(ggmu0 + ggm_add,
          ggmu1 + ggm_add,
          ggsd0 + ggsd_add,
          ggsd1 + ggsd_add)

 
ggcc0 <- mcmc_areas(
  Samples0,
  pars = rhonams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

ggcc1 <- mcmc_areas(
  Samples1,
  pars = rhonams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

ggcc_add <- geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                         data.frame(
                             x = ycorr[iil],
                             ya = length(iil):1-0.1,
                             yb = length(iil):1+0.5), color = 'red')

## it is clear that g1 is better and enough
ggarrange(ggcc0 + ggcc_add,
          ggcc1 + ggcc_add)
