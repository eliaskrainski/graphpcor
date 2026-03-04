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

## model idea for other than famhist and chd variables:
##     y[i,j] = Normal(mu[j] + x_sigmas[j] * x[i,j], y_sigmas[j])
##   model the latent x as
##     x ~ multi_normal(0, R), R a correlation matrix

## STAN model base code
Scode0 <- '
data {
  int<lower=1> n;
  int<lower=1> p;
  vector[p] y[n];
  real<lower=0> mu_sigma;
  real<lower=0> x_sigmas_lambda;
}
parameters {
  vector[p] mu;
  vector<lower=0>[p] x_sigmas;
}
transformed parameters {
  matrix[p,p] rho;
}
model {
  mu ~ normal(0, mu_sigma);
  x_sigmas ~ exponential(x_sigmas_lambda);
  y ~ multi_normal(mu, quad_form_diag(rho, x_sigmas));
}
'

## update the STAN code with code for the
## graphpcor prior for 'rho'
Sgrpc <- stan_add(
    Scode0, 'graphpcor',
    lambda = 1, name = "rho")

cat(Sgrpc)

## compile STAN code
library(rstan)
options(mc.cores = 4L)

system.time(
    Sgpc_cmpld <- stan_model(
        model_code = Sgrpc,
        model_name = "graphpcor"
    )
)

## STEP 2: prepare the data
data(SAheart, package = "msos")

str(SAheart)

## do not consider in this example the
## two discrete variables 'chd' and 'famhist'
dataf <- data.frame(
    SAheart[, c(1:4, 6:9)]
)
str(dataf)

(n <- nrow(dataf))
(p <- ncol(dataf))

## STAN (initial) data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(p),
    y = as.matrix(dataf),
    mu_sigma = 10,
    x_sigmas_lambda = 1,
    y_sigmas_lambda = 1
)

str(Sdata0)

## to update the data we need to define the graphpcor model

## observed correlation matrix
round((cc <- cor(dataf)) * 100)

vnams <- colnames(dataf)

lcc <- chol(cc)
qc <- chol2inv(lcc)

## observed partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(vnams, vnams)
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

## Define the graph from the MST
G0 <- matrix(0, p, p, dimnames = dimnames(pC))
for(i in 1:nrow(mst)) {
    G0[mst[i,1], mst[i,2]] <- 1
    G0[mst[i,2], mst[i,1]] <- 1
}
G0

## graphcor from the MST graph
g0 <- graphpcor(G0)
(dg0 <- dim(g0))

p*(p-1)/2

## define a graph consideing partial correlations (abs) > 0.09
tanh(qnorm(c(0.025, 0.975)) / sqrt(n-3))

## graphpcor from the (observed) partial correlations
g1 <- graphpcor(abs(pC)>0.05) ## actually >0.05
G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

## visualize both graphs
if(FALSE)
    png("g01SAheart.png", width = 1800, height = 900, res = 300)
par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(g0, Rgraphviz = TRUE)
plot(g1, Rgraphviz = TRUE)
if(FALSE)
    dev.off()

if(FALSE)
    system("eog g01SAheart.png &")

## base model for each graph
baseM0 <- basepcor(diag(p), iLtheta = g0)
baseM1 <- basepcor(diag(p), iLtheta = g1)

##   update STAN data for each graph
SdataM0 <- stan_add(Sdata0, baseM0, lambda = 1, name = 'rho')
SdataM1 <- stan_add(Sdata0, baseM1, lambda = 1, name = 'rho')

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
sxnams <- paste0("x_sigmas[", 1:p, "]")

thnams0 <- paste0("grpc_theta[", 1:dim(g0)[2], "]")
thnams1 <- paste0("grpc_theta[", 1:dim(g1)[2], "]")

rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

library(ggpubr)

## CI for the graphpcor model parameters (all edges are 'significant')
ggarrange(mcmc_intervals(Samples0, thnams0) +
         geom_abline(slope = 0, intercept = 0),
         mcmc_intervals(Samples1, thnams0) +
         geom_abline(slope = 0, intercept = 0))


ycorr <- cc
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
