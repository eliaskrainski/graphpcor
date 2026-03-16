setwd(here::here("examples"))
getwd()

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

## model idea:
##   Model chd (the response) is modeled as 
##     chd[i] ~ bernoulli(inv_logit(alpha_1 + beta[j] * x[i,j]))
##       for i = 1, ..., n and j = 1, ..., p 
##   where x is a latent variable modeled from the p covariates
##     x[i,] ~ N(0,R), R is a p = 9 dimensional correlation matrix
##   Model the first p-3 contiuous variables as 
##     y[i,j] = x[i,j], j = 1, ..., p-3
##   Model the next two (p-2 to p-1) contiuous variables as a Gamma(a, b=1)
##     y[i,j] ~ G1(a=exp(x[i,j]), j = p-2, p-1
##   Model famhist, coded as 0 (if 'Absent') or 1 (if 'Present') as
##     y[i,p] ~ bernoulli(inv_logit(alpha_2 + alpha_3 * x[i,p]))
## NOTE 1: y[,j], for j in 1, ..., p-2 were standardized
##         y[i,j] = x[i,j] for j in 1, ..., p-3
##         log(E(y[i,j])) = x[i,j] for j in p-2 and p-1
##         E(y[i,p]) = p_{ij} = 1/(1+exp(-alpha_2 - alpha_3 * x[i,p]))
## Note 2: I am not sure if we need alpha_3
##   (is there because famHist is the only covariate that is not scaled)
##   (fix it to one may be enough, its posterior is just above 1)
## Note 3: extension to allow error for the first p-1 covariates

## STAN model base code
## supply yx as the first p-1 covariates
Scode0 <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  array[n] int<lower=0,upper=1> chd;
  matrix[n,p-1] yx;
  array[n] int<lower=0,upper=1> famHist;
  real<lower=0> alphas_sigma;
  real<lower=0> betas_sigma;
}
parameters {
  vector[2] alphas;
  vector[p] betas;
  matrix[n,p] z; // latent for E() of the covariates
}
transformed parameters {
  corr_matrix[p] rho;
}
model {
  alphas ~ normal(0.0, alphas_sigma);
  betas ~ normal(0.0, betas_sigma);
  to_vector(z) ~ std_normal();
  matrix[p,p] L_rho;
  L_rho = cholesky_decompose(rho);
  matrix[n,p] x = z * (L_rho');
  to_vector(yx[,1:(p-3)]) ~ normal(to_vector(x[, 1:(p-3)]), 0.001);
  to_vector(yx[,(p-2):(p-1)]) ~ gamma(to_vector(exp(x[,(p-2):(p-1)])), 1.0);
  famHist ~ bernoulli(Phi(alphas[2] + x[,p]));
  chd ~ bernoulli_logit(alphas[1] + x * betas);
}
"

cat(Scode0, file = "Scode0_v2.txt")

## STEP 2: update the STAN code with code for the
## graphpcor prior for 'rho'
Sgrpc <- stan_add(
    Scode0, 'graphpcor',
    lambda = 1, name = "rho")

##Sgrpc <- Scode0

cat(Sgrpc)

## STEP 3:  compile STAN code
library(rstan)
options(mc.cores = 4L)
rstan_options(auto_write = TRUE)

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
    yx = as.matrix(xdata[, 1:(p-1)]),
    famHist = (SAheart$famhist=="Present")+0L, 
    alphas_sigma = 10,
    betas_sigma = 10
)

str(Sdata0)

round((cc <- cor(xdata)) * 100)

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

g0 <- graphpcor(G0)
(dg0 <- dim(g0))

p*(p-1)/2

## define a graph consideing partial correlations (abs) > 0.09
tanh(qnorm(c(0.025, 0.975)) / sqrt(n-3))

g1 <- graphpcor(abs(pC)>0.07)
G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

if(FALSE)
    png("g01SAheart.png", width = 1800, height = 900, res = 300)
par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(g0, Rgraphviz = TRUE)
plot(g1, Rgraphviz = TRUE)
if(FALSE)
    dev.off()

if(FALSE)
    system("eog g01SAheart.png &")

## base model
baseM0 <- basepcor(diag(p), iLtheta = g0)
baseM1 <- basepcor(diag(p), iLtheta = g1)

##   update STAN data
SdataM0 <- stan_add(Sdata0, baseM0, lambda = 1, name = 'rho')
SdataM1 <- stan_add(Sdata0, baseM1, lambda = 1, name = 'rho')

## STAN sampling
## STAN sampling for each graph
Samples0 <- sampling(
    Sgpc_cmpld, 
    data = SdataM0,
    iter = 3000,
    warmup = 500,
    thin = 1,
    chains = 4
)

Samples1 <- sampling(
    Sgpc_cmpld, 
    data = SdataM1,
    iter = 3000,
    warmup = 500,
    thin = 1,
    chains = 4
)

chd0fit <- glm(
    SAheart$chd ~ as.matrix(xdata),
    family = 'binomial'
)
chd0betas <- coef(summary(chd0fit))[-1,]
chd0betas

alphas0 <- c(qnorm(sum(xdata$famHist)/n),
             coef(chd0fit)[1])
alphas0

## PLOTS
library(coda)

anames <- paste0("alphas[", 1:2, "]")
bnames <- paste0("betas[", 1:p, "]")

thnams0 <- paste0("grpc_theta[", 1:dim(g0)[2], "]")
thnams1 <- paste0("grpc_theta[", 1:dim(g1)[2], "]")

rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

mcmc_trace(Samples0, c(anames, bnames))
mcmc_trace(Samples1, c(anames, bnames))

mcmc_trace(Samples0, thnams0)
mcmc_trace(Samples1, thnams1)

library(ggpubr)

## CI for the alphas and betas 
ggarrange(mcmc_intervals(Samples0, c(anames, bnames)) +
          geom_abline(slope = 0, intercept = 0) +
          geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                       data.frame(
                           x = c(alphas0, chd0betas[, 1]),
                           ya = (2+p):1-0.1,
                           yb = (2+p):1+0.5), color = 'red'),
         mcmc_intervals(Samples1, c(anames, bnames)) +
         geom_abline(slope = 0, intercept = 0) +
          geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                       data.frame(
                           x = c(alphas0, chd0betas[, 1]),
                           ya = (2+p):1-0.1,
                           yb = (2+p):1+0.5), color = 'red'))

## CI for the graphpcor model parameters (all edges are 'significant')
ggarrange(mcmc_intervals(Samples0, thnams0) +
         geom_abline(slope = 0, intercept = 0),
         mcmc_intervals(Samples1, thnams0) +
         geom_abline(slope = 0, intercept = 0))


ycorr <- cc
ylcorr <- cor(cbind(xdata[, 1:(p-3)],
                    log(xdata[, (p-2):(p-1)]), xdata[, p, drop = FALSE]))
iil <- which(lower.tri(ycorr))

## CI for the correlation pairs (observed in 'red')
## g0 not enough
ggarrange(
    mcmc_intervals(Samples0, rhonams) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb, color = ctype),
                 data.frame(
                     ctype = rep(c("y", "log(y)"), each = length(iil)),
                     x = c(ycorr[iil], ylcorr[iil]),
                     ya = length(iil):1-0.1,
                     yb = length(iil):1+0.5)), 
    mcmc_intervals(Samples1, rhonams) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb, color = ctype),
                 data.frame(
                     ctype = rep(c("y", "log(y)"), each = length(iil)),
                     x = c(ycorr[iil], ylcorr[iil]),
                     ya = length(iil):1-0.1,
                     yb = length(iil):1+0.5))
    )

## the latents
par(mfrow = c(3, 3), mar = c(4,4,1,1), mgp = c(2,0.5,0), bty = "n")
for(j in 1:(p-3)) {
    xnm <- paste0("z[", 1:n, ",", j, "]")
    xxsamples <- do.call('cbind', extract(Samples1, xnm))
    plot(colMeans(xxsamples), xdata[, j], 
         xlab = as.expression(bquote(x[.(j)])), ylab = xnames[j])
}
for(j in (p-2):(p-1)) {
    xnm <- paste0("z[", 1:n, ",", j, "]")
    xxsamples <- do.call('cbind', extract(Samples1, xnm))
    plot(xdata[, j]+0.01, exp(colMeans(xxsamples)), log = 'xy',
    xlab = as.expression(bquote(x[.(j)])), ylab = xnames[j])
}
j = p
    xnm <- paste0("z[", 1:n, ",", j, "]")
    xxsamples <- do.call('cbind', extract(Samples1, xnm))
plot(colMeans(xxsamples) ~ SAheart$famhist,
     xlab = "famHist",
     ylab = as.expression(bquote(x[.(j)])))
