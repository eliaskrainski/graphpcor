## EXAMPLE definition
## MODEL definition
## DATA sample
##   PARAMETER definition
##   SAMPLE drawn
## MCMC sampling
## PLOTS

## EXAMPLE definition
## y ~ N(\mu, \Sigma)
##  where \mu = [0,...,0] and \Sigma = C(\theta) 
## theta[1:m]: vector of parameters to be placed in L[iil]
##       where iil is an index vector for the lower triangle

## MODEL definition
pccor <- '
data {
  int n;
  int p;
  int m;
  vector[p] y[n];
  vector[p] ymu;
  real<lower=0> lambda;
  real<lower=0> logDetIhalf;
  vector[m] thetaBase;
  matrix[m,m] Ihalf;
}
parameters {
  vector[m] theta;     // correlation parameters
}
transformed parameters {
  vector[m] xi = Ihalf * (theta - thetaBase);
  real<lower=0> r = dot_self(xi);
  matrix[p,p] A;
  matrix[p,p] B;
  cholesky_factor_corr[p] L;
  L[1,1] = 1.0;
  for(i in 1:p) {
    A[i,i] = 1.0;
    B[i,i] = 1.0;
  }
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
        A[i,j] = 0.0;
        B[i,j] = 0.0;
        L[i,j] = 0.0;
    }
  }
  {
  int k = 0;
  for(j in 1:(p-1)) {
    for(i in (j+1):p) {
      k = k+1;
      A[i,j] = tanh(theta[k]);
      B[i,j] = sqrt(1-A[i,j]^2);
    }
  }
  }
  for(j in 2:(p-1)) {
    for(i in j:p) {
      B[i,j] *= B[i,j-1];
    }
  }
  for(i in 1:p) {
     L[i, 1] = A[i,1];
  }
  for(j in 2:p) {
    for(i in j:p) {
      L[i,j] = A[i,j]*B[i,j-1];
    }
  }
  corr_matrix[p] rho = tcrossprod(L);
}
model {
  y ~ multi_normal_cholesky(ymu, L);
  r ~ exponential(lambda);
  target += lgamma(m * 0.5) - m*0.5 * log(pi());
  target += logDetIhalf -(m-1.0)*log(r);
}
'

## strsplit(smodel, "\n", fixed = TRUE)[[1]][35]

## MCMC sampling
library(rstan)
options(mc.cores = 4L)

system.time(
    Spccor <- stan_model(
        model_code = pccor, 
        model_name = "pc_dense"
    )
)

## DATA sample
##   PARAMETER definition
corr <- matrix(c( 1.0, 0.7, 0.5,-0.3,
                  0.7, 1.0,-0.2, 0.4,
                  0.5,-0.2, 1.0,-0.8,
                 -0.3, 0.4,-0.8, 1.0), 4)
all.equal(corr,t(corr))
Uc <- chol(corr)

p <- ncol(corr)

## DATA sample
n <- 100
set.seed(1)
y <- matrix(rnorm(n*p), n) %*% Uc

corr
(ycorr <- cor(y))

library(graphpcor)
th.obs <- basecor(ycorr)$theta
(th.true <- basecor(corr)$theta)

## parameters for the PC-prior
m <- p * (p-1)/2
set.seed(2)
th.base <- rnorm(m)
baseC <- basecor(th.base, p = p)
I0 <- hessian(baseC)
I0dec <- graphpcor:::dspd(I0)

## STAN data
Sdata <- list(
    n = as.integer(n),
    p = as.integer(p),
    m = as.integer(m), 
    y = y,
    ymu = rep(0,p),
    lambda = 1,
    logDetIhalf = 0.5 * I0dec$logDeterminant,
    thetaBase = baseC$theta,
    Ihalf = I0dec$sqrt
)

## STAN sampling
Samples <- sampling(
    Spccor,
    data = Sdata,
    iter = 30000,
    warmup = 5000,
    chains = 4
)

## PLOTS
library(coda)

thnams <- paste0("theta[", 1:m, "]")

rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

th.true
th.base

mcmc_intervals(Samples, thnams)

mcmc_intervals(Samples, rhonams)

library(ggplot2)

ggcc <- mcmc_areas(
  Samples,
  pars = rhonams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = corr[lower.tri(corr)],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'red', lty = 2) +
 geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = ycorr[lower.tri(corr)],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'black', lty = 2) +
 geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = baseC$base[lower.tri(corr)],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'blue', lty = 2)

ggcc

## my own plots
Sth <- Reduce("cbind", extract(Samples, thnams))
colnames(Sth) <- thnams
Scorr <- Reduce("cbind", extract(Samples, rhonams))
colnames(Scorr) <- rhonams

summary(Scorr)
postCmean <- colMeans(Scorr)
postCmean

iil <- which(lower.tri(corr))

par(mfrow = c(p, p), mar = c(3.5,1.5,0.5,0.5),
    mgp = c(1.5,0.5,0), bty = "n")
k2 <- k1 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            if(j==1)
                legend("topleft", c("Observed", "TRUE", "base", "Posterior"),
                       bty = "n", col = c(1, 2, 4, 0),
                       lty = 2, lwd = 2, fill = c(0,0,0,gray(0.5)),
                       border = 'transparent')
        }
        if(j>i) {
            k1 <- k1 + 1
            thk1 <- c(th.obs[k1], th.true[k1], th.base[k1])
            h <- hist(Sth[, k1], 100, plot = FALSE)
            plot(h, main = '', freq = FALSE, border = 'transparent',
                 xlim = range(h$breaks, thk1),
                 xlab = as.expression(bquote(theta[.(k1)])))
            abline(v = thk1, col = c(1,2,4),
                   lty = 2, lwd = 2)
        }
        if(j<i) {
            k2 <- k2 + 1
            cc <- c(ycorr[iil[k2]], corr[iil[k2]],
                    baseC$base[iil[k2]])
            hS <- hist(Scorr[,k2], 100, plot = FALSE)
            plot(hS, freq = FALSE, main = '', border = 'transparent',
                 xlab = as.expression(bquote(rho[.(i)~","~.(j)])),
                 xlim = range(hS$breaks, cc, na.rm = TRUE))
            abline(v = cc, col = c(1,2,4), lty = 2, lwd = 2)
        }
    }
}

