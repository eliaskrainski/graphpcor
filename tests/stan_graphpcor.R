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
pcgraphpcor <- '
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
  vector[p] d0;
  array[m] int ii;
  array[m] int jj;
  int nfi;
  array[nfi] int ifi;
  array[nfi] int jfi;
}
parameters {
  vector[m] theta;     // correlation parameters
}
transformed parameters {
  vector[m] xi = Ihalf * (theta - thetaBase);
  real<lower=0> r = dot_self(xi);
  matrix[p,p] L;
  for(i in 1:p) {
    for(j in 1:p) {
      if(i==j) {
        L[i,i] = d0[i];
      } else {
        L[i,j] = 0.0;
      } 
    }
  }
  {int i,j;
  for(k in 1:m) {
    i = ii[k];
    j = jj[k];
    L[i,j] = theta[k];
  }
  if(nfi>0) {
    for(l in 1:nfi) {
      i = ifi[l];
      j = jfi[l];
      if(j>0) {
        for(k in 1:(j-1)) {
          L[i,j] -= L[i,k] * L[j,k];
        }
      }
    }
  }
  }
  matrix[p,p] V = chol2inv(L);
  vector[p] s0inv = inv_sqrt(diagonal(V));
  matrix[p,p] rho = quad_form_diag(V,s0inv);
}
model {
  y ~ multi_normal(ymu, rho);
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
    pcgmodel <- stan_model(
        model_code = pcgraphpcor,
        model_name = "pc_graphpcor"
    )
)

## Correlation model definition
library(graphpcor)
g <- graphpcor(x1~x2+x3, x4~x2+x3)
p <- dim(g)[1]
m <- dim(g)[2]
c(p, m)

## prior parameters definition
set.seed(1)
th.base <- rnorm(m)*0
baseC <- basepcor(th.base, p = p, iLtheta = g)

I0 <- hessian(baseC)
I0
I0dec <- graphpcor:::dspd(I0)
str(I0dec)

## DATA sample
##   PARAMETER definition
set.seed(2)
th.true <- rnorm(m)
corr <- basepcor(th.true, p, g)$base
Uc <- chol(corr)

## DATA sample
n <- 100
set.seed(3)
y <- matrix(rnorm(n*p), n) %*% Uc

(ycorr <- cor(y))

## prepare the prior definiton
## index for L[theta]
Q0 <- Laplacian(g)
iLth <- which((Q0 !=0) & (lower.tri(Q0)))
Q0
iLth

iijj <- list(
    ii = row(Q0)[iLth],
    jj = col(Q0)[iLth])

## fill-in index
iLfill <- setdiff(which(t(chol(Q0 + Diagonal(p))) != 0),
                  which(Q0 != 0))
iLfill

iijj$nfi <- length(iLfill)
if(iijj$nfi>0) {
    iijj$ifi <- row(Q0)[iLfill]
    iijj$jfi <- col(Q0)[iLfill]
}
iijj

## STAN data
Sdata <- list(
    n = as.integer(n),
    p = as.integer(p),
    m = as.integer(m), 
    y = y,
    ymu = rep(0,p),
    lambda = 1,
    logDetIhalf = 0.5 * abs(I0dec$logDeterminant),
    thetaBase = baseC$theta,
    Ihalf = I0dec$sqrt,
    d0 = baseC$d0,
    ii = as.array(iijj$ii),
    jj = as.array(iijj$jj),
    nfi = iijj$nfi,
    ifi = as.array(iijj$ifi),
    jfi = as.array(iijj$jfi)
)

## STAN sampling
Samples <- sampling(
    pcgmodel,
    data = Sdata,
    iter = 30000,
    warmup = 5000,
    chains = 4
)

## PLOTS
library(coda)

##str(Samples)

thnams <- paste0("theta[", 1:m, "]")
rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

th.base
th.true
mcmc_intervals(Samples, thnams)

baseC$base
corr
ycorr
mcmc_intervals(Samples, rhonams)

library(ggplot2)

iil <- which(lower.tri(corr))

ggcc <- mcmc_areas(
  Samples,
  pars = rhonams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = ycorr[iil],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'black', lty = 2) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = corr[iil],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'red', lty = 2) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = baseC$base[iil],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'blue', lty = 2) 

ggcc


## my own plots
Sth <- Reduce("cbind", extract(Samples, thnams))
colnames(Sth) <- thnams
Scorr <- Reduce("cbind", extract(Samples, rhonams))
colnames(Scorr) <- rhonams

summary(Scorr)
(ycorr <- cor(y))
postCmean <- colMeans(Scorr)
postCmean


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
            if(length(iijj$ii)>k1 && (i==iijj$jj[k1+1]) & (j==iijj$ii[k1+1])) {
                k1 <- k1 + 1
                thk1 <- c(NA, th.true[k1], th.base[k1])
                h <- hist(Sth[, k1], 100, plot = FALSE)
                plot(h, main = '', freq = FALSE, border = 'transparent',
                     xlim = range(h$breaks, thk1, na.rm = TRUE),
                     xlab = as.expression(bquote(theta[.(k1)])))
                abline(v = thk1, col = c(1,2,4),
                       lty = 2, lwd = 2)                
            } else {
                plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            }
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

if(FALSE) {  ## some checks at first drawn
    
    s1 <- sapply(Samples@sim$samples, sapply, head, 1)
    s1
    
    baseC$theta
    th.true
    s1[1:4, ]
    
    graphpcor:::Lprec0(s1[1:4, 1], p, iLtheta = g, d0 = baseC$d0)
    matrix(s1[p*2+1+1:(p^2), 1], p)
    
    matrix(s1[p*2 + 1 + p^2 + 1:(p^2), 1], p)
    chol2inv(t(matrix(s1[p*2+1 + 1:(p^2), 1], p)))
    
    sqrt(1/diag(matrix(s1[p*2+1 + p^2 + 1:(p^2), 1], p)))
    s1[p*2 + 1 + (p^2)*2 + 1:p, 1]
    
    matrix(s1[p*2+1+(p^2)*2+p + 1:(p^2), 1], p)
    basepcor(s1[1:4, 1], p, g)$base

    s1[p+1:p, ]
    Sdata$Ihalf %*% s1[1:p, ]

    colSums(s1[p+1:4,]^2)
    s1[p*2+1, ]

}
