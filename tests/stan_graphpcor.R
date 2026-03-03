## EXAMPLE definition
## STAN model definition
##   initial code
##   update code
##   compile code
## DATA sample
##   CORRELATION model definition
##   PARAMETER definition
##   SAMPLE drawn
##   initial STAN data
##   base model definition
##   update STAN data with prior parameters
## MCMC sampling
## PLOTS

## EXAMPLE definition
## y ~ N(\mu, \Sigma)
##  where \mu = [0,...,0] and \Sigma = C(\theta) 
## theta[1:m]: vector of parameters to be placed in L0[iil]
##       where iil is an index vector for the parameters in L0,
##       L0 the initial Cholesky, then compute its fill-in elements,
##       as it is a Cholesky of initial precision matrix Q0.
##       Then compute V0 and then C(\theta)

## MODEL definition
##   initial STAN code
grpc0 <- '
data {
  int<lower=1> n;
  int<lower=1> p;
  vector[p] y[n];
  vector[p] ymu;
}
transformed parameters {
  matrix[p,p] rho;
}
model {
  y ~ multi_normal(ymu, rho);
}
'

##   update STAN code 
library(graphpcor)

## add the code for the PC prior for a correlation matrix,
## which was named "rho" in the intial model code,
## considering a generic graph (not yet defined!!!)
grpc_code <- stan_add(grpc0, 'graphpcor', lambda = 1, "rho")

##   compile STAN code
library(rstan)
options(mc.cores = 4L)

system.time(
    pcgmodel <- stan_model(
        model_code = grpc_code, 
        model_name = "pc_graphpcor"
    )
)

## DATA sample
## Correlation model definition
g <- graphpcor(x1~x2+x3, x4~x2+x3)
p <- dim(g)[1]
m <- dim(g)[2]
c(p, m)

##   PARAMETER definition
set.seed(2)
th.true <- rnorm(m)
corr <- basepcor(th.true, p, g)$base
Uc <- chol(corr)

##   SAMPLE drawn
n <- 100
set.seed(3)
y <- matrix(rnorm(n*p), n) %*% Uc

(ycorr <- cor(y))

##   initial STAN data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(p),
    y = y,
    ymu = rep(0,p)
)

##   base model definition
## prior parameters definition
set.seed(1)
th.base <- rnorm(m)*0
baseC <- basepcor(th.base, p = p, iLtheta = g)

##   update STAN data
Sdata <- stan_add(Sdata0, baseC, lambda = 1, name = 'rho')

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

thnams <- paste0("grpc_theta[", 1:m, "]")
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
            if(length(Sdata$grpc_ii)>k1 &&
               (i==Sdata$grpc_jj[k1+1]) &
               (j==Sdata$grpc_ii[k1+1])) {
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

